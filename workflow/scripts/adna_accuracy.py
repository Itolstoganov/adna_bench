#!/usr/bin/env python3

import argparse
import gzip
import os
import sys
from collections import defaultdict
from dataclasses import dataclass
from itertools import zip_longest, groupby
from pathlib import Path
import re

from Bio import SeqIO
from pysam import AlignmentFile
from typing import Dict, List, Set, Tuple
from xopen import xopen

# Parsing functions from https://github.com/NBISweden/strobealign-evaluation

READ_ORIGIN_DICT = {'b': "bact", 'c': "cont", 'e': "endo"}

CIGAR_MATCH = 0  # M
CIGAR_INSERTION = 1  # I
CIGAR_DELETION = 2  # D
CIGAR_SOFTCLIP = 4  # S

class Scores:
    match = 2
    mismatch = 8
    gap_open = 12
    gap_extension = 1
    end_bonus = 10

@dataclass
class ReferenceInterval:
    name: str
    origin: str
    start: int
    end: int

    def __iter__(self):
        return iter((self.name, self.origin, self.start, self.end))


@dataclass
class AccuracyPercentages:
    aligned: float
    correct: float
    score_correct: float
    jaccard_correct: float
    overmapped_bact: int
    overmapped_cont: int

    def row(self):
        return (
            self.aligned,
            self.correct,
            self.overmapped_bact,
            self.overmapped_cont
            # self.unaligned,
            # self.jaccard_correct,
            # self.score_correct if self.score_correct is not None else "",
        )


@dataclass
class Accuracy:
    n: int
    aligned: int
    correct: int
    score_correct: int | None  # unavailable for PAF
    jaccard_correct: int
    overmapped_bact: int
    overmapped_cont: int
    origin_count: dict

    def percentages(self) -> AccuracyPercentages:
        return AccuracyPercentages(
            aligned=100 * self.aligned / self.origin_count["e"],
            correct=100 * self.correct / self.origin_count["e"],
            overmapped_bact=100 * self.overmapped_bact / self.origin_count["b"] if self.origin_count["b"] > 0 else 0,
            overmapped_cont=100 * self.overmapped_cont / self.origin_count["c"] if self.origin_count["c"] > 0 else 0,
            jaccard_correct=round(100 * self.jaccard_correct / self.n, 5),
            score_correct=100 * self.score_correct / self.n if self.score_correct is not None else "",
        )


@dataclass
class AccuracyResults:
    thresholds: List[int]
    results: Dict[int, Accuracy]

    def rows(self):
        for threshold in self.thresholds:
            yield (threshold, *self.results[threshold].percentages().row())

    def to_tsv(self):
        lines = []
        for threshold in self.thresholds:
            row = [str(threshold)] + [str(x) for x in self.results[threshold].percentages().row()]
            lines.append('\t'.join(row))
        return '\n'.join(lines)
    

def parse_gargammel_name(read_name: str):
    ref_name, orientation, start_pos, end_pos, origin_string = read_name.split(":")
    match = re.match(r'(\d+)(a|b|c|d|e)(.*)', origin_string)

    if not match:
        print(origin_string)
        raise ValueError(f"Could not parse read name: {read_name}, please make sure that the reads were simulated using gargammel")
    
    read_len, read_origin, allele = match.groups()
    return ref_name, int(start_pos), int(end_pos), read_origin


def read_alignments(bam_path):
    bam_path = AlignmentFile(bam_path, check_sq=False)
    read_positions = {}

    bam_iter = bam_path.fetch(until_eof=True)
    for read in bam_iter:
        if read.is_secondary:
            continue
        if read.is_supplementary:
            continue
        if read.flag == 0 or read.flag == 16:  # single end
            # print('Single end')
            query_name = read.query_name
            if not query_name.endswith("/1"):
                query_name += "/1"
            if read.reference_end is None:
                read_positions[query_name] = False
            else:
                read_positions[query_name] = ReferenceInterval(
                    read.reference_name,
                    "endo",
                    read.reference_start,
                    read.reference_end,
                )
        elif read.is_paired:
            query_name = read.query_name
            if read.is_read1 and not query_name.endswith("/1"):
                query_name += "/1"
            if read.is_read2 and not query_name.endswith("/2"):
                query_name += "/2"

            if read.is_unmapped:
                read_positions[query_name] = False
            else:
                read_positions[query_name] = ReferenceInterval(
                    read.reference_name,
                    "endo",
                    read.reference_start,
                    read.reference_end,
                )
        elif read.is_unmapped:  # single and unmapped
            assert not read.is_paired
            query_name = read.query_name
            if not query_name.endswith("/1"):
                query_name += "/1"
            read_positions[query_name] = False

    return read_positions


def parse_int_or_not(s):
    try:
        return int(s)
    except ValueError:
        return s


def recompute_alignment_score(segment, scores) -> int:
    md_tag = segment.get_tag("MD")
    md = [
        parse_int_or_not(e)
        for e in re.split("([A-Z]|^[A-Z]+|[0-9]+)", md_tag)
        if e != ""
    ]
    score = 0
    # To compute the score, it is not necessary to reconstruct the full
    # alignment.
    # - The numbers in the MD tag give us the number of matches
    # - The letters in the MD tag give us the number of mismatches
    for item in md:
        if isinstance(item, int):
            # Match
            score += item * scores.match
        elif not item.startswith("^"):
            # Mismatch
            assert len(item) == 1
            score -= scores.mismatch

    # - I or D operations in the CIGAR string give us indels and their lengths
    cigartuples = segment.cigartuples
    for op, length in cigartuples:
        if op == CIGAR_INSERTION or op == CIGAR_DELETION:
            score -= scores.gap_open + (length - 1) * scores.gap_extension

    if cigartuples:
        if cigartuples[0][0] != CIGAR_SOFTCLIP:
            score += scores.end_bonus
        if cigartuples[-1][0] != CIGAR_SOFTCLIP:
            score += scores.end_bonus
    return score


def overlap(q_a, q_b, p_a, p_b):
    assert q_a <= q_b and p_a <= p_b
    # if (q_a == q_b) or (p_a == p_b):
    #     print("Cigar bug")
    return (
        (p_a <= q_a <= p_b)
        or (p_a <= q_b <= p_b)
        or (q_a <= p_a <= q_b)
        or (q_a <= p_b <= q_b)
    )


def jaccard_overlap(a_start, a_end, b_start, b_end):
    assert a_start < a_end
    assert b_start < b_end

    intersect = min(a_end, b_end) - max(a_start, b_start)
    if intersect < 0:
        return 0
    union = max(a_end, b_end) - min(a_start, b_start)
    result = intersect / union
    assert 0 <= result <= 1.0
    return result


assert jaccard_overlap(5, 10, 5, 10) == 1
assert jaccard_overlap(10, 20, 20, 30) == 0
assert jaccard_overlap(20, 30, 10, 20) == 0
assert jaccard_overlap(0, 4, 1, 3) == 0.5
assert jaccard_overlap(1, 3, 0, 4) == 0.5


def filter_bam(alignment_file):
    for record in alignment_file:
        if not record.is_supplementary and not record.is_secondary:
            yield record
            

def get_iter_stats(gt_path, predicted, score_thresholds) -> AccuracyResults:
    """
    Compute accuracy statistics for multiple score thresholds in a single pass.

    Args:
        gt_path: Path to ground truth BED file
        predicted: BAM file iterator
        score_thresholds: List of score thresholds to evaluate

    Returns:
        AccuracyResults object with stats for each threshold
    """
    stats = {}
    for threshold in score_thresholds:
        stats[threshold] = {
            'n': 0,
            'nr_aligned': 0,
            'overmapped_bact': 0,
            'overmapped_cont': 0,
            'correct': 0,
            'correct_jaccard': 0.0,
            'correct_score': 0
        }

    unscored = 0
    ground_truth = {}
    origin_count = {"b": 0, "c": 0, "e": 0}

    with open(gt_path) as gt_handle:
        for line in gt_handle:
            ref_name, start_pos, end_pos, read_name, _, orientation = line.strip().split()
            _, _, _, read_origin = parse_gargammel_name(read_name)
            origin_count[read_origin] += 1
            ground_truth[read_name] = ReferenceInterval(ref_name, read_origin, int(start_pos), int(end_pos))

    for p in filter_bam(predicted):
        if p.query_name not in ground_truth:
            print(f"{p.query_name} is not present in ground truth")
            continue

        truth = ground_truth[p.query_name]
        read_origin = truth.origin

        try:
            score = p.get_tag("AS") if p.has_tag("AS") else recompute_alignment_score(p, Scores)
        except:
            unscored += 1
            score = min(score_thresholds) - 1  

        jacc = 0.0
        is_overlap = False
        if not p.is_unmapped and truth.name == p.reference_name:
            try:
                jacc = jaccard_overlap(p.reference_start, p.reference_end, truth.start, truth.end)
                is_overlap = overlap(p.reference_start, p.reference_end, truth.start, truth.end)
            except:
                print(f"Malformed ground truth interval for read {p.query_name}", file=sys.stderr)
                print(p.reference_start, p.reference_end, file=sys.stderr)
                print(truth, file=sys.stderr)

        for threshold in score_thresholds:
            stats[threshold]['n'] += 1

            if not p.is_unmapped and score >= threshold:
                if read_origin == "b":
                    stats[threshold]['overmapped_bact'] += 1
                elif read_origin == "c":
                    stats[threshold]['overmapped_cont'] += 1
                elif read_origin == "e":
                    stats[threshold]['nr_aligned'] += 1

                    if is_overlap:
                        stats[threshold]['correct'] += 1
                    stats[threshold]['correct_jaccard'] += jacc
                
    results = {}
    for threshold in score_thresholds:
        s = stats[threshold]
        results[threshold] = Accuracy(
            n=s['n'],
            aligned=s['nr_aligned'],
            correct=s['correct'],
            jaccard_correct=s['correct_jaccard'],
            score_correct=(s['correct_score'] + s['correct']),
            overmapped_bact=s['overmapped_bact'],
            overmapped_cont=s['overmapped_cont'],
            origin_count=origin_count
        )

    return AccuracyResults(thresholds=score_thresholds, results=results)

def measure_accuracy(
    ground_truth: Path,
    predicted: Path,
    score_thresholds: List[int]
) -> AccuracyResults:
    with (AlignmentFile(predicted) as predicted):
        result = get_iter_stats(ground_truth, predicted, score_thresholds)

    return result


if __name__ == "__main__":
    # random.seed(0)
    parser = argparse.ArgumentParser(
        description="Calc identity",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--skip-r2", "--only-r1", default=False, action="store_true", help="Skip R2 reads in truth BAM")
    parser.add_argument("--recompute-score", default=False, action="store_true", help="Recompute score in *predicted* BAM. Default: Use score from AS tag")
    parser.add_argument("--multiple-primary", default=False, action="store_true", help="Allow multiple primary alignments (violates SAM specification) and pick one randomly")
    parser.add_argument("--synthesize-unmapped", default=False, action="store_true", help="If an alignment is missing from predicted, assume the read is unmapped")
    parser.add_argument("--ground-truth", type=Path, help="Path to ground truth alignments (in .bed format)")
    parser.add_argument("--predicted", "--predicted_sam", "--predicted_paf", type=Path, help="Predicted SAM/BAM/PAF")
    parser.add_argument("--score-thresholds", type=int, nargs='+', default=[25, 30, 35, 40, 45], help="Score thresholds to evaluate")
    parser.add_argument("--outfile", help="Path to file")
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    accuracy_results = measure_accuracy(args.ground_truth, args.predicted, args.score_thresholds)

    print(accuracy_results.to_tsv())