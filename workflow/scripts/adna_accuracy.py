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
            # self.correct,
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
            overmapped_bact=100 * self.overmapped_bact / self.origin_count["b"],
            overmapped_cont=100 * self.overmapped_cont / self.origin_count["c"],
            jaccard_correct=round(100 * self.jaccard_correct / self.n, 5),
            score_correct=100 * self.score_correct / self.n if self.score_correct is not None else "",
        )
    

def parse_gargammel_name(read_name: str):
    match = re.match(r'(.*):(\+|\-):(\d+):(\d+):(\d+)(a|b|c|d|e)(\d+)(.*)', read_name)
    if not match:
        raise ValueError(f"Could not parse read name: {read_name}, please make sure that the reads were simulated using gargammel")

    ref_name, orientation, start_pos, end_pos, read_len, read_origin, allele, _ = match.groups()
    return ref_name, start_pos, end_pos, read_origin


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


def filter_bam(alignment_file):
    for record in alignment_file:
        if not record.is_supplementary and not record.is_secondary:
            yield record
            

def get_iter_stats(fastq_path, predicted) -> Accuracy:
    n = 0
    unaligned = 0
    nr_aligned = 0
    overmapped_bact = 0
    overmapped_cont = 0
    correct = 0
    correct_jaccard = 0.0
    correct_score = 0  # Same or better alignment score

    score_threshold = 20
    unscored = 0

    origin_count = {"b": 0, "c": 0, "e": 0}
    with xopen(fastq_path) as fastq_handle:
        for read in SeqIO.parse(fastq_handle, "fastq"):
            _, _, _, read_origin = parse_gargammel_name(read.name)
            origin_count[read_origin] += 2
    # iterator = zip_longest(SeqIO.parse(truth, "fastq"), filter_bam(predicted))
    # print(origin_count)

    for p in filter_bam(predicted):

        ref_name, start_pos, end_pos, read_origin = parse_gargammel_name(p.query_name)

        n += 1
        try:
            score = p.get_tag("AS") if p.has_tag("AS") else recompute_alignment_score(p, Scores)
        except:
            unscored += 1
            score = score_threshold
        if not p.is_unmapped and score >= score_threshold:  # TODO and not p.is_unmapped:
            # print(p.query_name)
            if read_origin == "b":
                overmapped_bact += 1
            elif read_origin == "c":
                overmapped_cont += 1
            else:
                nr_aligned += 1
        # if t.is_unmapped and p.is_unmapped:
        #     continue

        is_correct = False
        if ref_name == p.reference_name:
            # TODO obtain true coordinates for variated genomes
            correct += 1
            correct_jaccard += 1
        #     jacc = jaccard_overlap(
        #         p.reference_start, p.reference_end, t.reference_start, t.reference_end
        #     )
        #     correct_jaccard += jacc
        #     if overlap(
        #         p.reference_start, p.reference_end, t.reference_start, t.reference_end
        #     ):
        #         correct += 1
        #         is_correct = True

        # if not is_correct:
        #     if recompute_predicted_score:
        #         predicted_score = recompute_alignment_score(p, Scores)
        #     else:
        #         predicted_score = p.get_tag("AS")
        #     truth_score = recompute_alignment_score(t, Scores)
        #     # print(f"true: {t.reference_name} {t.reference_start} {t.cigarstring} AS:{truth_score}  -- actual: {p.reference_name} {p.reference_start} {p.cigarstring} AS:{predicted_score}")
        #     correct_score += predicted_score >= truth_score

    # print(f"{unscored} alignments without a score")
    return Accuracy(
        n=n,
        aligned=nr_aligned,
        correct=correct,
        jaccard_correct=correct_jaccard,
        score_correct=(correct_score + correct),
        overmapped_bact=overmapped_bact,
        overmapped_cont=overmapped_cont,
        origin_count=origin_count
    )

def measure_accuracy(
    fastq: Path,
    predicted: Path
    # outfile: Path = None,
    # skip_r2: bool = False,
    # recompute_score: bool = False,
    # multiple_primary: bool = False,
) -> Accuracy:
    with (AlignmentFile(predicted) as predicted):
        # if skip_r2:
        #     truth = skip_r2_iter(truth)
        # if multiple_primary:
        #     if skip_r2:
        #         predicted = pick_random_primary_single_end_iter(predicted)
        #     else:
        #         predicted = pick_random_primary_paired_end_iter(predicted)
        result = get_iter_stats(fastq, predicted)

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
    parser.add_argument("--fastq", type=Path, help="Path to gargammel simulated reads")
    parser.add_argument("--predicted", "--predicted_sam", "--predicted_paf", type=Path, help="Predicted SAM/BAM/PAF")
    # parser.add_argument("--paf", dest="force_paf", action="store_true", help="Assume PAF input for predicted (usually autodetected, only needed if reading PAF from stdin)")
    parser.add_argument("--outfile", help="Path to file")
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    accuracy = measure_accuracy(args.fastq, args.predicted)

    # print("\% Mapped endogenous reads\t\% Mapped bacterial reads\t\% Mapped contaminated reads")
    print(*accuracy.percentages().row(), sep="\t")