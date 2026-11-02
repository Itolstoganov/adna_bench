#!/usr/bin/env python
"""Find reads correct in BAM A but incorrect in BAM B, relative to ground truth.

Default (primary) correctness reuses the pipeline's authoritative definition from
adna_accuracy.py: a read is correct iff it is mapped, on the same reference as the
truth, its alignment interval overlaps the truth interval, and the read is
endogenous (origin "e").

With --secondary, correctness instead matches the "% Accuracy with secondary"
column of accuracy_table.tsv: a read is correct iff its primary passes the score
threshold and either the primary overlaps the truth OR a secondary alignment tied
on the primary's AS overlaps the truth. As in the pipeline (Benchmark.snake), the
threshold score is recomputed from MD/CIGAR for bwa but read straight from the AS
tag for strobealign (inferred from "strobealign" in the BAM path).
"""
import argparse
import sys
import pysam

# adna_accuracy.py lives in the same directory; reuse its truth parsing, overlap
# check, and score recomputation so this stays in sync with accuracy.tsv.
from adna_accuracy import (
    ReferenceInterval, overlap, parse_gargammel_name,
    recompute_alignment_score, Scores,
)

parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("truth_bed")
parser.add_argument("bam_a", help="expected-correct BAM")
parser.add_argument("bam_b", help="candidate-regressed BAM")
parser.add_argument("out_tsv")
parser.add_argument("--secondary", action="store_true",
                    help="use secondary-alignment correctness (matches accuracy_table.tsv)")
parser.add_argument("--score-threshold", type=int, default=60,
                    help="primary score threshold for --secondary mode (default: 60)")
args = parser.parse_args()

THR = args.score_threshold

truth = {}
with open(args.truth_bed) as fh:
    for line in fh:
        ref_name, start, end, name, _, orientation = line.rstrip("\n").split("\t")
        _, _, _, read_origin, is_deaminated = parse_gargammel_name(name)
        truth[name] = ReferenceInterval(
            ref_name, read_origin, int(start), int(end), is_deaminated
        )


def load(path, with_secondary):
    """Index primaries by query name; collect secondaries only when needed."""
    prim, secs = {}, {}
    with pysam.AlignmentFile(path, "rb") as bam:
        for r in bam:
            if r.is_supplementary:
                continue
            if r.is_secondary:
                if with_secondary:
                    secs.setdefault(r.query_name, []).append(r)
                continue
            prim[r.query_name] = r
    return prim, secs


def recompute_for(path):
    # Benchmark.snake: recompute_as = False for strobealign, True otherwise.
    return "strobealign" not in path.lower()


print(f"Loading BAM_A {args.bam_a}", file=sys.stderr)
A_prim, A_sec = load(args.bam_a, args.secondary)
print(f"Loading BAM_B {args.bam_b}", file=sys.stderr)
B_prim, B_sec = load(args.bam_b, args.secondary)
A_recompute = recompute_for(args.bam_a)
B_recompute = recompute_for(args.bam_b)


def as_tag(rec):
    if rec is None or rec.is_unmapped or not rec.has_tag("AS"):
        return None
    return rec.get_tag("AS")


def primary_overlaps(rec, t):
    return (rec is not None and not rec.is_unmapped
            and rec.reference_name == t.name
            and overlap(rec.reference_start, rec.reference_end, t.start, t.end))


def score_of(rec, recompute):
    a = as_tag(rec)
    if recompute or a is None:
        try:
            return recompute_alignment_score(rec, Scores)
        except Exception:
            return THR - 1
    return a


def correct_primary(rec, t):
    if rec is None or rec.is_unmapped or t.origin != "e":
        return False
    return primary_overlaps(rec, t)


def correct_secondary(prim, secs, name, t, recompute):
    if t.origin != "e":
        return False
    p = prim.get(name)
    if p is None or p.is_unmapped:
        return False
    if score_of(p, recompute) < THR:
        return False
    if primary_overlaps(p, t):
        return True
    p_as = as_tag(p)
    if p_as is None:
        return False
    for r in secs.get(name, []):
        if primary_overlaps(r, t) and as_tag(r) == p_as:
            return True
    return False


def is_correct(prim, secs, name, t, recompute):
    if args.secondary:
        return correct_secondary(prim, secs, name, t, recompute)
    return correct_primary(prim.get(name), t)


candidates = []
for name, t in truth.items():
    if is_correct(A_prim, A_sec, name, t, A_recompute) and \
       not is_correct(B_prim, B_sec, name, t, B_recompute):
        b = B_prim.get(name)
        if b is None or b.is_unmapped:
            dist = -1
        elif b.reference_name == t.name:
            dist = abs(b.reference_start - t.start)
        else:
            dist = 10**9
        candidates.append((name, dist, t, A_prim.get(name), b))

print(f"Found {len(candidates)} candidate reads", file=sys.stderr)
# Sort by distance descending (largest misalignments first, unmapped=-1 sorts last)
candidates.sort(key=lambda x: x[1], reverse=True)

with open(args.out_tsv, "w") as out:
    out.write("read_name\tb_dist\ttruth_chrom\ttruth_start\ttruth_strand\ta_chrom\ta_start\ta_strand\ta_mapq\ta_as\tb_chrom\tb_start\tb_strand\tb_mapq\tb_as\n")
    for name, dist, t, a, b in candidates:
        def fmt(r):
            if r is None or r.is_unmapped:
                return ("*", "-1", "*", "0", "0")
            strand = "-" if r.is_reverse else "+"
            try:
                asv = r.get_tag("AS")
            except KeyError:
                asv = 0
            return (r.reference_name, str(r.reference_start), strand, str(r.mapping_quality), str(asv))
        af = fmt(a)
        bf = fmt(b)
        # truth strand column kept for triage context; it no longer drives correctness
        out.write("\t".join([name, str(dist), t.name, str(t.start), "."] + list(af) + list(bf)) + "\n")

print(f"Wrote {args.out_tsv}", file=sys.stderr)
