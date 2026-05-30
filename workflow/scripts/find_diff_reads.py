#!/usr/bin/env python
"""Find reads correct in BAM A but incorrect in BAM B, relative to ground truth.

Correctness reuses the pipeline's authoritative definition from adna_accuracy.py:
a read is correct iff it is mapped, on the same reference as the truth, its alignment
interval overlaps the truth interval, and the read is endogenous (origin "e").
"""
import sys
import pysam

# adna_accuracy.py lives in the same directory; reuse its truth parsing and the
# overlap-based correctness check so this stays in sync with accuracy.tsv.
from adna_accuracy import ReferenceInterval, overlap, parse_gargammel_name

TRUTH_BED = sys.argv[1]
BAM_A = sys.argv[2]  # expected correct
BAM_B = sys.argv[3]  # candidate regressed
OUT_LIST = sys.argv[4]

truth = {}
with open(TRUTH_BED) as fh:
    for line in fh:
        ref_name, start, end, name, _, orientation = line.rstrip("\n").split("\t")
        _, _, _, read_origin, is_deaminated = parse_gargammel_name(name)
        truth[name] = ReferenceInterval(
            ref_name, read_origin, int(start), int(end), is_deaminated
        )


def bam_index(path):
    d = {}
    with pysam.AlignmentFile(path, "rb") as bam:
        for r in bam:
            if r.is_secondary or r.is_supplementary:
                continue
            d[r.query_name] = r
    return d


print(f"Loading BAM_A {BAM_A}", file=sys.stderr)
A = bam_index(BAM_A)
print(f"Loading BAM_B {BAM_B}", file=sys.stderr)
B = bam_index(BAM_B)


def correct(rec, truth_rec):
    if rec is None or rec.is_unmapped:
        return False
    if truth_rec.origin != "e":
        return False
    if rec.reference_name != truth_rec.name:
        return False
    return overlap(rec.reference_start, rec.reference_end, truth_rec.start, truth_rec.end)


candidates = []
for name, truth_rec in truth.items():
    a = A.get(name)
    b = B.get(name)
    if correct(a, truth_rec) and not correct(b, truth_rec):
        # distance of b from truth start
        if b is None or b.is_unmapped:
            dist = -1  # unmapped
        elif b.reference_name == truth_rec.name:
            dist = abs(b.reference_start - truth_rec.start)
        else:
            dist = 10**9
        candidates.append((name, dist, truth_rec, a, b))

print(f"Found {len(candidates)} candidate reads", file=sys.stderr)
# Sort by distance descending (largest misalignments first, unmapped=-1 sorts last)
candidates.sort(key=lambda x: x[1], reverse=True)

with open(OUT_LIST, "w") as out:
    out.write("read_name\tb_dist\ttruth_chrom\ttruth_start\ttruth_strand\ta_chrom\ta_start\ta_strand\ta_mapq\ta_as\tb_chrom\tb_start\tb_strand\tb_mapq\tb_as\n")
    for name, dist, truth_rec, a, b in candidates:
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
        out.write("\t".join([name, str(dist), truth_rec.name, str(truth_rec.start), "."] + list(af) + list(bf)) + "\n")

print(f"Wrote {OUT_LIST}", file=sys.stderr)
