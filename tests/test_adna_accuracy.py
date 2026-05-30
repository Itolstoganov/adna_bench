#!/usr/bin/env python3
"""Unit tests for secondary-alignment accuracy in workflow/scripts/adna_accuracy.py.

A read is "secondary-correct" if its primary overlaps the truth OR some secondary
alignment tied on the AS score with the primary overlaps the truth.

Run directly with the project Python (has pysam + adna_accuracy's deps):

    /home/itolstoganov/miniforge3/envs/strobealign-eval/bin/python tests/test_adna_accuracy.py

or, if pytest is available:  pytest tests/test_adna_accuracy.py
"""
import os
import sys
import tempfile

import pysam

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "workflow", "scripts")))
import adna_accuracy as A

CHROM = "chr1"
THRESHOLDS = [25]

# Column order produced by AccuracyResults.to_tsv() with emit_secondary=True:
#   threshold, aligned, correct, aligned_deam, correct_deam,
#   overmapped_bact, overmapped_cont, accuracy_secondary, accuracy_secondary_deam
I_CORRECT = 2
I_ACC_SECONDARY = 7
I_ACC_SECONDARY_DEAM = 8

# gargammel-style read names: ref:orientation:start:end:<len><origin><allele>
READ = "chr1:+:100:140:40e0"            # endogenous, not deaminated, truth 100-140
READ_DEAM = "chr1:+:300:340:40e0_DEAM:5"  # endogenous, deaminated, truth 300-340


def _write_bed(path, entries):
    with open(path, "w") as fh:
        for chrom, start, end, name in entries:
            fh.write(f"{chrom}\t{start}\t{end}\t{name}\t.\t+\n")


def _build_bam(path, records):
    """records: list of (query_name, reference_start, is_secondary, as_tag)."""
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": CHROM, "LN": 1000}]}
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        for qname, pos, secondary, as_tag in records:
            seg = pysam.AlignedSegment(bam.header)
            seg.query_name = qname
            seg.query_sequence = "A" * 40
            seg.flag = 0x100 if secondary else 0
            seg.reference_id = 0
            seg.reference_start = pos
            seg.mapping_quality = 60
            seg.cigartuples = [(0, 40)]
            seg.set_tag("AS", as_tag)
            bam.write(seg)


def _fields(tsv):
    """First (only) threshold row of a to_tsv() output, split into fields.

    Note: do not strip() — trailing tabs are meaningful (blank secondary columns).
    """
    return tsv.split("\n")[0].split("\t")


def _make_dataset(records, deam=False):
    """Write a truth BED + BAM into a temp dir; return (bed_path, bam_path)."""
    d = tempfile.mkdtemp()
    bed = os.path.join(d, "truth.bed")
    bam = os.path.join(d, "pred.bam")
    if deam:
        _write_bed(bed, [(CHROM, 300, 340, READ_DEAM)])
    else:
        _write_bed(bed, [(CHROM, 100, 140, READ)])
    _build_bam(bam, records)
    return bed, bam


def test_secondary_off_is_unchanged():
    # primary wrong (@500), correct secondary tied on AS (@100) — but mode off
    bed, bam = _make_dataset([(READ, 500, False, 80), (READ, 100, True, 80)])
    res = A.measure_accuracy(bed, bam, THRESHOLDS, recompute_score=False)
    fields = _fields(res.to_tsv())
    assert len(fields) == 7, fields                 # no secondary columns emitted
    assert float(fields[I_CORRECT]) == 0.0, fields  # primary is wrong


def test_secondary_tied_rescues_read():
    # primary wrong (@500), correct secondary tied on AS (@100)
    bed, bam = _make_dataset([(READ, 500, False, 80), (READ, 100, True, 80)])
    res = A.measure_accuracy(bed, bam, THRESHOLDS, recompute_score=False,
                             compute_secondary=True, emit_secondary=True)
    fields = _fields(res.to_tsv())
    assert len(fields) == 9, fields
    assert float(fields[I_CORRECT]) == 0.0, fields          # primary still wrong
    assert float(fields[I_ACC_SECONDARY]) == 100.0, fields  # tied secondary rescues it


def test_secondary_tied_rescues_deaminated_read():
    bed, bam = _make_dataset([(READ_DEAM, 700, False, 80), (READ_DEAM, 300, True, 80)], deam=True)
    res = A.measure_accuracy(bed, bam, THRESHOLDS, recompute_score=False,
                             compute_secondary=True, emit_secondary=True)
    fields = _fields(res.to_tsv())
    assert float(fields[I_CORRECT]) == 0.0, fields
    assert float(fields[I_ACC_SECONDARY]) == 100.0, fields
    assert float(fields[I_ACC_SECONDARY_DEAM]) == 100.0, fields


def test_lower_scoring_secondary_not_counted():
    # correct secondary, but its AS (70) is below the primary's (80): not a tie
    bed, bam = _make_dataset([(READ, 500, False, 80), (READ, 100, True, 70)])
    res = A.measure_accuracy(bed, bam, THRESHOLDS, recompute_score=False,
                             compute_secondary=True, emit_secondary=True)
    fields = _fields(res.to_tsv())
    assert float(fields[I_ACC_SECONDARY]) == 0.0, fields


def test_correct_primary_stays_correct():
    # primary already correct (@100): both accuracy and secondary accuracy are 100
    bed, bam = _make_dataset([(READ, 100, False, 80)])
    res = A.measure_accuracy(bed, bam, THRESHOLDS, recompute_score=False,
                             compute_secondary=True, emit_secondary=True)
    fields = _fields(res.to_tsv())
    assert float(fields[I_CORRECT]) == 100.0, fields
    assert float(fields[I_ACC_SECONDARY]) == 100.0, fields


def test_emit_without_compute_is_blank():
    # non-applicable tools (bwaaln/safari): column emitted but not computed -> blank
    bed, bam = _make_dataset([(READ, 500, False, 80), (READ, 100, True, 80)])
    res = A.measure_accuracy(bed, bam, THRESHOLDS, recompute_score=False,
                             compute_secondary=False, emit_secondary=True)
    fields = _fields(res.to_tsv())
    assert len(fields) == 9, fields
    assert fields[I_ACC_SECONDARY] == "", fields
    assert fields[I_ACC_SECONDARY_DEAM] == "", fields


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for t in tests:
        t()
        print(f"PASS {t.__name__}")
    print(f"\nAll {len(tests)} tests passed.")
