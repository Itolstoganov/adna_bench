#!/usr/bin/env python3
"""Parse gargammel simulation reads into per-read ground truth.

Gargammel/fragSim encodes each read's origin in the FASTQ header as
``<ref_id>:<orientation>:<start>:<end>:<len><class>...`` where the class char is
one of ``b`` (bacterial), ``c`` (contaminant) or ``e`` (endogenous). The
Itolstoganov gargammel fork additionally appends ``_DEAM:<signed positions>``
(via ``deamSim -name``) listing deaminated bases (1-based from the 5' end, or
negative counting from the 3' end of the fragment).

Outputs (all optional except the fragment BED):
  -o/--output         fragment BED (all reads): ref_id start end read_id 0 strand
  -D/--damage-output  per-base deamination BED, restricted to --damage-classes
  -t/--truth          per-read TSV: read_id class taxon taxid ref_id start end strand n_deam

Both paired-end (``-1``/``-2``) and single-end (``-1`` only) inputs are supported.
For bacterial reads, ``--contig2taxon`` maps the source contig (ref_id) to a
species/taxid; endogenous and contaminant reads are labelled with the
configurable ``--endo-taxon``/``--cont-taxon`` (default human).
"""

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import dnaio


READ_CLASS_DICT = {'b': "bact", 'c': "cont", 'e': "endo"}

# Matches the suffix deamSim appends when invoked with -name.
# Format: _DEAM:<csv-of-signed-ints>, optionally followed by art_illumina's
# amplicon counter "-<int>" at the very end of the (pair-stripped) header.
DEAM_RE = re.compile(r'_DEAM:(-?\d+(?:,-?\d+)*)(?:-\d+)?$')

READ_TYPE_RE = re.compile(r'(\d+)(a|b|c|d|e)(.*)')


def parse_fastq_header(header: str) -> Tuple[str, str, int, int, str, List[int]]:
    header = header.lstrip('@').strip()

    parts = header.split(':')
    if len(parts) < 5:
        raise ValueError(f"Invalid header format: {header}")

    ref_id = parts[0]
    orientation = parts[1]
    start = int(parts[2])
    end = int(parts[3])
    read_type_str = ':'.join(parts[4:])

    deam_positions: List[int] = []
    deam_match = DEAM_RE.search(read_type_str)
    if deam_match:
        deam_positions = [int(p) for p in deam_match.group(1).split(',')]
        read_type_str = read_type_str[:deam_match.start()]

    match = READ_TYPE_RE.search(read_type_str)
    if not match:
        raise ValueError(f"Cannot parse read type from header: {header}")
    _, read_type, _ = match.groups()

    return (
        ref_id,
        orientation,
        start,
        end,
        READ_CLASS_DICT.get(read_type, "unknown"),
        deam_positions,
    )


def normalize_read_name(read_name: str) -> str:
    if read_name.endswith('/1') or read_name.endswith('/2'):
        return read_name[:-2]
    return read_name


def damage_positions_to_ref(
    frag_start: int,
    frag_end: int,
    orientation: str,
    deam_positions: List[int],
) -> List[int]:
    """Map damage positions (1-based from 5', or negative from 3' of the
    fragment as emitted by deamSim) to 0-based reference coordinates."""
    frag_len = frag_end - frag_start
    ref_positions: List[int] = []
    for pos in deam_positions:
        if pos > 0:
            frag_offset = pos - 1
        else:
            frag_offset = frag_len + pos  # pos negative; -1 = last base
        if frag_offset < 0 or frag_offset >= frag_len:
            continue
        if orientation == '+':
            ref_positions.append(frag_start + frag_offset)
        else:
            ref_positions.append(frag_end - 1 - frag_offset)
    return ref_positions


def _record_from_name(name: str) -> Dict:
    normalized_id = normalize_read_name(name.split()[0])
    ref_id, orientation, start, end, read_class, deam_positions = \
        parse_fastq_header(normalized_id)
    return {
        'read_id': normalized_id,
        'reference_id': ref_id,
        'start': start,
        'end': end,
        'orientation': orientation,
        'reference_class': read_class,
        'deam_positions': deam_positions,
    }


def parse_paired(fastq1_path: Path, fastq2_path: Path) -> List[Dict]:
    reads_data: List[Dict] = []
    with dnaio.open(fastq1_path, fastq2_path, mode="r") as reader:
        for record1, record2 in reader:
            id1 = normalize_read_name(record1.name.split()[0])
            id2 = normalize_read_name(record2.name.split()[0])
            if id1 != id2:
                raise ValueError(
                    f"Read pair names do not match: {record1.name} / {record2.name}"
                )
            try:
                reads_data.append(_record_from_name(record1.name))
            except Exception as exc:
                print(f"Could not parse header for read {record1.name}", file=sys.stderr)
                raise exc
    return reads_data


def parse_single(fastq_path: Path) -> List[Dict]:
    reads_data: List[Dict] = []
    with dnaio.open(fastq_path, mode="r") as reader:
        for record in reader:
            try:
                reads_data.append(_record_from_name(record.name))
            except Exception as exc:
                print(f"Could not parse header for read {record.name}", file=sys.stderr)
                raise exc
    return reads_data


def load_contig2taxon(path: Optional[Path]) -> Dict[str, Tuple[str, str]]:
    """contig_id -> (taxon, taxid) from a TSV with columns contig, taxon, taxid."""
    mapping: Dict[str, Tuple[str, str]] = {}
    if path is None:
        return mapping
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            contig = fields[0]
            taxon = fields[1] if len(fields) > 1 else contig
            taxid = fields[2] if len(fields) > 2 else "NA"
            mapping[contig] = (taxon, taxid)
    return mapping


def resolve_taxon(
    read: Dict,
    contig2taxon: Dict[str, Tuple[str, str]],
    endo_taxon: str,
    endo_taxid: str,
    cont_taxon: str,
    cont_taxid: str,
) -> Tuple[str, str]:
    cls = read['reference_class']
    if cls == 'bact':
        return contig2taxon.get(read['reference_id'], (read['reference_id'], "NA"))
    if cls == 'endo':
        return (endo_taxon, endo_taxid)
    if cls == 'cont':
        return (cont_taxon, cont_taxid)
    return (read['reference_id'], "NA")


def output_bed(reads_data: List[Dict], output_path: Path) -> None:
    with open(output_path, 'w') as f:
        for read in reads_data:
            f.write('\t'.join([
                str(read['reference_id']),
                str(read['start']),
                str(read['end']),
                str(read['read_id']),
                "0",
                str(read['orientation']),
            ]) + '\n')


def output_damage_bed(
    reads_data: List[Dict], output_path: Path, damage_classes: set
) -> int:
    count = 0
    with open(output_path, 'w') as f:
        for read in reads_data:
            if read['reference_class'] not in damage_classes or not read['deam_positions']:
                continue
            ref_positions = damage_positions_to_ref(
                read['start'], read['end'], read['orientation'], read['deam_positions']
            )
            for ref_pos in ref_positions:
                f.write('\t'.join([
                    str(read['reference_id']),
                    str(ref_pos),
                    str(ref_pos + 1),
                    str(read['read_id']),
                    "0",
                    str(read['orientation']),
                ]) + '\n')
                count += 1
    return count


def output_truth(
    reads_data: List[Dict],
    output_path: Path,
    contig2taxon: Dict[str, Tuple[str, str]],
    endo_taxon: str,
    endo_taxid: str,
    cont_taxon: str,
    cont_taxid: str,
    pass_label: Optional[str],
) -> None:
    header = ["read_id", "class", "taxon", "taxid", "ref_id",
              "start", "end", "strand", "n_deam"]
    if pass_label is not None:
        header.append("pass")
    with open(output_path, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for read in reads_data:
            taxon, taxid = resolve_taxon(
                read, contig2taxon, endo_taxon, endo_taxid, cont_taxon, cont_taxid
            )
            row = [
                read['read_id'],
                read['reference_class'],
                taxon,
                taxid,
                read['reference_id'],
                str(read['start']),
                str(read['end']),
                read['orientation'],
                str(len(read['deam_positions'])),
            ]
            if pass_label is not None:
                row.append(pass_label)
            f.write('\t'.join(row) + '\n')


def parse_damage_classes(value: str) -> set:
    if value.strip().lower() == "all":
        return {"bact", "cont", "endo"}
    return {c.strip() for c in value.split(',') if c.strip()}


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Parse gargammel simulation reads into per-read ground truth',
    )
    parser.add_argument('-d', '--dataset-dir', type=Path, default=None,
                        help='Optional dataset dir (existence check only)')
    parser.add_argument('-1', '--read1', type=Path, required=True)
    parser.add_argument('-2', '--read2', type=Path, default=None,
                        help='R2 for paired-end; omit for single-end')
    parser.add_argument('-o', '--output', type=Path, required=True,
                        help='Output fragment BED (all reads)')
    parser.add_argument('-D', '--damage-output', type=Path, default=None,
                        help='Optional per-base deamination BED')
    parser.add_argument('-t', '--truth', type=Path, default=None,
                        help='Optional per-read truth TSV (taxon-labelled)')
    parser.add_argument('--damage-classes', default='endo',
                        help='Comma list or "all"; classes to emit damage for '
                             '(default: endo, for sim_bench back-compat)')
    parser.add_argument('--contig2taxon', type=Path, default=None,
                        help='TSV contig<TAB>taxon<TAB>taxid for bacterial reads')
    parser.add_argument('--endo-taxon', default='Homo_sapiens_ancient')
    parser.add_argument('--endo-taxid', default='9606')
    parser.add_argument('--cont-taxon', default='Homo_sapiens_modern')
    parser.add_argument('--cont-taxid', default='9606')
    parser.add_argument('--pass-label', default=None,
                        help='If set, add a "pass" column to the truth TSV')
    args = parser.parse_args()

    if args.dataset_dir is not None and not args.dataset_dir.exists():
        print(f"Error: Dataset directory {args.dataset_dir} does not exist", file=sys.stderr)
        sys.exit(1)
    if not args.read1.exists():
        print(f"Error: R1 file {args.read1} does not exist", file=sys.stderr)
        sys.exit(1)
    if args.read2 is not None and not args.read2.exists():
        print(f"Error: R2 file {args.read2} does not exist", file=sys.stderr)
        sys.exit(1)

    if args.read2 is not None:
        print(f"Parsing paired-end files:\n  R1: {args.read1}\n  R2: {args.read2}",
              file=sys.stderr)
        reads_data = parse_paired(args.read1, args.read2)
    else:
        print(f"Parsing single-end file:\n  R1: {args.read1}", file=sys.stderr)
        reads_data = parse_single(args.read1)
    print(f"  Extracted {len(reads_data)} reads", file=sys.stderr)

    contig2taxon = load_contig2taxon(args.contig2taxon)
    damage_classes = parse_damage_classes(args.damage_classes)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    output_bed(reads_data, args.output)

    if args.damage_output is not None:
        args.damage_output.parent.mkdir(parents=True, exist_ok=True)
        damage_count = output_damage_bed(reads_data, args.damage_output, damage_classes)
        print(f"  Wrote {damage_count} damage positions ({sorted(damage_classes)}) "
              f"to {args.damage_output}", file=sys.stderr)

    if args.truth is not None:
        args.truth.parent.mkdir(parents=True, exist_ok=True)
        output_truth(reads_data, args.truth, contig2taxon,
                     args.endo_taxon, args.endo_taxid,
                     args.cont_taxon, args.cont_taxid, args.pass_label)
        print(f"  Wrote per-read truth to {args.truth}", file=sys.stderr)

    if reads_data:
        class_counts: Dict[str, int] = defaultdict(int)
        damaged_reads = 0
        total_damage_sites = 0
        for read in reads_data:
            class_counts[read['reference_class']] += 1
            if read['deam_positions']:
                damaged_reads += 1
                total_damage_sites += len(read['deam_positions'])

        print("\nReference class distribution:", file=sys.stderr)
        for ref_class, count in sorted(class_counts.items()):
            percentage = 100 * count / len(reads_data)
            print(f"  {ref_class}: {count} ({percentage:.2f}%)", file=sys.stderr)
        print(
            f"\nReads with deamination: {damaged_reads} "
            f"({100 * damaged_reads / len(reads_data):.2f}%), "
            f"total damage sites (incl. duplicate positions): {total_damage_sites}",
            file=sys.stderr,
        )


if __name__ == '__main__':
    main()
