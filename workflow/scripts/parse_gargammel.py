#!/usr/bin/env python3

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


def parse_paired(fastq1_path: Path, fastq2_path: Path) -> List[Dict]:
    reads_data: List[Dict] = []
    with dnaio.open(fastq1_path, fastq2_path, mode="r") as reader:
        for record1, record2 in reader:
            normalized_id1 = normalize_read_name(record1.name.split()[0])
            normalized_id2 = normalize_read_name(record2.name.split()[0])
            if normalized_id1 != normalized_id2:
                raise ValueError(
                    f"Read pair names do not match: {record1.name} / {record2.name}"
                )
            try:
                ref_id, orientation, start, end, read_class, deam_positions = \
                    parse_fastq_header(normalized_id1)
            except Exception as exc:
                print(
                    f"Could not parse header for read {record1.name}", file=sys.stderr
                )
                raise exc

            reads_data.append({
                'read_id': normalized_id1,
                'reference_id': ref_id,
                'start': start,
                'end': end,
                'orientation': orientation,
                'reference_class': read_class,
                'deam_positions': deam_positions,
            })

    return reads_data


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


def output_damage_bed(reads_data: List[Dict], output_path: Path) -> int:
    count = 0
    with open(output_path, 'w') as f:
        for read in reads_data:
            if read['reference_class'] != 'endo' or not read['deam_positions']:
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


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Parse Gargammel paired-end simulation reads and extract metadata',
    )
    parser.add_argument('-d', '--dataset-dir', type=Path, required=True)
    parser.add_argument('-1', '--read1', type=Path, required=True)
    parser.add_argument('-2', '--read2', type=Path, required=True)
    parser.add_argument('-o', '--output', type=Path, required=True,
                        help='Output BED with endogenous fragment positions')
    parser.add_argument('-D', '--damage-output', type=Path, default=None,
                        help='Optional BED listing per-base deamination positions')
    args = parser.parse_args()

    if not args.dataset_dir.exists():
        print(f"Error: Dataset directory {args.dataset_dir} does not exist", file=sys.stderr)
        sys.exit(1)
    if not args.read1.exists():
        print(f"Error: R1 file {args.read1} does not exist", file=sys.stderr)
        sys.exit(1)
    if not args.read2.exists():
        print(f"Error: R2 file {args.read2} does not exist", file=sys.stderr)
        sys.exit(1)

    print(f"Parsing paired-end files:", file=sys.stderr)
    print(f"  R1: {args.read1}", file=sys.stderr)
    print(f"  R2: {args.read2}", file=sys.stderr)

    reads_data = parse_paired(args.read1, args.read2)
    print(f"  Extracted {len(reads_data)} valid paired reads", file=sys.stderr)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    output_bed(reads_data, args.output)

    if args.damage_output is not None:
        args.damage_output.parent.mkdir(parents=True, exist_ok=True)
        damage_count = output_damage_bed(reads_data, args.damage_output)
        print(f"  Wrote {damage_count} damage positions to {args.damage_output}",
              file=sys.stderr)

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
