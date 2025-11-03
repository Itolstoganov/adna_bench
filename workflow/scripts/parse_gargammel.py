#!/usr/bin/env python3

import argparse
import gzip
import os
import sys
from collections import defaultdict
from pathlib import Path
import re
from typing import Dict, List, Set, Tuple
from xopen import xopen

from Bio import SeqIO


def parse_fastq_header(header: str) -> Tuple[str, str, int, int, str]:
    read_class_dict = {'b': "bact", 'c': "cont", 'e': "endo"}

    header = header.lstrip('@').strip()

    parts = header.split(':')
    if len(parts) < 5:
        raise ValueError(f"Invalid header format: {header}")

    ref_id = parts[0]
    orientation = parts[1]
    start = int(parts[2])
    end = int(parts[3])
    read_type_str = ':'.join(parts[4:])
    
    read_type_pattern = r'(\d+)(a|b|c|d|e)(.*)'
    match = re.search(read_type_pattern, read_type_str)
    read_len, read_type, _ = match.groups()


    return ref_id, orientation, start, end, read_class_dict.get(read_type, "unknown")


def build_reference_map(dataset_dir: Path) -> Dict[str, Set[str]]:
    ref_map = defaultdict(set)

    for ref_class in ['endo', 'cont', 'bact']:
        class_dir = dataset_dir / ref_class

        if not class_dir.exists():
            print(f"Warning: Directory {class_dir} does not exist", file=sys.stderr)
            continue

        fasta_files = list(class_dir.glob('*.fa')) + list(class_dir.glob('*.fasta')) + list(class_dir.glob('*.fna'))

        for fasta_file in fasta_files:
            with open(fasta_file, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        contig_id = line[1:].split()[0].strip()
                        ref_map[contig_id].add(ref_class)

    return ref_map


def normalize_read_name(read_name: str) -> str:
    if not (read_name.endswith('/1') or read_name.endswith('/2')):
        raise ValueError("Read name format incorrect")
    return read_name[:-2]


def parse_paired(fastq1_path: Path, fastq2_path: Path) -> List[Dict[str, any]]:
    reads_data = []
    read_count = 0

    with xopen(fastq1_path) as f1, xopen(fastq2_path) as f2:
        for record1, record2 in zip(SeqIO.parse(f1, "fastq"), SeqIO.parse(f2, "fastq")):
            read_count += 1

            read_id1 = record1.id
            read_id2 = record2.id

            normalized_id1 = normalize_read_name(read_id1)
            normalized_id2 = normalize_read_name(read_id2)

            try:
                if normalized_id1 != normalized_id2:
                    raise ValueError("Read pair names do not match")
                ref_id, orientation, start, end, read_class = parse_fastq_header(record1.id)

                # Check if the reference contig is present in the dataset
                # if ref_id not in ref_map:
                #     raise ValueError(f"Contig {ref_id} not found in the dataset")
                # elif read_class not in ref_map[ref_id]:
                #     raise ValueError(f"Reference {read_class} does not contain contig {ref_id}")

                reads_data.append({
                    'read_id': normalized_id1,
                    'reference_id': ref_id,
                    'start': start,
                    'end': end,
                    'orientation': orientation,
                    'reference_class': read_class
                })

            except Exception as e:
                print(f"Could not parse header for read {record1.id}, {record2.id}",
                      file=sys.stderr)
                raise e

    return reads_data


def write_output(reads_data: List[Dict[str, any]], output_path: Path):
    delimiter = '\t' 

    with open(output_path, 'w') as f:
        f.write(delimiter.join([
            'read_id',
            'reference_id',
            'start',
            'end',
            'orientation',
            'reference_class'
        ]) + '\n')

        for read in reads_data:
            f.write(delimiter.join([
                str(read['read_id']),
                str(read['reference_id']),
                str(read['start']),
                str(read['end']),
                str(read['orientation']),
                str(read['reference_class'])
            ]) + '\n')


def main():
    parser = argparse.ArgumentParser(
        description='Parse Gargammel paired-end simulation reads and extract metadata',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # Parse paired-end read files from a dataset
  python parse_gargammel.py -d datasets/my_dataset \\
      -1 datasets/my_dataset/simulated_s1.fq.gz \\
      -2 datasets/my_dataset/simulated_s2.fq.gz \\
      -o output/read_metadata.tsv
        """
    )

    parser.add_argument(
        '-d', '--dataset-dir',
        type=Path,
        required=True,
        help='Path to dataset directory containing endo/, cont/, and bact/ subdirectories'
    )

    parser.add_argument(
        '-1', '--read1',
        type=Path,
        required=True,
        help='Input R1 FASTQ file (can be gzipped)'
    )

    parser.add_argument(
        '-2', '--read2',
        type=Path,
        required=True,
        help='Input R2 FASTQ file (can be gzipped)'
    )

    parser.add_argument(
        '-o', '--output',
        type=Path,
        required=True,
        help='Output file path'
    )

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

    # print(f"Building reference map from {args.dataset_dir}...", file=sys.stderr)
    # ref_map = build_reference_map(args.dataset_dir)
    # print(f"Found {len(ref_map)} reference contigs across {len(set(ref_map.values()))} classes",
    #       file=sys.stderr)

    print(f"Parsing paired-end files:", file=sys.stderr)
    print(f"  R1: {args.read1}", file=sys.stderr)
    print(f"  R2: {args.read2}", file=sys.stderr)

    reads_data = parse_paired(args.read1, args.read2)
    print(f"  Extracted {len(reads_data)} valid paired reads", file=sys.stderr)

    print(f"Writing output to {args.output}...", file=sys.stderr)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_output(reads_data, args.output)
    print(f"Done! Total reads processed: {len(reads_data)}", file=sys.stderr)

    if reads_data:
        class_counts = defaultdict(int)
        for read in reads_data:
            class_counts[read['reference_class']] += 1

        print("\nReference class distribution:", file=sys.stderr)
        for ref_class, count in sorted(class_counts.items()):
            percentage = 100 * count / len(reads_data)
            print(f"  {ref_class}: {count} ({percentage:.2f}%)", file=sys.stderr)


if __name__ == '__main__':
    main()
