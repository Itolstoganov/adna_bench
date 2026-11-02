#!/usr/bin/env python3

import argparse
import re
import sys

import dnaio


def sanitize_fasta(input_path, output_path):
    pattern = re.compile(r"[^ACGTNacgtn]")
    total_changed = 0

    with dnaio.open(input_path) as reader, dnaio.open(output_path, mode="w") as writer:
        for record in reader:
            changed = len(pattern.findall(record.sequence))
            total_changed += changed
            if changed:
                record.sequence = pattern.sub("N", record.sequence)
            writer.write(record)

    print(f"Sanitized {input_path}: {total_changed} characters changed to N", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Replace non-ACGTN characters in FASTA sequences with N"
    )
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file")
    args = parser.parse_args()

    sanitize_fasta(args.input, args.output)


if __name__ == "__main__":
    main()
