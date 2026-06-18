#!/usr/bin/env python3
"""Parse rymer/strobealign --trace logs into a per-read dataset.

The aligner, when run with --trace (single threaded), emits one block per read
on stderr:

    Query: <name>
    length=<len>
    we have <f> + <r> randstrobes
    Chaining <n> anchors
    Found <h> hits (<rescued> rescued, <filtered> filtered):
    querypos count (p=partial, F=filtered)
    <qpos> <p|space><count> <F|space>      # one per seed hit
    ...
    Found <m> NAMs
    - Nam(ref_id=<N>, query: <a>..<b>, ref: <c>..<d>, rc=<0/1>, score=<s>)
    ...

For every read this script records:
  1. number of full seed hits with >=1 index occurrence
  2. number of partial (rymer) seed hits with >=1 index occurrence
  3. full matches as a list of `query_pos:occurrence_count` pairs
  4. partial matches as a list of `query_pos:occurrence_count` pairs
  5. read deamination positions (decoded from the gargammel read name)
  6. whether a chain (NAM) overlaps the read's ground-truth position, and whether
     that correct chain is one of the highest-scoring chains for the read
     (correct_chain_is_top -> the read is "accurately mapped").

It then prints aggregate statistics -- the joint distribution of (n_full,
n_partial) per read plus mean full/partial occurrence -- for several read
categories (all / deaminated / accurately-vs-incorrectly mapped splits).

Ground-truth coordinates come from the dataset's ref_ground_truth.bed (parsed
with the same helpers as adna_accuracy.py). NAM ref_id is a numeric contig
index = fasta order, so it is mapped to a chromosome name via the ref.fa header
order.
"""

import argparse
import re
import sys
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Reuse the existing gargammel/ground-truth helpers living next to this script.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from parse_gargammel import normalize_read_name, parse_fastq_header  # noqa: E402
from adna_accuracy import overlap, parse_gargammel_name  # noqa: E402


# --- trace line patterns ----------------------------------------------------

QUERY_RE = re.compile(r"^Query: (.+)$")
LENGTH_RE = re.compile(r"^length=(\d+)$")
HIT_HEADER = "querypos count (p=partial, F=filtered)"
# per-hit line: "{:6} {}{:6} {}" -> qpos, (p|space), count, (F|space)
HIT_RE = re.compile(r"^\s*(\d+)\s+(p)?\s*(\d+)\s*(F)?\s*$")
NAM_RE = re.compile(
    r"^- Nam\(ref_id=(\d+), query: (\d+)\.\.(\d+), "
    r"ref: (\d+)\.\.(\d+), rc=(\d+), score=([\d.]+)\)"
)


@dataclass
class Hit:
    query_start: int
    count: int
    is_partial: bool
    is_filtered: bool


@dataclass
class Nam:
    ref_id: int
    query_start: int
    query_end: int
    ref_start: int
    ref_end: int
    is_revcomp: bool
    score: float


@dataclass
class ReadBlock:
    name: str
    length: int = 0
    hits: List[Hit] = field(default_factory=list)
    nams: List[Nam] = field(default_factory=list)


def iter_read_blocks(trace_path: Path):
    """Yield ReadBlock objects, one per `Query:` block in the trace log."""
    block: Optional[ReadBlock] = None
    in_hit_table = False
    with open(trace_path, encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")

            m = QUERY_RE.match(line)
            if m:
                if block is not None:
                    yield block
                block = ReadBlock(name=m.group(1).strip())
                in_hit_table = False
                continue

            if block is None:
                continue  # pre-read noise (indexing/info lines)

            m = LENGTH_RE.match(line)
            if m:
                block.length = int(m.group(1))
                continue

            if line == HIT_HEADER:
                in_hit_table = True
                continue

            if in_hit_table:
                m = HIT_RE.match(line)
                if m:
                    block.hits.append(
                        Hit(
                            query_start=int(m.group(1)),
                            count=int(m.group(3)),
                            is_partial=m.group(2) == "p",
                            is_filtered=m.group(4) == "F",
                        )
                    )
                    continue
                in_hit_table = False  # table ended

            m = NAM_RE.match(line)
            if m:
                block.nams.append(
                    Nam(
                        ref_id=int(m.group(1)),
                        query_start=int(m.group(2)),
                        query_end=int(m.group(3)),
                        ref_start=int(m.group(4)),
                        ref_end=int(m.group(5)),
                        is_revcomp=m.group(6) == "1",
                        score=float(m.group(7)),
                    )
                )

    if block is not None:
        yield block


def load_ref_names(ref_fa: Path) -> List[str]:
    """Contig names in fasta order; index == aligner ref_id."""
    names: List[str] = []
    with open(ref_fa, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                names.append(line[1:].split()[0])
    return names


def load_ground_truth(gt_path: Path) -> Dict[str, Tuple[str, int, int, str, str]]:
    """read_name -> (chrom, start, end, strand, origin). Keyed by normalized name."""
    truth: Dict[str, Tuple[str, int, int, str, str]] = {}
    with open(gt_path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 6:
                continue
            chrom, start, end, read_name, _, strand = fields[:6]
            _, _, _, origin, _ = parse_gargammel_name(read_name)
            key = normalize_read_name(read_name.split()[0])
            truth[key] = (chrom, int(start), int(end), strand, origin)
    return truth


def read_deam_positions(read_name: str) -> List[int]:
    """0-based deamination offsets along the fragment (5'->3'), from the read name.

    deamSim emits positions 1-based from the 5' end (or negative from the 3'
    end); convert to 0-based fragment offsets. Returns [] if the name has no
    _DEAM tag or cannot be parsed.
    """
    try:
        _, _, start, end, _, deam = parse_fastq_header(read_name)
    except (ValueError, IndexError):
        return []
    frag_len = end - start
    offsets: List[int] = []
    for pos in deam:
        off = pos - 1 if pos > 0 else frag_len + pos
        if 0 <= off < frag_len:
            offsets.append(off)
    return sorted(offsets)


def join_ints(values) -> str:
    return ",".join(str(v) for v in values)


def join_matches(hits) -> str:
    """List of `query_pos:occurrence_count` pairs, comma-separated."""
    return ",".join(f"{h.query_start}:{h.count}" for h in hits)


def find_correct_chain(
    nams: List[Nam], ref_names: List[str], gt: Tuple[str, int, int, str, str]
) -> Optional[Nam]:
    """Highest-scoring NAM whose contig+interval overlaps the ground truth."""
    gt_chrom, gt_start, gt_end = gt[0], gt[1], gt[2]
    best: Optional[Nam] = None
    for nam in nams:
        if nam.ref_id >= len(ref_names):
            continue
        if ref_names[nam.ref_id] != gt_chrom:
            continue
        if not overlap(nam.ref_start, nam.ref_end, gt_start, gt_end):
            continue
        if best is None or nam.score > best.score:
            best = nam
    return best


@dataclass
class ReadStat:
    """Per-read summary used for the aggregate statistics report."""
    is_deaminated: bool
    is_accurate: bool  # correct chain exists AND ties the highest NAM score
    full_counts: List[int]  # index occurrence count of each full seed hit
    partial_counts: List[int]  # index occurrence count of each partial seed hit

    @property
    def n_full(self) -> int:
        return len(self.full_counts)

    @property
    def n_partial(self) -> int:
        return len(self.partial_counts)


COLUMNS = [
    "read_name", "read_len", "origin", "is_deaminated",
    "gt_chrom", "gt_start", "gt_end", "gt_strand",
    "n_full", "n_partial",
    "full_matches", "partial_matches",
    "deam_positions",
    "has_correct_chain", "correct_chain_is_top", "chain_ref", "chain_score",
]


def build_row(
    block: ReadBlock,
    ref_names: List[str],
    truth: Dict[str, Tuple[str, int, int, str, str]],
) -> Tuple[Dict[str, object], ReadStat]:
    full_hits = [h for h in block.hits if not h.is_partial and h.count >= 1]
    partial_hits = [h for h in block.hits if h.is_partial and h.count >= 1]

    key = normalize_read_name(block.name.split()[0])
    gt = truth.get(key)
    deam = read_deam_positions(block.name)
    is_deaminated = bool(deam)

    row = {
        "read_name": block.name,
        "read_len": block.length,
        "origin": gt[4] if gt else "",
        "is_deaminated": int(is_deaminated),
        "gt_chrom": gt[0] if gt else "",
        "gt_start": gt[1] if gt else "",
        "gt_end": gt[2] if gt else "",
        "gt_strand": gt[3] if gt else "",
        "n_full": len(full_hits),
        "n_partial": len(partial_hits),
        "full_matches": join_matches(full_hits),
        "partial_matches": join_matches(partial_hits),
        "deam_positions": join_ints(deam),
        "has_correct_chain": 0,
        "correct_chain_is_top": 0,
        "chain_ref": "", "chain_score": "",
    }

    is_top = False
    if gt is not None:
        chain = find_correct_chain(block.nams, ref_names, gt)
        if chain is not None:
            # Is the correct (ground-truth-overlapping) chain among the
            # highest-scoring chains for this read? Scores are parsed from the
            # 2-decimal trace, so a small epsilon covers tie comparisons.
            max_score = max(n.score for n in block.nams)
            is_top = chain.score >= max_score - 1e-6
            row.update({
                "has_correct_chain": 1,
                "correct_chain_is_top": int(is_top),
                "chain_ref": f"{ref_names[chain.ref_id]}:{chain.ref_start}..{chain.ref_end}",
                "chain_score": chain.score,
            })

    stat = ReadStat(
        is_deaminated=is_deaminated,
        is_accurate=is_top,
        full_counts=[h.count for h in full_hits],
        partial_counts=[h.count for h in partial_hits],
    )
    return row, stat


# --- aggregate statistics report -------------------------------------------

# Read categories for the statistics report. "Accurately mapped" == the correct
# chain exists and is one of the highest-scoring chains (is_accurate).
CATEGORIES = [
    ("all reads", lambda s: True),
    ("deaminated reads", lambda s: s.is_deaminated),
    ("non-deam, top correct chain", lambda s: not s.is_deaminated and s.is_accurate),
    ("deam, top correct chain", lambda s: s.is_deaminated and s.is_accurate),
    ("non-deam, unmapped/incorrect", lambda s: not s.is_deaminated and not s.is_accurate),
    ("deam, unmapped/incorrect", lambda s: s.is_deaminated and not s.is_accurate),
]

JOINT_CAP = 10  # collapse seed counts >= JOINT_CAP into a single "N+" bucket


def _mean(values: List[int]) -> str:
    return f"{sum(values) / len(values):.3f}" if values else "n/a"


def format_joint(records: List[ReadStat]) -> str:
    """Joint distribution of (n_full, n_partial) per read, as a capped matrix.

    Rows are the number of full matches, columns the number of partial matches;
    counts >= JOINT_CAP are bucketed into "N+". Cell = number of reads.
    """
    counter = Counter(
        (min(s.n_full, JOINT_CAP), min(s.n_partial, JOINT_CAP)) for s in records
    )
    if not counter:
        return "    (no reads)"
    max_f = max(f for f, _ in counter)
    max_p = max(p for _, p in counter)

    def label(v: int) -> str:
        return f"{v}+" if v >= JOINT_CAP else str(v)

    col_labels = [label(p) for p in range(max_p + 1)]
    width = max(6, max(len(c) for c in col_labels) + 1)
    lines = []
    header = "full\\part".ljust(9) + "".join(c.rjust(width) for c in col_labels)
    lines.append("    " + header)
    for f in range(max_f + 1):
        cells = "".join(str(counter.get((f, p), 0)).rjust(width) for p in range(max_p + 1))
        lines.append("    " + label(f).ljust(9) + cells)
    return "\n".join(lines)


def print_report(records: List[ReadStat]) -> None:
    print("\n========== seed-match statistics ==========")
    for name, predicate in CATEGORIES:
        subset = [s for s in records if predicate(s)]
        full_occ = [c for s in subset for c in s.full_counts]
        partial_occ = [c for s in subset for c in s.partial_counts]
        print(f"\n### {name}  (n={len(subset)})")
        print(f"  mean full-match occurrence:    {_mean(full_occ)}  (over {len(full_occ)} hits)")
        print(f"  mean partial-match occurrence: {_mean(partial_occ)}  (over {len(partial_occ)} hits)")
        print("  joint distribution of (n_full rows x n_partial cols):")
        print(format_joint(subset))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--trace", type=Path, required=True,
                        help="Aligner --trace log (stderr)")
    parser.add_argument("--ground-truth", type=Path, required=True,
                        help="ref_ground_truth.bed for the dataset")
    parser.add_argument("--ref", type=Path, required=True,
                        help="ref.fa (for ref_id -> contig name mapping)")
    parser.add_argument("--out", type=Path, required=True,
                        help="Output per-read TSV")
    args = parser.parse_args()

    ref_names = load_ref_names(args.ref)
    truth = load_ground_truth(args.ground_truth)
    print(f"Loaded {len(ref_names)} contigs, {len(truth)} ground-truth reads",
          file=sys.stderr)

    records: List[ReadStat] = []
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as out:
        out.write("\t".join(COLUMNS) + "\n")
        for block in iter_read_blocks(args.trace):
            row, stat = build_row(block, ref_names, truth)
            out.write("\t".join(str(row[c]) for c in COLUMNS) + "\n")
            records.append(stat)

    print(f"Wrote {len(records)} reads to {args.out}", file=sys.stderr)
    print_report(records)


if __name__ == "__main__":
    main()
