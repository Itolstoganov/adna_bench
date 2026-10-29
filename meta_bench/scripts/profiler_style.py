"""Display names and plot styles of the metagenomic profilers.

The aligners take theirs from shared/scripts/tool_style.py, the same table the
simx benchmark uses. This benchmark profiles a single strobealign version, which
takes strobealign's first color and the name set as `strobealign_name` in
meta_bench/config/config.yaml. Only plots use these: TSVs, directories and rule
names keep the raw profiler names (bwa_aln, strobealign, ...).
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "shared", "scripts"))
import tool_style  # noqa: E402

CONFIG = os.path.join(HERE, "..", "config", "config.yaml")


def _strobealign_name():
    import yaml
    with open(CONFIG) as fh:
        return (yaml.safe_load(fh) or {}).get("strobealign_name", "strobealign")


STROBEALIGN_NAME = _strobealign_name()


def display(profiler):
    """Plot label for a profiler directory/rule name."""
    return STROBEALIGN_NAME if profiler == "strobealign" else tool_style.display(profiler)


def styles(profilers):
    """{profiler: (color, marker)}, keyed by the raw profiler names."""
    by_name = tool_style.styles([display(p) for p in profilers], [STROBEALIGN_NAME])
    return {p: by_name[display(p)] for p in profilers}
