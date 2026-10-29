"""Display names, colors and markers shared by the simx and metagenomic plots.

Aligners have a fixed name and style. Strobealign is the tool under test, so it
can appear as several versions: those keep the names they are given (the simx
config's `versions[].name`, the metagenomic config's `strobealign_name`) and take
their style from STROBEALIGN_STYLES by position. The first two strobealign
colors are kept clear of the aligners' colors.

Names are matched ignoring case, spaces, '-' and '_', so a table's raw program
name ("bwaaln", "bwa mem", "bwa_aln") and its display name ("BWA-aln") resolve
to the same tool. A " (k=N)" suffix is kept on display and ignored for style.
"""

import re

# canonical key -> (display name, color, marker)
ALIGNERS = {
    "bowtie2": ("Bowtie2", "#1f77b4", "o"),
    "bwaaln": ("BWA-aln", "#ff7f0e", "s"),
    "bwamem": ("BWA-mem", "#2ca02c", "^"),
    "malt": ("MALT", "#d62728", "D"),
}

# (color, marker), one per strobealign version in config order
STROBEALIGN_STYLES = [
    ("#9467bd", "v"),   # purple
    ("#17becf", "P"),   # cyan
    ("#bcbd22", "X"),   # olive
    ("#8c564b", "*"),   # brown
    ("#e377c2", "h"),   # pink
]
MAX_STROBEALIGN_VERSIONS = len(STROBEALIGN_STYLES)

# for tools that are neither, in order of appearance
FALLBACK_STYLES = [
    ("#7f7f7f", "<"),
    ("#393b79", ">"),
    ("#637939", "p"),
    ("#843c39", "8"),
]

_K_SUFFIX = re.compile(r" \(k=\d+\)$")


def _base(name):
    return _K_SUFFIX.sub("", str(name))


def _key(name):
    return re.sub(r"[\s_-]", "", _base(name)).lower()


def display(name):
    """Display name for `name`: an aligner's fixed name (keeping a " (k=N)"
    suffix), anything else unchanged."""
    name = str(name)
    entry = ALIGNERS.get(_key(name))
    if entry is None:
        return name
    return entry[0] + name[len(_base(name)):]


def styles(names, strobealign=()):
    """{name: (color, marker)} for `names`.

    `strobealign` lists the strobealign version names in config order; each takes
    the STROBEALIGN_STYLES entry at its position, so a version keeps its color
    whichever subset of tools a plot shows.
    """
    strobealign = list(dict.fromkeys(_base(s) for s in strobealign))
    if len(strobealign) > MAX_STROBEALIGN_VERSIONS:
        raise ValueError(f"at most {MAX_STROBEALIGN_VERSIONS} strobealign versions "
                         f"have a color, got {len(strobealign)}: {strobealign}")
    out, extra = {}, 0
    for n in dict.fromkeys(names):
        entry = ALIGNERS.get(_key(n))
        if entry is not None:
            out[n] = entry[1:]
        elif _base(n) in strobealign:
            out[n] = STROBEALIGN_STYLES[strobealign.index(_base(n))]
        else:
            out[n] = FALLBACK_STYLES[extra % len(FALLBACK_STYLES)]
            extra += 1
    return out
