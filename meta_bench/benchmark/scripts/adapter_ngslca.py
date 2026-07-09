#!/usr/bin/env python3
"""Adapter STUB: aligner + ngsLCA -> normalized taxon-abundance table.

Reserved for the future aMeta-ngslca profiler (replacing the MALT branch with a
read aligner + ngsLCA, https://github.com/miwipe/ngsLCA). Not implemented.

To add it: emit the same normalized TSV contract as adapter_malt.py
    sample <TAB> taxid <TAB> taxon <TAB> count
from ngsLCA's per-read LCA assignments (ngsLCA reports taxids directly, so this
adapter should populate the taxid column), then register "aMeta-ngslca" in
run.sh's profiler dispatch. No changes to score.py are required.
"""

import sys


def main() -> None:
    sys.exit(
        "adapter_ngslca is not implemented yet. The aMeta-ngslca profiler "
        "(aligner + ngsLCA) is reserved; implement this adapter to emit the "
        "normalized sample/taxid/taxon/count TSV. See the module docstring."
    )


if __name__ == "__main__":
    main()
