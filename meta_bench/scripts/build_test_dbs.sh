#!/usr/bin/env bash
# Build the small KrakenUniq + MALT databases for `run_full_benchmark.sh --test`
# from the test genomes, in place of the ~750 GB downloaded ones. Run inside the
# aMeta conda env (krakenuniq-build, jellyfish 1.x).
#
#   build_test_dbs.sh GENOMES_DIR ACCESSIONS_TSV KRAKENUNIQ_DB MALT_FASTA \
#                     MALT_SEQID2TAXID MALT_ACC2TAXID [THREADS]
#
# The output paths are the ones the test aMeta config names. Every genome listed
# in ACCESSIONS_TSV (microbes and host) goes into both databases, each sequence
# mapped to that row's taxid. The NCBI taxonomy (~60 MB download) is shared by
# KrakenUniq, and through its taxDB, by ngsLCA and the summary. It is pinned to
# the June 2022 archive, the age of the full benchmark's KrakenUniq DB: ngsLCA's
# built-in -fix_ncbi patches assume taxids that newer dumps no longer have
# (assertion `parent.count(out[i])==1` in mod_db).
set -euo pipefail

[[ $# -ge 6 ]] || { sed -n '2,15p' "$0" >&2; exit 1; }
genomes="$1" accessions="$2" ku_db="$3" malt_fna="$4" seqid2taxid="$5" acc2taxid="$6"
threads="${7:-4}"
TAXDUMP_URL="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2022-06-01.zip"

# species -> taxid, from the accessions table
mapfile -t rows < <(awk -F'\t' '!/^#/ && $1 != "species" && NF >= 3 { print $1 "\t" $3 }' "$accessions")
(( ${#rows[@]} )) || { echo "error: no genomes in $accessions" >&2; exit 1; }

# 1. MALT reference: all genomes in one FASTA + seqid/accession -> taxid maps.
if [[ -s "$malt_fna" && -s "$seqid2taxid" && -s "$acc2taxid" ]]; then
  echo ">> MALT test reference present, skipping" >&2
else
  mkdir -p "$(dirname "$malt_fna")" "$(dirname "$seqid2taxid")" "$(dirname "$acc2taxid")"
  : > "$malt_fna.tmp"; : > "$seqid2taxid.tmp"
  printf 'accession\taccession.version\ttaxid\tgi\n' > "$acc2taxid.tmp"
  for row in "${rows[@]}"; do
    species="${row%%$'\t'*}" taxid="${row#*$'\t'}"
    fa="$genomes/$species.fa"
    [[ -s "$fa" ]] || { echo "error: missing genome $fa (run the test simulation's downloads first)" >&2; exit 1; }
    [[ "$taxid" =~ ^[0-9]+$ ]] || { echo "error: no taxid for $species in $accessions" >&2; exit 1; }
    cat "$fa" >> "$malt_fna.tmp"
    # NCBI accession2taxid layout (accession without version first), as the
    # full nucl_gb.accession2taxid, so Malt_Acc2Taxa and ngsLCA read it the same way.
    grep '^>' "$fa" | sed -E 's/^>([^[:space:]]+).*/\1/' | while read -r id; do
      printf '%s\t%s\n' "$id" "$taxid" >> "$seqid2taxid.tmp"
      printf '%s\t%s\t%s\t0\n' "${id%.*}" "$id" "$taxid" >> "$acc2taxid.tmp"
    done
  done
  mv "$malt_fna.tmp" "$malt_fna"; mv "$seqid2taxid.tmp" "$seqid2taxid"; mv "$acc2taxid.tmp" "$acc2taxid"
  echo ">> wrote $malt_fna ($(grep -c '^>' "$malt_fna") sequences)" >&2
fi

# 2. KrakenUniq: library = the MALT FASTA, the same seqid -> taxid map, NCBI taxonomy.
if [[ -s "$ku_db/database.kdb" && -s "$ku_db/database.idx" && -s "$ku_db/taxDB" ]]; then
  echo ">> KrakenUniq test database present, skipping" >&2
  exit 0
fi
mkdir -p "$ku_db/library" "$ku_db/taxonomy"
if [[ ! -s "$ku_db/taxonomy/nodes.dmp" || ! -s "$ku_db/taxonomy/names.dmp" ]]; then
  echo ">> downloading the NCBI taxonomy (2022-06-01 archive)" >&2
  curl -SL --fail --retry 5 -o "$ku_db/taxonomy/taxdmp.zip" "$TAXDUMP_URL"
  python -c 'import sys, zipfile; zipfile.ZipFile(sys.argv[1]).extractall(sys.argv[2], ["nodes.dmp", "names.dmp"])' \
    "$ku_db/taxonomy/taxdmp.zip" "$ku_db/taxonomy"
  rm "$ku_db/taxonomy/taxdmp.zip"
fi
ln -sfn "$(readlink -f "$malt_fna")" "$ku_db/library/library.fna"
cp "$seqid2taxid" "$ku_db/seqid2taxid.map"
echo ">> building the KrakenUniq test database in $ku_db" >&2
# Minimizer length 11 as in aMeta's .test: the default 15 makes an 8 GB
# database.idx regardless of the library size (4^15 buckets x 8 bytes).
krakenuniq-build --db "$ku_db" --threads "$threads" --minimizer-len 11 \
  --jellyfish-bin "$(command -v jellyfish)"
# Drop build intermediates and the raw taxonomy dump; queries need only
# database.kdb/.idx and taxDB.
rm -f "$ku_db/database0.kdb" "$ku_db/database.jdb" "$ku_db/taxonomy/nodes.dmp" "$ku_db/taxonomy/names.dmp"
