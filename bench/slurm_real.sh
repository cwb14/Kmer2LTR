#!/bin/bash
#SBATCH -A CHANGEME
#SBATCH -p shared
#SBATCH -n 20
#SBATCH -t 12:00:00
#SBATCH -J ltrk2p_real
#SBATCH -o bench/out/real_%j.log
set -euo pipefail
# Set REPO and DATA for your site, or export them before submitting.
REPO=${REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
DATA=${DATA:-$(dirname "$REPO")}
cd "$REPO"
PY=${PY:-python}
D="$DATA"

# Run all seven real datasets end to end with the shipped defaults. No CLI flags
# are needed: every default is read from src/ltrk2p/align.py at import time, so
# this measures exactly what a user gets.
for f in arab_ltr_all_clean.fa.gz poa_ltr_all_clean.fa.gz human_ltr_all_clean.fa.gz \
         repbase.fa Dfam-RepeatMasker.lib MTEC/maizeTE04092026 riceTElib/rice7.0.0.liban; do
  b=$(basename "$f" | sed 's/\..*//')
  echo "=== $f -> bench/out/real_${b}.tsv ($(date)) ==="
  t0=$SECONDS
  $PY -m ltrk2p "$D/$f" -o "bench/out/real_${b}.tsv" -t 20 -v
  echo "=== $f done in $((SECONDS - t0))s ==="
done
