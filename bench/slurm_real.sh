#!/bin/bash
#SBATCH -A bio250178
#SBATCH -p shared
#SBATCH -n 20
#SBATCH -t 12:00:00
#SBATCH -J ltrk2p_real
#SBATCH -o bench/out/real_%j.log
set -euo pipefail
cd /anvil/projects/x-bio250178/chris/LTRRT4/ltrk2p
PY=/anvil/projects/x-bio250178/conda/envs/ltrrt4/bin/python
D=/anvil/projects/x-bio250178/chris/LTRRT4

# Task 16: run all seven real datasets end to end with the final,
# benchmark-set defaults (T_BITS=10.0, MAX_EVALUE calibrated against
# bench/out/negatives.fa -- see bench/out/memo_real.md). No CLI flags are
# needed: --flank-bits defaults to align.T_BITS and align.classify's
# max_evalue default follows align.MAX_EVALUE, both read at import time.
for f in arab_ltr_all_clean.fa.gz poa_ltr_all_clean.fa.gz human_ltr_all_clean.fa.gz \
         repbase.fa Dfam-RepeatMasker.lib MTEC/maizeTE04092026 riceTElib/rice7.0.0.liban; do
  b=$(basename "$f" | sed 's/\..*//')
  echo "=== $f -> bench/out/real_${b}.tsv ($(date)) ==="
  t0=$SECONDS
  $PY -m ltrk2p "$D/$f" -o "bench/out/real_${b}.tsv" -t 20 -v
  echo "=== $f done in $((SECONDS - t0))s ==="
done
