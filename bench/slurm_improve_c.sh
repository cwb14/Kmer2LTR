#!/bin/bash
#SBATCH -A bio250178
#SBATCH -p shared
#SBATCH -n 20
#SBATCH -t 02:00:00
#SBATCH -J ltrk2p_impC
#SBATCH -o bench/out/improve_c_%j.log
set -euo pipefail
cd /anvil/projects/x-bio250178/chris/LTRRT6/ltrk2p
PY=/anvil/projects/x-bio250178/conda/envs/ltrrt4/bin/python
export PYTHONPATH=/anvil/projects/x-bio250178/chris/LTRRT6/ltrk2p/src
D=/anvil/projects/x-bio250178/chris/LTRRT6
O=bench/out/improve
mkdir -p "$O"

# Job C: final validation of the SHIPPED defaults -- no flags, so every value
# comes from src/ltrk2p/align.py exactly as a user would get it.
for f in arab_ltr_all_clean.fa.gz poa_ltr_all_clean.fa.gz human_ltr_all_clean.fa.gz \
         repbase.fa Dfam-RepeatMasker.lib MTEC/maizeTE04092026 riceTElib/rice7.0.0.liban; do
  b=$(basename "$f" | sed 's/\..*//')
  echo "=== $f -> $O/real_${b}.tsv ($(date)) ==="
  t0=$SECONDS
  $PY -m ltrk2p "$D/$f" -o "$O/real_${b}.tsv" -t 20 -v
  echo "=== $f done in $((SECONDS - t0))s ==="
done

# Negative controls: real non-LTR TEs, and a composition-matched shuffled null.
# The significance model is pinned to GENERIC_MATRIX + SIG_GAPS, so MAX_EVALUE
# should be untouched by this campaign -- "should be" is what this checks.
for n in negatives shuffled_truth; do
  echo "=== $n ($(date)) ==="
  $PY -m ltrk2p "bench/out/${n}.fa" -o "$O/neg_${n}.tsv" -t 20 -v
done

# Final runs of both benchmark grids at the shipped defaults.
$PY -m ltrk2p bench/out/homology_grid.fa -o "$O/final_homology.tsv" -t 20 -v
$PY bench/homology_grid.py score --truth bench/out/homology_truth.tsv \
    --pred "$O/final_homology.tsv" --out "$O/cells_hom_FINAL.json"
$PY -m ltrk2p bench/out/gold_perturbed.fa -o "$O/final_gold.tsv" -t 20 -v

echo "=== Job C complete: $(date) ==="
