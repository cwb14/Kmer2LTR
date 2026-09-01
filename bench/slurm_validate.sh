#!/bin/bash
#SBATCH -A CHANGEME
#SBATCH -p shared
#SBATCH -n 20
#SBATCH -t 02:00:00
#SBATCH -J ltrk2p_validate
#SBATCH -o bench/out/validate_%j.log
set -euo pipefail
# Set REPO and DATA for your site, or export them before submitting.
REPO=${REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
DATA=${DATA:-$(dirname "$REPO")}
cd "$REPO"
PY=${PY:-python}
export PYTHONPATH="$REPO/src"
D="$DATA"
O=bench/out/configs
mkdir -p "$O"

# Final validation of the SHIPPED defaults -- no flags, so every value
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
$PY -c "
import sys; sys.path.insert(0,'bench'); sys.path.insert(0,'src')
from bench.run_bench import score_gold_grid, cells_to_json
cells_to_json(score_gold_grid('bench/out/gold_truth.tsv', '$O/final_gold.tsv'),
              '$O/cells_gold_FINAL.json')
print('scored final_gold')"

# Row-count invariant: one output row per input record, on every dataset.
bash bench/out/_scratch6/verify_rowcounts.sh "$O"

# TG..CA and TSD -- signals the tool never uses, so they are unbiased external
# accuracy proxies. Secondary evidence: the homology grid is the primary check
# precisely because it needs no proxy at all.
$PY bench/out/_scratch6/tgca_tsd_check.py "$(cat <<'JSON'
[["$DATA"/arab_ltr_all_clean.fa.gz","bench/out/configs/real_arab_ltr_all_clean.tsv","arabidopsis"],
 ["$DATA"/human_ltr_all_clean.fa.gz","bench/out/configs/real_human_ltr_all_clean.tsv","human"],
 ["$DATA"/poa_ltr_all_clean.fa.gz","bench/out/configs/real_poa_ltr_all_clean.tsv","poa"],
 ["$DATA"/MTEC/maizeTE04092026","bench/out/configs/real_maizeTE04092026.tsv","mtec"]]
JSON
)" "$O/tgca_tsd.json"

echo "=== validation complete: $(date) ==="
