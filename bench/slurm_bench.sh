#!/bin/bash
#SBATCH -A CHANGEME
#SBATCH -p shared
#SBATCH -n 20
#SBATCH -t 08:00:00
#SBATCH -J Kmer2LTR_bench
#SBATCH -o bench/out/bench_%j.log
set -euo pipefail
# Set REPO and DATA for your site, or export them before submitting.
REPO=${REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
DATA=${DATA:-$(dirname "$REPO")}
cd "$REPO"
PY=${PY:-python}

# --- spec 6.4 ablation grid on a synthetic
# make_dataset() grid built from the library-consensus truth set. Kept
# modest (--n-truth-sample) -- the headline real-data numbers this task
# exists to produce come from the gold-perturbed grid below, not from here.
$PY bench/run_bench.py --truth-fa bench/out/truth.fa --truth-tsv bench/out/truth.tsv \
    --outdir bench/out --threads 20 --all-ablations --n-truth-sample 1500 --seed 0

# --- T_BITS sweep + ablation grid measured on the REAL gold-perturbed
# data built by gold_robustness.py (bench/out/gold_perturbed.fa /
# gold_truth.tsv), the d_hat pass the divergence-aware-threshold analysis
# needs, and the Stage 4 real-data diff (with vs without stage4 on the raw,
# unperturbed real datasets). t_bits=5.0 and the calibrated/trim_0 configs
# reuse bench/out/gold_pred.tsv (the existing default-config run)
# instead of recomputing it.
$PY bench/run_bench.py \
    --outdir bench/out --threads 20 \
    --gold-fa bench/out/gold_perturbed.fa --gold-truth bench/out/gold_truth.tsv \
    --tbits-sweep 2,5,8,10,15,20,30 \
    --gold-ablations --dhat \
    --stage4-diff \
    --raw-dataset "arabidopsis="$DATA"/arab_ltr_all_clean.fa.gz" \
    --raw-dataset "human="$DATA"/human_ltr_all_clean.fa.gz" \
    --raw-dataset "truth=bench/out/truth.fa"
