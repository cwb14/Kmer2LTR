#!/bin/bash
#SBATCH -A bio250178
#SBATCH -p shared
#SBATCH -n 20
#SBATCH -t 08:00:00
#SBATCH -J ltrk2p_bench
#SBATCH -o bench/out/bench_%j.log
set -euo pipefail
cd /anvil/projects/x-bio250178/chris/LTRRT4/ltrk2p
PY=/anvil/projects/x-bio250178/conda/envs/ltrrt4/bin/python

# --- Task 14 brief, literal: spec-6.4 ablation grid on a synthetic
# make_dataset() grid built from the library-consensus truth set. Kept
# modest (--n-truth-sample) -- the headline real-data numbers this task
# exists to produce come from the gold-perturbed grid below, not from here.
$PY bench/run_bench.py --truth-fa bench/out/truth.fa --truth-tsv bench/out/truth.tsv \
    --outdir bench/out --threads 20 --all-ablations --n-truth-sample 1500 --seed 0

# --- Priority-change deliverable: T_BITS sweep + ablation grid measured on
# Task 18's REAL gold-perturbed data (bench/out/gold_perturbed.fa /
# gold_truth.tsv), the d_hat pass the divergence-aware-threshold analysis
# needs, and the Stage 4 real-data diff (with vs without stage4 on the raw,
# unperturbed real datasets). t_bits=5.0 and the calibrated/trim_0 configs
# reuse bench/out/gold_pred.tsv (Task 18's existing default-config run)
# instead of recomputing it.
$PY bench/run_bench.py \
    --outdir bench/out --threads 20 \
    --gold-fa bench/out/gold_perturbed.fa --gold-truth bench/out/gold_truth.tsv \
    --tbits-sweep 2,5,8,10,15,20,30 \
    --gold-ablations --dhat \
    --stage4-diff \
    --raw-dataset "arabidopsis=/anvil/projects/x-bio250178/chris/LTRRT4/arab_ltr_all_clean.fa.gz" \
    --raw-dataset "human=/anvil/projects/x-bio250178/chris/LTRRT4/human_ltr_all_clean.fa.gz" \
    --raw-dataset "truth=bench/out/truth.fa"
