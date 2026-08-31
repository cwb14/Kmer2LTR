#!/bin/bash
#SBATCH -A bio250178
#SBATCH -p shared
#SBATCH -n 20
#SBATCH -t 01:00:00
#SBATCH -J ltrk2p_gold
#SBATCH -o bench/out/gold_robustness_%j.log
set -euo pipefail
PY=/anvil/projects/x-bio250178/conda/envs/ltrrt4/bin/python
cd /anvil/projects/x-bio250178/chris/LTRRT4/ltrk2p

# Mandatory per task-18-brief.md: arabidopsis + truth.fa. Human included too
# (brief: "add human if time allows" -- cheap to include in the same job).
# poa (75,180 records, ~15ms/record) is deliberately skipped: brief calls it
# optional given its size, and arab+truth+human already exercise every
# dataset family (real genomic, library-consensus-constructed, real genomic
# again) the benchmark needs.
$PY bench/gold_robustness.py \
    --dataset "arabidopsis=/anvil/projects/x-bio250178/chris/LTRRT4/arab_ltr_all_clean.fa.gz" \
    --dataset "truth=bench/out/truth.fa" \
    --dataset "human=/anvil/projects/x-bio250178/chris/LTRRT4/human_ltr_all_clean.fa.gz" \
    --outdir bench/out --threads 20 --max-gold 800 --seed 0 -v
