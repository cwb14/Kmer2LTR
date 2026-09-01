#!/bin/bash
#SBATCH -A CHANGEME
#SBATCH -p shared
#SBATCH -n 20
#SBATCH -t 01:00:00
#SBATCH -J ltrk2p_gold
#SBATCH -o bench/out/gold_robustness_%j.log
set -euo pipefail
PY=${PY:-python}
cd "$DATA"/ltrk2p

# arabidopsis + truth.fa + human. poa (75,180 records, ~15 ms/record) is skipped
# deliberately: the other three already exercise both dataset families the
# benchmark needs (real genomic, and library-consensus construction), and poa
# would dominate the runtime without adding a family.
$PY bench/gold_robustness.py \
    --dataset "arabidopsis="$DATA"/arab_ltr_all_clean.fa.gz" \
    --dataset "truth=bench/out/truth.fa" \
    --dataset "human="$DATA"/human_ltr_all_clean.fa.gz" \
    --outdir bench/out --threads 20 --max-gold 800 --seed 0 -v
