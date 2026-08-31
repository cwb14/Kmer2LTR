#!/bin/bash
#SBATCH -A bio250178
#SBATCH -p shared
#SBATCH -n 20
#SBATCH -t 02:00:00
#SBATCH -J ltrk2p_grad
#SBATCH -o bench/out/improve_graded_%j.log
set -euo pipefail
cd /anvil/projects/x-bio250178/chris/LTRRT6/ltrk2p
PY=/anvil/projects/x-bio250178/conda/envs/ltrrt4/bin/python
export PYTHONPATH=/anvil/projects/x-bio250178/chris/LTRRT6/ltrk2p/src
$PY bench/run_improve.py \
    --outdir bench/out/improve --threads 20 \
    --gold-fa bench/out/gold_perturbed.fa --gold-truth bench/out/gold_truth.tsv \
    --hom-fa bench/out/homology_grid.fa --hom-truth bench/out/homology_truth.tsv \
    --dhat-for baseline --configs graded
echo "=== graded complete: $(date) ==="
