#!/bin/bash
#SBATCH -A CHANGEME
#SBATCH -p shared
#SBATCH -n 20
#SBATCH -t 02:00:00
#SBATCH -J ltrk2p_tsweep
#SBATCH -o bench/out/tsweep_%j.log
set -euo pipefail
# Set REPO and DATA for your site, or export them before submitting.
REPO=${REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
DATA=${DATA:-$(dirname "$REPO")}
cd "$REPO"
PY=${PY:-python}
export PYTHONPATH="$REPO/src"
BASE=${1:?usage: sbatch bench/slurm_threshold_sweep.sh <config-name>}

# The t_bits sweep on a chosen pipeline, plus the d_hat pass that sweep must be
# binned by. The schedule is then derived by bench/calibrate_flank_threshold.py,
# whose rule is fixed before these numbers exist.
$PY bench/run_configs.py \
    --outdir bench/out/configs --threads 20 \
    --gold-fa bench/out/gold_perturbed.fa --gold-truth bench/out/gold_truth.tsv \
    --hom-fa bench/out/homology_grid.fa --hom-truth bench/out/homology_truth.tsv \
    --dhat-for "$BASE" \
    --configs "$BASE" \
    --tbits-sweep 2,5,8,10,15,20 --sweep-base "$BASE"

echo "=== sweep complete: $(date) ==="
