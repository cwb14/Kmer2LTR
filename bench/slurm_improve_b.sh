#!/bin/bash
#SBATCH -A bio250178
#SBATCH -p shared
#SBATCH -n 20
#SBATCH -t 02:00:00
#SBATCH -J ltrk2p_impB
#SBATCH -o bench/out/improve_b_%j.log
set -euo pipefail
cd /anvil/projects/x-bio250178/chris/LTRRT6/ltrk2p
PY=/anvil/projects/x-bio250178/conda/envs/ltrrt4/bin/python
export PYTHONPATH=/anvil/projects/x-bio250178/chris/LTRRT6/ltrk2p/src
BASE=${1:?usage: sbatch bench/slurm_improve_b.sh <config-name>}

# Job B: the t_bits sweep on the pipeline Job A selected, plus the d_hat pass
# that sweep must be binned by. The schedule is then derived by
# bench/derive_schedule.py under the rule committed in 5ac6c89 -- before any of
# these numbers existed.
$PY bench/run_improve.py \
    --outdir bench/out/improve --threads 20 \
    --gold-fa bench/out/gold_perturbed.fa --gold-truth bench/out/gold_truth.tsv \
    --hom-fa bench/out/homology_grid.fa --hom-truth bench/out/homology_truth.tsv \
    --dhat-for "$BASE" \
    --configs "$BASE" \
    --tbits-sweep 2,5,8,10,15,20,30 --sweep-base "$BASE"

echo "=== Job B complete: $(date) ==="
