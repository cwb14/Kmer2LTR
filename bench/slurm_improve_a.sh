#!/bin/bash
#SBATCH -A bio250178
#SBATCH -p shared
#SBATCH -n 20
#SBATCH -t 02:00:00
#SBATCH -J ltrk2p_impA
#SBATCH -o bench/out/improve_a_%j.log
set -euo pipefail
cd /anvil/projects/x-bio250178/chris/LTRRT6/ltrk2p
PY=/anvil/projects/x-bio250178/conda/envs/ltrrt4/bin/python
# The conda env's editable install of ltrk2p points at a DIFFERENT checkout
# (LTRRT4). Without this every worker would silently benchmark that tree.
export PYTHONPATH=/anvil/projects/x-bio250178/chris/LTRRT6/ltrk2p/src

# Job A: one knob at a time, from the shipped pre-campaign baseline, over both
# grids -- the 260,876-record gold-perturbed grid (continuity with every prior
# ltrk2p measurement) and the 70,000-record homology grid (perfect elements,
# known mutations/flanks/indels, exact true alignment, no motif or TSD prior).
$PY bench/run_improve.py \
    --outdir bench/out/improve --threads 20 \
    --gold-fa bench/out/gold_perturbed.fa --gold-truth bench/out/gold_truth.tsv \
    --hom-fa bench/out/homology_grid.fa --hom-truth bench/out/homology_truth.tsv \
    --configs baseline,stage4_rerun,keep_weak,graded,inner_joint,comp_core,gaps_static,gaps_adaptive,schedule

echo "=== Job A complete: $(date) ==="
