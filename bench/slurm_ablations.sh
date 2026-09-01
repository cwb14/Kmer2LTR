#!/bin/bash
#SBATCH -A CHANGEME
#SBATCH -p shared
#SBATCH -n 20
#SBATCH -t 02:00:00
#SBATCH -J ltrk2p_ablations
#SBATCH -o bench/out/ablations_%j.log
set -euo pipefail
# Set REPO and DATA for your site, or export them before submitting.
REPO=${REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
DATA=${DATA:-$(dirname "$REPO")}
cd "$REPO"
PY=${PY:-python}
# Ensure workers import this checkout rather than any installed copy.
export PYTHONPATH="$REPO/src"

# One knob at a time from an explicit baseline, over both
# grids -- the 260,876-record gold-perturbed grid (continuity with every prior
# ltrk2p measurement) and the 70,000-record homology grid (perfect elements,
# known mutations/flanks/indels, exact true alignment, no motif or TSD prior).
#
# The t_bits sweep is run here rather than reusing an earlier one: that sweep was
# measured with MAX_EVALUE at 1e-3 while the tool ships 1e-10. Verified directly --
# records the old sweep reports as `pass` are `no_pair` at 1e-10 -- so its numbers
# describe a significance regime the tool does not use.
$PY bench/run_configs.py \
    --outdir bench/out/configs --threads 20 \
    --gold-fa bench/out/gold_perturbed.fa --gold-truth bench/out/gold_truth.tsv \
    --hom-fa bench/out/homology_grid.fa --hom-truth bench/out/homology_truth.tsv \
    --dhat-for baseline \
    --configs inner_joint,comp_core,gaps_static,gaps_adaptive,schedule \
    --tbits-sweep 2,5,8,10,15,20,30 --sweep-base baseline

echo "=== ablations complete: $(date) ==="
