#!/bin/bash
#SBATCH -A CHANGEME
#SBATCH -p shared
#SBATCH -n 20
#SBATCH -t 01:00:00
#SBATCH -J ltrk2p_flankbeta
#SBATCH -o bench/out/flankbeta_%j.log
set -euo pipefail
# Set REPO and DATA for your site, or export them before submitting.
REPO=${REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
DATA=${DATA:-$(dirname "$REPO")}
cd "$REPO"
PY=${PY:-python}
export PYTHONPATH="$REPO/src"
O=bench/out/configs

$PY - <<'PYEOF'
import sys, time
sys.path.insert(0, "bench"); sys.path.insert(0, "src")
from bench.run_bench import ablate
from bench.homology_grid import score_homology, cells_to_json
import ltrk2p.align as A

FA, TRUTH = "bench/out/homology_grid.fa", "bench/out/homology_truth.tsv"
DHAT = "bench/out/configs/dhat_hom_candidate.tsv"
for name, beta in (("strict", None), ("b0.25", 0.25), ("b0.5", 0.5),
                   ("b0.75", 0.75), ("b1.0", 1.0), ("b1.5", 1.5)):
    A.FLANK_SENSITIVITY[name] = beta          # register the sweep point
    pred = f"bench/out/configs/pred_hom_beta_{name}.tsv"
    t0 = time.time()
    s = ablate("calibrated", FA, pred, threads=20, flank_sensitivity=name)
    print(f"beta {name}: {s['n']} records in {s['elapsed_s']:.0f}s", flush=True)
    cells_to_json(score_homology(TRUTH, pred, DHAT),
                  f"bench/out/configs/cells_hom_beta_{name}.json")
    print(f"beta {name}: scored", flush=True)
PYEOF
echo "=== beta sweep complete: $(date) ==="
