#!/bin/bash
#SBATCH -A CHANGEME
#SBATCH -p shared
#SBATCH -n 20
#SBATCH -t 00:30:00
#SBATCH -J Kmer2LTR_genome_sweep
#SBATCH -o bench/out/genome_sweep_%j.log
set -euo pipefail
# Set REPO and DATA for your site, or export them before submitting. WORK is the
# directory holding the annotator's `*.work/` dirs; without it every element is
# reported as `unmatched` and the per-source split -- the whole point -- is lost.
REPO=${REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
DATA=${DATA:-$(dirname "$REPO")}
WORK=${WORK:-$DATA}
cd "$REPO"
PY=${PY:-python}
export PYTHONPATH="$REPO/src"
O=bench/out/genome
mkdir -p "$O"

ARAB_G="$DATA/arab_all_genome.fa.gz"
HUMAN_G="$DATA/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"

# 1. --tsd-anchor. This is the only dimension that can move a boundary, so it is
#    the only one that needs Kmer2LTR re-run rather than swept after the fact.
#    0 is the shipped default and the reference point every other value is read
#    against; inf is the hard veto the flag's range is bounded by.
for a in 0 4 8 12 20 1e9; do
  for s in arab human; do
    g=$ARAB_G; [ "$s" = human ] && g=$HUMAN_G
    echo "=== anchor=$a $s ($(date)) ==="
    $PY -m kmer2ltr "$DATA/${s}_ltr_all_clean.fa.gz" -o "$O/${s}_anchor${a}.tsv" \
        -t 20 --genome "$g" --tsd-anchor "$a"
  done
done

# 2. Everything that only reads the answer: TSD length set, shift window,
#    mismatch tolerance, low-complexity floor, and the orientation probe. One
#    Kmer2LTR run each, swept afterwards against shift-budget-matched controls.
$PY bench/genome_sweep.py --elements "$DATA/arab_ltr_all_clean.fa.gz" \
    --genome "$ARAB_G" --tsv "$O/arab_anchor0.tsv" \
    --scn "harvest=$WORK/*.work/*.ltrharvest.stitched.scn" \
          "finder=$WORK/*.work/*.ltrfinder.stitched.scn" \
    --out "$O/arab_sweep.json"

# The human set has no `*.work/` dirs here, so it carries no source split and is
# reported as one population -- read it as a portability check on the parameters
# the arabidopsis split chooses, not as evidence about boundaries.
$PY bench/genome_sweep.py --elements "$DATA/human_ltr_all_clean.fa.gz" \
    --genome "$HUMAN_G" --tsv "$O/human_anchor0.tsv" \
    --out "$O/human_sweep.json"

# 3. Read the anchor grid back per source, against the terminal motif -- the one
#    signal --tsd-anchor cannot see, and therefore the only honest check on it.
for s in arab human; do
  args=(); for a in 0 4 8 12 20 1e9; do args+=("$a=$O/${s}_anchor${a}.tsv"); done
  $PY bench/genome_sweep.py --anchors "${args[@]}" \
      --scn "harvest=$WORK/*.work/*.ltrharvest.stitched.scn" \
            "finder=$WORK/*.work/*.ltrfinder.stitched.scn" > "$O/${s}_anchors.md"
done
$PY bench/genome_sweep.py --report "$O/arab_sweep.json"  > "$O/arab_sweep.md"
$PY bench/genome_sweep.py --report "$O/human_sweep.json" > "$O/human_sweep.md"

echo "=== genome sweep complete: $(date) ==="
