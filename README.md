# Kmer2LTR

Estimate when LTR retrotransposons inserted, from the divergence between their two LTRs.

An LTR-RT is born with two identical LTRs (the long terminal repeats at each end).
After insertion the two copies mutate independently, so the divergence between them is a
molecular clock: the more they differ, the older the element. Kmer2LTR finds the two LTRs
in each element, aligns them, measures their divergence, and converts it to an age with a
neutral mutation rate.

Give it a FASTA of intact LTR-RTs. It returns, per element, the LTR divergence (p-distance,
JC69, K2P) and the matching age estimate.

## How it works

For each sequence in the input FASTA:

1. **Find the two LTRs.**
   - *Fast path* (with `-D`/`--domains`): you already know the LTR length, so it just cuts
     the first and last `LTR_len + extension` bp.
   - *kmer path* (default): find kmers (length 8–12) that occur in both halves of the
     element, ≥ `-d` bp apart; keep the collinear set; the kmer spread marks the 5′ and 3′
     LTR boundaries, symmetrized and padded by `-e` bp.
2. **Align** the 5′ LTR against the 3′ LTR (MAFFT by default, or WFA with `--wfa-align`).
3. **Trim** the alignment with trimal, then re-detect the reliable core with WFA (`-k 5 -K`:
   needs 5 matching bp to enter/leave the trustworthy region, and drops those anchor columns
   so seeded matches don't inflate identity). This is how ragged, uncertain LTR termini get
   chopped before counting mutations.
4. **Score** the divergence (substitutions, transitions/transversions) and convert to age:
   `age = distance / (2 × mutation_rate)`.

## Install

Kmer2LTR is a single Python script. It shells out to a few standard tools:

| Tool | Needed for | Install |
|------|-----------|---------|
| **MAFFT** | alignment (default) | `mamba install -c bioconda mafft` |
| **trimal** | alignment cleanup (always) | `mamba install -c bioconda trimal` |
| **mmseqs2** | clustering only (`--ltr-cluster`, `--cluster-only`) | `mamba install -c bioconda mmseqs2` |
| **numpy + matplotlib** | the density plot (default; skip with `--no-plot`) | `mamba install numpy matplotlib` |
| **WFA aligner** | divergence scoring (always) | **ships prebuilt** — see below |

Python 3.10+.

```bash
mamba create -n kmer2ltr -c bioconda -c conda-forge python=3.11 mafft trimal mmseqs2 numpy matplotlib
mamba activate kmer2ltr
```

The WFA aligner is bundled as prebuilt binaries (`wfa_mac`, `wfa_linux1`, `wfa_linux2`) sitting
next to `Kmer2LTR.py`; the right one is picked automatically for your platform. If none run on
your machine you'll need to compile `wfa.cpp` against [WFA2-lib](https://github.com/smarco/WFA2-lib)
(it's not a plain `g++ wfa.cpp` — it links `libwfa2.a`). The shipped binaries cover macOS and Linux,
so most users never touch this.

## Quick start

```bash
python Kmer2LTR/Kmer2LTR.py -i elements.fa -u 1.8e-8
```

`-u` is the neutral mutation rate (default `3e-8`; use one appropriate for your species).

Writes, using the input prefix (text before the first `.`):

```
elements.LTRs.alns.results          # one row per element (the main table)
elements.LTRs.alns.results.summary  # pooled, genome-wide divergence
elements.LTRs.alns.density.pdf      # divergence / age distribution
elements_temp/                      # scratch (deleted unless -k)
```

## Output

`*.LTRs.alns.results` — one element per line, tab-separated, no header:

```
# LTR-RT  LTR_LEN  ALN_LEN  subs  Ti  Tv  p-dist  p-time  JC69-dist  JC69-time  K2P-dist  K2P-time  left_trim  right_trim  end5p  start3p
Gypsy1#LTR_Ty3  584  574  51  40  11  0.088850  1480833  0.094572  1576200  0.096055  1600917  5  5  584  8210
Copia3#LTR_Ty1  301  298  19  15   4  0.063758  1062633  0.066666  1111100  0.067383  1123050  4  6  301  4502
```

- `LTR_LEN` — LTR length used (kmer-discovered, or from the domains file).
- `ALN_LEN` — length of the aligned, scored region. `LTR_LEN` can be larger (gaps in the LTR)
  or smaller (the extension expanded the kmer-called boundary).
- `subs / Ti / Tv` — substitutions, transitions, transversions.
- `*-dist` — LTR divergence under p-distance, JC69, K2P. `*-time` — the matching age in years.
- `left_trim / right_trim` — bp the WFA boundary detector chopped off each end.
- `end5p / start3p` — 1-based coordinates of the last bp of the 5′ LTR and the first bp of the
  3′ LTR (extension excluded). These mark where the internal region sits.

`*.summary` — pooled over all elements (one molecular clock for the whole set):

```
total_length        22408270
total_transitions   532951
total_transversions 217017
raw_d               0.033468
JC69_d              0.034238
K2P_d               0.034368
```

A results file can be fed straight back in as a domains file (`-D`) — column 1 is the name,
column 2 is the LTR length.

## Fast path: you already know the LTR lengths

If you have the LTR lengths (e.g. from LTR_retriever), skip the kmer search with a **domains
file** — two tab-separated columns, element name and LTR length:

```bash
head -1 species1.domains
CMHA_chr1:90368..96317   172     # this element's LTRs are each 172 bp

python Kmer2LTR/Kmer2LTR.py -i species1.fa -D species1.domains
```

Much faster, and it removes any uncertainty in the LTR boundaries.

## Many files at once

```bash
python Kmer2LTR/Kmer2LTR.py -i species1.fa species2.fa species3.fa -p 50
```

With multiple inputs each file is processed independently and `-o`/`-t` are ignored — outputs
are named per input (`species1.LTRs.alns.results`, …) and the combined plot is
`kmer2ltr_density.pdf`. Pass one domains file per input; they're matched by prefix:

```bash
python Kmer2LTR/Kmer2LTR.py -i species1.fa species2.fa species3.fa -p 50 -D *.domains
```

`-p` is workers **per input file**.

## Dropping dubious elements

Both filters are **off by default** — turn them on if you want to be strict:

- `--max-win-overdisp` — drop elements whose mutations are clumped unevenly along the LTR
  (e.g. one half clean, the other half divergent), which usually means the boundaries are
  wrong. Start at `6` and lower for more stringency.
- `--min-retained-fraction` — drop elements where trimal threw away too much of the LTR.
  Start at `0.5`–`0.6`.

## Consensus LTRs & clustering

- `--ltr-consensus` — write one IUPAC consensus LTR per element to `*.consensus.fa`.
- `--ltr-cluster` — build the consensus FASTA **and** cluster it with mmseqs across a sweep of
  identity thresholds (0.70 → 0.98), writing one `*.consensus_id<id>_cluster.tsv` per threshold.
  Low identity lumps elements at the lineage level, high identity splits toward families/recent
  bursts — different analyses want different depths, so all are kept.
- `--cluster-only consensus.fa` — just run the clustering sweep on an existing consensus FASTA.

The clustering parameters (`mmseqs easy-cluster`, `-c 0.5 -s 7.5`, etc.) were grid-searched on
*Arabidopsis* LTR annotations to jointly minimize singletons and cross-family mixing. There,
mmseqs on the **consensus LTR** beat cd-hit-est, a wavefront pipeline, and both full-length and
5′-LTR-only inputs. Needs mmseqs in `PATH`; not compatible with `--wfa-align`.

## Extra outputs

- `--internal-fasta` — also write each element's internal (between-LTR) sequence to
  `*.internal.fa`. When the header carries a `chrom:start-end` locus, the internal record gets
  the internal region's genomic interval too.
- `--make-perfect-ltr-rt {5p,3p,consensus}` — write "perfect" (unmutated) LTR-RTs: the internal
  sequence flanked by two identical LTR copies, as at insertion. One output per mode
  (`*.perfect_5p.fa`, `*.perfect_3p.fa`, `*.perfect_consensus.fa`). Headers gain `~LTRlen:<len>`.

## All options

Run `python Kmer2LTR/Kmer2LTR.py -h` for the authoritative list. The essentials:

```
input / output
  -i, --input-fastas   input LTR-RT FASTA(s)           (required)
  -o                   results table (single input)    (default: <prefix>.LTRs.alns.results)
  -t                   temp dir (single input)         (default: <prefix>_temp; point at fast scratch)

LTR detection & divergence
  -D, --domains        domains TSV(s): name + LTR length (fast path)
  -u                   mutation rate for age           (default: 3e-8)
  -e                   bp kept past each LTR end        (default: 120)
  --wfa-align          align with WFA, not MAFFT (~30-50x faster; slightly different divergence)

quality filters (off by default)
  --max-win-overdisp   drop elements with clumped mutations   (try 6)
  --min-retained-fraction  require this fraction to survive trimming  (try 0.6)

consensus & clustering
  --ltr-consensus      write IUPAC consensus LTRs
  --ltr-cluster        consensus + mmseqs identity sweep (0.70-0.98)
  --cluster-only FA    only cluster an existing consensus FASTA
  --internal-fasta     write internal (between-LTR) sequences
  --make-perfect-ltr-rt {5p,3p,consensus}  write unmutated LTR-RTs

kmer boundary tuning (advanced)
  --kmer-range MIN MAX kmer lengths for boundary anchoring   (default: 8 12)
  -d                   min bp between the two kmer copies     (default: 80)

performance & temp
  -p                   workers per input                (default: 20)
  -k                   keep the temp directory
  --purge-subdirs [N]  delete per-element temp as you go (helps with huge inputs)
  --reuse-existing     resume: keep existing results, only process what's missing
  --assume-duplicate-same-ltr   fast-path shortcut for duplicate headers (use with care)

diagnostics
  -v                   echo each command
  --debug              full per-element intermediates under <temp>/<seq>/debug/ (slow)
  --no-plot            skip the density plot
```

## Helper scripts

Standalone utilities for prepping inputs and post-processing — none are needed for a normal run:

```bash
# LTR_retriever pass.list  ->  domains file (for -D)
for f in *.pass.list; do
  python pass_list_domians.py "$f" > "${f%.pass.list}.domains"
done

# LTR_retriever pass.list + genome  ->  LTR-RT FASTA
for tsv in *.pass.list; do
  python pass_list_fa_extract.py -fa "${tsv%.pass.list}.fa" -tsv "$tsv" > "${tsv}.fa" &
done
wait

# Kmer2LTR output  ->  RepeatMasker-style library (>NAME#Class/Superfamily)
python kmer2ltrfa_to_RMfa.py --fasta elements.fa --tsv elements.tsv \
    --clean-non-bases --make-perfect-repeat > library.fa
```

## Benchmarking & developer notes

**Make a mutated test library** — `lib_mutator.py` mutates only the LTRs (internal sequence left
alone) at a set rate and Ti/Tv, so you can check recovered divergence against truth. It reads the
`~LTRlen:<len>` headers that `--make-perfect-ltr-rt` writes:

```bash
python lib_mutator.py -i perfect.fa -o perfect.15mp.titv2.fa -mp 15 -TiTv 2 --seed 51
```

**Runtime** (100 workers): ~20,700 LTR-RTs in **25m**. With a domains file (fast path):
**3.5m**.

**Debugging** — `--debug` keeps every intermediate (kmer pairs, dot plots, alignments, WFA
boundary calls) under `<temp>/<seq>/debug/`, with a per-element narrative. Use it to see why a
domains run disagrees with a kmer run, why MAFFT and WFA differ, or why an element was skipped.

**Legacy** — `extract_ltrs.py`, `filter_kmers.py`, and `map_kmers_to_fasta.py` are the original
standalone kmer/extraction steps. That logic now lives inside `Kmer2LTR.py` (in-process, no
jellyfish), so these are kept for reference only — they aren't part of the pipeline.
