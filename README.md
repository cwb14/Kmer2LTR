# Kmer2LTR

> **This is the legacy branch.** Kmer2LTR has since been rewritten; the
> current version is on [`main`](https://github.com/cwb14/Kmer2LTR). This
> branch is kept so earlier results stay reproducible, and is still the one
> to use if you need a domains file (`-D`), JC69 distances, or the pooled
> `*.summary` output. Both carry the same command name.

Estimate when LTR retrotransposons inserted, from the divergence between their two LTRs.

An LTR-RT starts life with two identical LTRs; after insertion they mutate independently, so
their divergence is a molecular clock - the more they differ, the older the element. Give
Kmer2LTR a FASTA of intact LTR-RTs and it returns, per element, the LTR divergence (p-distance,
JC69, K2P) and the matching age.

## All options

Run `python Kmer2LTR/Kmer2LTR.py -h` for the authoritative list.

```
input / output
  -i, --input-fastas   input LTR-RT FASTA(s)           (required)
  -o                   results table (single input)    (default: <prefix>.LTRs.alns.results)
  -t                   temp dir (single input)         (default: <prefix>_temp; point at fast scratch)

LTR detection & divergence
  -D, --domains        domains TSV(s): name + LTR length - skips the kmer search when you already
                       know the LTR length
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

## Install

Single Python script (3.10+). It calls out to **MAFFT** (alignment) and **trimal** (cleanup), plus
**mmseqs2** for clustering; the density plot needs **numpy** + **matplotlib**.

```bash
mamba create -n kmer2ltr -c bioconda -c conda-forge python=3.11 mafft trimal mmseqs2 numpy matplotlib
mamba activate kmer2ltr
```

The WFA aligner ships prebuilt (`wfa_mac`, `wfa_linux1`, `wfa_linux2`) and is picked automatically;
if none run, compile `wfa.cpp` against [WFA2-lib](https://github.com/smarco/WFA2-lib).

## Quick start

```bash
python Kmer2LTR/Kmer2LTR.py -i elements.fa -u 1.8e-8
```

`-u` is the neutral mutation rate (default `3e-8`; use one appropriate for your species). Outputs
are named from the input prefix (text before the first `.`):

```
elements.LTRs.alns.results          # one row per element (the main table)
elements.LTRs.alns.results.summary  # pooled, genome-wide divergence
elements.LTRs.alns.density.pdf      # divergence / age distribution
```

## Output

`*.LTRs.alns.results` - one element per line, tab-separated, no header:

```
# LTR-RT  LTR_LEN  ALN_LEN  subs  Ti  Tv  p-dist  p-time  JC69-dist  JC69-time  K2P-dist  K2P-time  left_trim  right_trim  end5p  start3p
Gypsy1#LTR_Ty3  584  574  51  40  11  0.088850  1480833  0.094572  1576200  0.096055  1600917  5  5  584  8210
Copia3#LTR_Ty1  301  298  19  15   4  0.063758  1062633  0.066666  1111100  0.067383  1123050  4  6  301  4502
```

- `LTR_LEN` - LTR length used (kmer-discovered, or from the domains file); `ALN_LEN` - length of
  the aligned, scored region.
- `*-dist` / `*-time` - LTR divergence (p-distance, JC69, K2P) and the matching age in years.
- `left_trim` / `right_trim` - bp the WFA boundary detector chopped off each end.
- `end5p` / `start3p` - 1-based coords of the last bp of the 5′ LTR and first bp of the 3′ LTR.

`*.summary` - pooled over all elements (one molecular clock for the whole set):

```
total_length        22408270
total_transitions   532951
total_transversions 217017
raw_d               0.033468
JC69_d              0.034238
K2P_d               0.034368
```

A results file can be fed straight back in as a domains file (`-D`) - column 1 is the name,
column 2 is the LTR length.

## Consensus LTRs & clustering

- `--ltr-consensus` - one IUPAC consensus LTR per element → `*.consensus.fa`.
- `--ltr-cluster` - build the consensus FASTA and cluster it with mmseqs across an identity sweep
  (0.70 → 0.98), one `*.consensus_id<id>_cluster.tsv` per threshold. Low identity lumps elements at
  the lineage level, high identity splits toward families - all are kept.
- `--cluster-only consensus.fa` - just run the sweep on an existing consensus FASTA.

Clustering params were grid-searched on *Arabidopsis*: mmseqs on the consensus LTR beat cd-hit-est,
a wavefront pipeline, and both full-length and 5′-LTR-only inputs. Needs mmseqs; not compatible
with `--wfa-align`.

## Extra outputs

- `--internal-fasta` - write each element's internal (between-LTR) sequence to `*.internal.fa`.
- `--make-perfect-ltr-rt {5p,3p,consensus}` - write "perfect" (unmutated) LTR-RTs: the internal
  sequence flanked by two identical LTR copies, as at insertion. One file per mode
  (`*.perfect_5p.fa`, etc.); headers gain `~LTRlen:<len>`.

## Many files at once

```bash
python Kmer2LTR/Kmer2LTR.py -i species1.fa species2.fa species3.fa -p 50
```

Each file is processed independently; `-o`/`-t` are ignored, outputs are named per input, and the
combined plot is `kmer2ltr_density.pdf`. `-p` is workers **per input file**. Pass one domains file
per input and they're matched by prefix: `-D *.domains`.

## Dropping dubious elements

Both filters are off by default:

- `--max-win-overdisp` - drop elements whose mutations are clumped unevenly along the LTR (often a
  sign the boundaries are wrong). Start at `6`, lower for more stringency.
- `--min-retained-fraction` - drop elements where trimal threw away too much of the LTR. Start at
  `0.5`–`0.6`.

## Benchmarking & developer notes

**Helper scripts** - standalone prep/post-processing, not needed for a normal run:

```bash
# LTR_retriever pass.list  ->  domains file (for -D)
for f in *.pass.list; do python pass_list_domians.py "$f" > "${f%.pass.list}.domains"; done

# pass.list + genome  ->  LTR-RT FASTA
for tsv in *.pass.list; do
  python pass_list_fa_extract.py -fa "${tsv%.pass.list}.fa" -tsv "$tsv" > "${tsv}.fa" &
done; wait

# Kmer2LTR output  ->  RepeatMasker-style library (>NAME#Class/Superfamily)
python kmer2ltrfa_to_RMfa.py --fasta elements.fa --tsv elements.tsv \
    --clean-non-bases --make-perfect-repeat > library.fa
```

**Mutated test library** - `lib_mutator.py` mutates only the LTRs (internal left alone) at a set
rate and Ti/Tv, so you can check recovered divergence against truth. Reads the `~LTRlen:<len>`
headers from `--make-perfect-ltr-rt`:

```bash
python lib_mutator.py -i perfect.fa -o perfect.15mp.titv2.fa -mp 15 -TiTv 2 --seed 51
```

**Runtime** (100 workers): ~20,700 LTR-RTs in **25m**; **~3.5m** with a domains file.

**Debugging** - `--debug` keeps every intermediate (kmer pairs, dot plots, alignments, WFA
boundary calls) under `<temp>/<seq>/debug/`, with a per-element narrative.

**Legacy** - `extract_ltrs.py`, `filter_kmers.py`, and `map_kmers_to_fasta.py` are the original
standalone kmer/extraction steps, now folded into `Kmer2LTR.py` (in-process, no jellyfish). Kept
for reference only.
