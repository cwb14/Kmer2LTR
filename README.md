# Kmer2LTR

Find the two LTRs in an LTR retrotransposon and report how far they have
diverged.

Give it a FASTA of putative intact elements; it gives you back a TSV with one
row per input record: where each LTR starts and ends, the Kimura 2-parameter
distance between them, and the alignment that produced it.

```bash
Kmer2LTR elements.fa.gz -o elements.tsv -t 20
```

## What it assumes

One thing only: **the 5' and 3' ends of each input sequence are homologous to
each other.** No LTR length parameter, no terminal-motif requirement, no family
model, no strand handling (a direct terminal repeat stays one under
reverse-complement). Everything else is estimated from the element itself.

That matters in practice: because the tool never looks for `TG`…`CA` termini or
target-site duplications, those signals stay available as *independent* checks on
its output. They are used that way in `docs/benchmarks.md` and nowhere in the
code. The `motif` column reports whichever two dinucleotides the boundary landed
on — it is read off the answer, never used to find it, which is exactly what
makes it worth reading.

## What it is for

Two things, and the second comes free with the first.

1. **Dating insertions.** An LTR retrotransposon inserts with two identical
   LTRs, which then diverge as neutral sequence. Their divergence is a molecular
   clock for the insertion, and K2P is the standard correction for multiple hits
   at the same site.

2. **Detecting overextension.** Structure-based annotators often hand you a bit
   of flanking genome along with the element. Since the coordinates are relative
   to the sequence you supplied, a perfectly bounded element reports
   `ltr5_start = 1` and `ltr3_end = seq_len`; anything else is flank, quantified
   directly by the `flank5_len` and `flank3_len` columns.

## Install

```bash
mamba env create -f environment.yml
mamba activate kmer2ltr
pip install -e .
```

Dependencies: `parasail` (SIMD Smith-Waterman), `pywfa` (wavefront alignment),
`numpy`. Python 3.10+. Two more are needed only by the extra outputs, and only
when you ask for them: `matplotlib` for `--plot`, and `mmseqs2` for the
clustering flags. Both are in `environment.yml`.

## Usage

```
Kmer2LTR [-o OUT] [-u RATE] [--cs] [-t THREADS] [--resume] [-v] input.fa[.gz]
```

| flag | meaning |
|---|---|
| `-o, --output` | output TSV (default: stdout) |
| `-u, --mutation-rate` | neutral substitution rate per site per year; fills the `k2p_time` column |
| `-t, --threads` | worker processes (default: 1) |
| `--cs` | emit a minimap2-style `cs` string instead of an extended CIGAR |
| `--resume` | skip records already in the output file and append |
| `-v, --verbose` | per-step progress |

`-u` has no default on purpose. The rate is a property of your species, not of
the software, and a wrong one silently rescales every age in the table; without
it `k2p_time` is `NA` and nothing else changes. For *Arabidopsis thaliana*,
`7e-9` is the usual choice.

Input may be plain or gzipped, wrapped or unwrapped, and may contain IUPAC
ambiguity codes, mixed case, whitespace inside sequence lines, and duplicate IDs.
Everything outside `ACGT` becomes `N`, which scores zero against everything and
is excluded from the substitution counts.

**Output is one row per input record, in input order, always.** Nothing is
silently dropped, so the TSV doubles as a complete accounting of the input and as
a filter: `awk -F'\t' '$3=="pass"'` gives you the clean intact elements.

### Extra outputs

None of these run by default. All of them are named after `-o`, so `-o` is
required, and all of them take only `status == "pass"` records — a pair the
significance gate declined to support is still in the table, but it is not an
element these outputs should be asserting things about.

| flag | writes |
|---|---|
| `--trim-flanks` | `<out>.trimmed.fa` — your elements with the called flanks cut off |
| `--ltr-cluster` | `<out>.consensus.fa` and `<out>.consensus_id<ID>_cluster.tsv` |
| `--internal-cluster` | `<out>.internal_id<ID>_cluster.tsv` |
| `--perfect-ltr-rt {5p,3p,consensus}` | `<out>.perfect_<mode>.fa` |
| `--plot` | `<out>.density.pdf` |
| `--min-seq-id` | cluster at one identity instead of the 0.70–0.98 sweep |

Chaining is safe: `Kmer2LTR` refuses to start if any output it is about to open
would land on the input file, rather than truncating the file it is about to
read. Identities are used to two decimal places, and a finer `--min-seq-id` is
refused rather than silently rounded. If a clustering does not produce every
table you asked for, `Kmer2LTR` says so, keeps the internal FASTA it would
otherwise have deleted, and **exits non-zero** so a pipeline notices.

**`--trim-flanks` gives you back your own input, corrected.** Structure-based
annotators often hand you flanking genome along with the element; the
`flank5_len`/`flank3_len` columns say how much, and this writes the elements
with that much removed. When a header carries a locus the coordinates are
shifted inwards to match. Both common placements are recognised —
`chr1:1000-2000#LTR/Gypsy` and `bedtools getfasta`'s
`TE_1#LTR/Copia::chr1:1000-2000` — and `..` is accepted for `-`. Headers it
cannot parse exactly are left alone.

Every sequence these outputs cut from your input is reproduced **character for
character**: soft-masking and ambiguity codes survive, even though the alignment
itself ran on an uppercased ACGT/N copy. Only the IUPAC consensus is
synthesised, and it is uppercase.

**`--ltr-cluster` groups elements into families.** It builds one IUPAC consensus
LTR per element and clusters those with mmseqs across an identity sweep, keeping
one table per identity — low identities lump at the lineage level, high ones
split toward recent bursts, and which you want depends on the question.
`--internal-cluster` does the same over the internal regions, which is a largely
independent opinion on the same grouping and so a useful cross-check; its FASTA
is scratch and is deleted once the tables exist. Clustering parameters were
grid-searched on *Arabidopsis* against singleton count and family mixing:
mmseqs beat cd-hit-est, and the consensus LTR beat both the full-length element
and the 5' LTR alone as input.

The consensus costs nothing. It is read straight off the alignment `Kmer2LTR`
already made to measure K2P, so a consensus and the divergence beside it are two
readings of one alignment and cannot disagree.

**`--perfect-ltr-rt`** rewinds each element's LTRs to insertion: the internal
region flanked by one identical LTR copy on both ends, so the pair is identical
again. `5p` and `3p` use that copy verbatim; `consensus` uses the IUPAC
consensus. Headers gain `~LTRlen:<n>`.

**`--plot`** draws the K2P density with a second axis in millions of years when
`-u` is given. The axis is an exact relabelling of the same numbers, not a
second estimate.

`--resume` cannot be combined with any of these: the table resumes on its
data-line count, but these files hold only the passing subset, so that count
cannot say where they stopped.

One caveat on duplicate record IDs. The TSV is positional, so duplicates there
are harmless — but these files are keyed by ID, and a cluster table naming a
duplicated ID cannot be joined back to one element. `Kmer2LTR` warns once on
stderr if that happens.

### Advanced flags

Every default is set from the benchmarks in `docs/benchmarks.md` rather than by
hand, so you should rarely need these.

| flag | meaning |
|---|---|
| `--flank-bits` | pin the evidence required to call a flank. The default is a schedule keyed on each element's own estimated divergence. |
| `--flank-sensitivity` | `strict` (default), `balanced`, `sensitive` — see below |
| `--min-bitscore` | additional floor on the reported alignment score |
| `--max-window` | cap the prefix/suffix search window |

**`--flank-sensitivity` is the one knob that encodes something about your data
rather than about the sequence.** A flank of length *k* can only ever supply
`k × α` bits of evidence that it is a flank, where α is the per-base cost of
aligning unrelated sequence under that element's own matrix. Since α falls with
divergence (≈4.2 bits/bp on a young element, ≈0.7 on an old one), a fixed
threshold is unreachable below a few bases on young elements and below a dozen on
old ones — short flanks are undetectable by construction, not by accident.

`strict` keeps the fixed threshold and suits tightly-extracted structural
predictions. The looser settings cap the demand at what the flank could supply,
which recovers short-flank detection at a real cost in false flanks:

| setting | 5 bp flanks found | false flanks on truly-unflanked elements |
|---|---|---|
| `strict` | 30% | 0.6% |
| `balanced` | 56% | 7.8% |
| `sensitive` | 79% | 10.7% |

Use the looser settings when your input is padded with genomic context. Note the
looser settings also reduce K2P bias — 37% lower at `balanced` — because the
boundary stops absorbing flank bases into the reported LTR.

## Output columns

```
1  seq_id        7  ltr3_end     13 n_sites     19 k2p                25 k2p_time
2  seq_len       8  ltr5_len     14 n_ts        20 k2p_se
3  status        9  ltr3_len     15 n_tv        21 bitscore
4  ltr5_start   10  flank5_len   16 n_gapcols   22 flank_margin_bits
5  ltr5_end     11  flank3_len   17 identity    23 cigar
6  ltr3_start   12  aln_len      18 p_dist      24 motif
```

`motif` and `k2p_time` sit *after* `cigar` rather than before it, so every
column that predates them keeps its number and no existing `cut -f` shifts.

Coordinates are **1-based inclusive**, relative to the sequence as you supplied
it, so `seq[ltr5_start-1:ltr5_end]` is the 5' LTR.

The six fields most people want are `cut -f1,4-7,19,23`.

`flank_margin_bits` says how decisively each boundary was called. Values near
zero mark elements sitting at the detection floor, which a cautious analysis may
prefer to exclude.

`motif` is the two dinucleotides the called boundary landed on, lowercased —
`tg...ca` for canonical termini. Nothing in the tool looks for it, so the rate at
which it comes back canonical is a check on the boundaries rather than a
restatement of them.

`k2p_time` is `round(k2p / (2 × µ))` in years, and `NA` without `-u`. The
divergence spans two branches, one per LTR copy, which is where the factor of
two comes from.

### `status`

| value | meaning | data columns |
|---|---|---|
| `pass` | a significant LTR pair, K2P defined | all populated |
| `weak_pair` | pair located and measured, but not significant | all populated |
| `k2p_undefined` | pair located, divergence saturated | all but `k2p`/`k2p_se` |
| `no_pair` | no terminal repeat found | `NA` |
| `too_short` | under 100 bp | `NA` |
| `all_ambiguous` | entirely `N` | `NA` |

`weak_pair` exists because a failed significance test is not a reason to throw
away a measurement. The coordinates, counts and divergence are still real; only
the claim "this is an LTR retrotransposon" is not supported. Filter on
`status == "pass"` and you get exactly the old behaviour.

## How it works

Five stages, all per-element, no training and no reference database.

1. **Discovery.** Smith-Waterman of a prefix window against a suffix window,
   score-only, with the window doubling until the hit neither touches an edge nor
   looks like noise. Exact dynamic programming rather than seed-and-extend: at
   30% divergence a shared exact 15-mer occurs at only 1.3% of positions, so
   k-mer seeding develops a sensitivity cliff precisely in the regime this tool
   exists to serve.

2. **Calibration.** Estimate this element's own divergence, transition/
   transversion ratio, base composition and indel rate from the pair itself, and
   build the log-odds substitution matrix and affine gap penalties they imply.
   Scores come out in bits, so "is there an LTR pair here?" becomes an E-value
   rather than an arbitrary identity cutoff.

3. **Boundary model selection.** Ask directly whether homology reaches each
   terminus, as an exact model comparison rather than a greedy alignment
   endpoint. This is the stage that makes the tool usable: plain Smith-Waterman
   trims the true terminus whenever the terminal bases are diverged, and would
   falsely report overextension on roughly three quarters of perfectly bounded
   inputs at 25% divergence.

4. **Outermost pair.** If a flank is still called, look for a significant repeat
   pair outside it, so a retained nested element reports the outer element rather
   than the younger nested one. Stages 2–3 are re-run on whatever it finds.

5. **Refinement.** Global wavefront alignment of the final pair, then count
   transitions and transversions over ungapped unambiguous columns and apply the
   Kimura (1980) correction with its variance estimator.

`docs/design.md` has the full rationale, including the measured evidence for
every default and the alternatives that were tried and rejected.

## Performance

Streaming throughout — one record in memory at a time, so whole-genome-scale
input is fine. On 20 cores: 256,776 records across seven repeat libraries and
genomic call sets in 6 min 52 s; the largest single set (75,180 records,
674 Mbp) in 3 min 56 s.

Output is byte-identical regardless of thread count, and identical whether or
not any extra output was requested. `--resume` counts the data lines already
written and appends, which stays correct even with duplicate IDs; a partial
final line left by a killed run is discarded first, and resuming onto a file
that does not exist yet still writes the header.

## Validation

`docs/benchmarks.md`. In short: boundaries are scored against elements whose
truth is a construction rather than an annotation, so there is no annotation
error to confound the measurement; against real elements under known
perturbation; and against negative controls of real non-LTR transposons plus a
composition-matched shuffled null. The `TG`…`CA` and target-site-duplication
signals the tool never uses are then checked as independent confirmation.

## Tests

```bash
pytest
```

310 tests on small synthetic fixtures — no real data required. The handful
that shell out to `mmseqs` or import `matplotlib` skip themselves when those are
absent.

## Citation

If this is useful in published work, please cite the repository.

## Licence

MIT. See `LICENSE`.
