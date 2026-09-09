# Kmer2LTR

Find the two LTRs in an LTR retrotransposon (LTR-RTs) and report how far they have
diverged.

Provide a fasta of putative intact elements and it classifies where each LTR starts and ends, the Kimura 2-parameter (K2P)
distance between them, and the alignment that produced it.

```bash
Kmer2LTR elements.fa.gz -o elements.tsv -t 20
```

> This is a rewrite, and the current version. The original k-mer/MAFFT
> implementation is preserved unchanged on the
> [`legacy`](https://github.com/cwb14/Kmer2LTR/tree/legacy) branch.

## What it assumes

**the 5' and 3' ends of each input sequence are homologous to
each other.** 

That matters in practice because the tool never looks for `TG…CA` termini or
target-site duplications (TSD), those signals stay available as *independent* checks on
its output. They are used that way in `docs/benchmarks.md` and nowhere in the
code. The `motif` and `tsd` columns report whichever dinucleotides and whichever
duplication the boundary landed on, never used to find it,
which is exactly what makes them worth reading.

One flag can change that, and only if you ask for it. `--tsd-anchor` lets a
target-site duplication argue against calling a flank; turning it on spends the
TSD's independence to buy a boundary correction. It is off by default, and
`docs/benchmarks.md` shows why the evidence does not currently support turning it
on.

Note: your input LTR-RTs may have been detected using pipelines that are tuned for boundary classification using `TG…CA` termini and/or TSDs. 
LTRharvest and LTR_finder both have parameters for this, and the latter tool uses `TG…CA`, which cannot be turned off. 

## What it is for

Two things, and the second comes free with the first.

1. **Dating insertions.** An LTR-RT inserts with two identical
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
git clone https://github.com/cwb14/Kmer2LTR.git
cd Kmer2LTR
mamba env create -f environment.yml
mamba activate kmer2ltr
pip install -e .
cd ..
```

## Usage

```
Kmer2LTR [-o OUT] [-u RATE] [--cs] [-t THREADS] [--genome REF] [--resume] [-v] input.fa[.gz]
```

| flag | meaning |
|---|---|
| `-o, --output` | output TSV (default: stdout) |
| `-u, --mutation-rate` | neutral substitution rate per site per year; fills the `k2p_time` column |
| `-t, --threads` | worker processes (default: 1) |
| `--genome` | reference the input was cut from; fills the orientation and TSD columns |
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

**Give it `--genome` if your input might be stored reverse-complemented.** A
reversed record's 5' terminus sits at the header's `end`, so the two trims have
to be applied to the opposite coordinates; nothing in the sequence says which
case you are in. Without a reference, forward storage is assumed.

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

`--resume` also cannot tell that you have changed a flag since the run it is
continuing. Resuming a `--genome` run onto a table written without one appends
populated columns after `NA` ones, and the file gives no sign of it — the same
hazard `-u` has always had, now four columns wider. Resume the command you ran,
not a different one.

One caveat on duplicate record IDs. The TSV is positional, so duplicates there
are harmless — but these files are keyed by ID, and a cluster table naming a
duplicated ID cannot be joined back to one element. `Kmer2LTR` warns once on
stderr if that happens.

### Genomic context

```bash
Kmer2LTR elements.fa.gz --genome ref.fa.gz -o elements.tsv -t 20
```

Three things a record cannot tell you about itself, and one reference can.

**Which strand it is stored on.** Annotation pipelines routinely
reverse-complement an element while leaving forward genomic coordinates in its
header, and no amount of staring at the sequence reveals that — a direct
terminal repeat stays one under reverse-complement, which is the same property
that lets `Kmer2LTR` ignore strand everywhere else. It matters in exactly one
place: `--trim-flanks` shifts those coordinates inwards, and on a reversed
record the 5' trim belongs to the header's `end`. Applying it the forward way
translates the interval by `flank5_len - flank3_len` while leaving its *length*
correct, so a length check does not catch it. On a three-genome *Arabidopsis*
call set 30.7% of records are stored reversed and 7.9% of `--trim-flanks`
headers were wrong by a mean of 26 bp; on a human set, 9.9% and 2.6%. Without
`--genome` forward storage is assumed, which is what every earlier version did.

**Target-site duplications.** A TSD lies *outside* the element, so without the
reference it can only be looked for when a flank was called — precisely the
records whose boundaries are least certain. With one it is available at every
boundary, and it is the sharpest signal in the table: at the boundaries
`Kmer2LTR` calls perfectly bounded it is present 22–40× more often than at a
control position ten bases away. The range is not noise — it is which detector
called the element, and `docs/benchmarks.md` is largely about why that matters.

**Whether the input was already perfectly bounded**, which is what `tsd_input`
reports and what `--tsd-anchor` can act on.

Headers are matched as `chr1:1000-2000#LTR/Gypsy` or `bedtools getfasta -name`'s
`TE_1#LTR/Copia::chr1:1000-2000`; `..` is accepted for `-`. A record whose
header carries no locus, or names a contig the reference lacks, gets `NA` in all
four columns and is otherwise untouched. Several references may be given —
`--genome a.fa b.fa` — so a multi-species call set needs no concatenation.

Each terminus is anchored **independently**, the 5' end from `start` and the 3'
from `end`. That is not fussiness: annotation sets routinely contain records
shorter than their header span, because a nested inner element was excised while
the header kept the outer interval (14.6% of the *Arabidopsis* set, 29.1% of the
human one). The middles do not correspond; both termini still do. An end that
fails to anchor reports `NA` rather than a value read from the wrong place.

The reference is read **once, streaming**, in bounded blocks rather than by
line, so a contig is never held whole even when the reference is written
unwrapped: GRCh38.p14 passes in 27 s at a 38 MB peak, and a 50 Mbp single-line
contig costs nothing measurable. No index, no `bgzip`, no extra dependency;
plain gzip is fine, and both coordinate conventions are accepted, so a
`bedtools getfasta` header works as well as a 1-based one. What memory the run
does use is the harvest, which scales with your record count and not with the
reference — about 700 bytes each. Measured cost of the extra pass:
*Arabidopsis* 26 s → 31 s, GRCh38.p14 18 s → 45 s.

### Advanced flags

Every default is set from the benchmarks in `docs/benchmarks.md` rather than by
hand, so you should rarely need these.

| flag | meaning |
|---|---|
| `--flank-bits` | pin the evidence required to call a flank. The default is a schedule keyed on each element's own estimated divergence. |
| `--flank-sensitivity` | `strict` (default), `balanced`, `sensitive` — see below |
| `--period-rule` | `best-score` (default) or `outermost` — which pair wins when a record offers more than one; see below |
| `--min-bitscore` | additional floor on the reported alignment score |
| `--max-window` | cap the prefix/suffix search window |
| `--tsd-anchor` | bits of credit a TSD at a record's own termini gets against calling a flank there. Needs `--genome`; `0` (off) by default. If you want to use, try setting 4-8 to be conservative |

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

**`--period-rule outermost` is for elements whose LTRs carry a tandem repeat.**
By default the tool reports the highest-scoring alignment it can find between the
start and the end of a record. That is the LTR pair almost always. It is not the
LTR pair when the record contains an internal duplication that is *longer* than
the LTRs. The aligner locks onto that register instead, and since the two
registers differ by a whole number of repeat units, the pair it reports spans
less of the record than the true pair does.

You can recognise it in the output without knowing the right answer, as long as
your input is tightly extracted. The called pair stops short of the ends of the
record, so `flank5_len` or `flank3_len` comes back non-zero on an element that
should have no flank at all. Across the 40 elements below, 25 of which the
default rule gets wrong, that symptom caught 23 of the 25 and raised one false
alarm among the 15 it gets right. Two weaker signs usually travel with it, an
inflated `ltr5_len` and an internal region shorter than the family's.

`outermost` reads the same two windows as a set of candidate periods and takes
the pair reaching furthest towards both termini, among those that stay
significant and still leave at least 100 bp between the two copies. A candidate
only wins if it contains the pair the default would have picked, so switching the
flag can only widen a located pair, never trade one side for the other. Where the
first alignment already runs from the first base to the last, it returns
immediately, which is the common case.

On 40 real elements drawn from seven source accessions, 20 chosen because the
boundaries were wrong and 20 length-matched ones that were right:

| set | `best-score` correct | `outermost` correct |
|---|---|---|
| known-bad (n=20) | 0 | 20 |
| controls (n=20) | 15 | 20 |

The five controls that move were mislabelled: they carry the same displacement,
just less obviously. Nothing regressed. Reproduce with
`python bench/period_fixtures.py <fixture_dir>`.

The default stays `best-score` because every calibrated constant in this tool —
the flank threshold schedule, `MAX_EVALUE`, the divergence-aware `T_BITS` — was
measured under it. Switch deliberately, per run.

## Output columns

```
1  seq_id        7  ltr3_end     13 n_sites     19 k2p                25 k2p_time
2  seq_len       8  ltr5_len     14 n_ts        20 k2p_se             26 orientation
3  status        9  ltr3_len     15 n_tv        21 bitscore           27 tsd
4  ltr5_start   10  flank5_len   16 n_gapcols   22 flank_margin_bits  28 tsd_offset
5  ltr5_end     11  flank3_len   17 identity    23 cigar              29 tsd_input
6  ltr3_start   12  aln_len      18 p_dist      24 motif
```

Everything after `cigar` was *appended* rather than inserted, so every column
that predates it keeps its number and no existing `cut -f` shifts. Columns 26–29
need `--genome` and are `NA` without it.

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

`orientation` is `+` or `-`: whether the record is stored as the reference holds
it or reverse-complemented, relative to its own header coordinates. `-` is not a
statement about the element's biology — an LTR-RT has no preferred strand — only
about how the file you were given was written.

`tsd` is the target-site duplication at the boundary `Kmer2LTR` called,
uppercase, or `.` if there is none. It is the companion to `motif` and is read
the same way: off the answer, never used to find it. `tsd_input` is the same
measurement at the record's termini *as you supplied them*, so it describes the
annotator rather than the tool. The two are the same measurement whenever no
flank was called, and they separate exactly where the two disagree — which is
the interesting population, and the one `docs/benchmarks.md` splits by detector.

`tsd_offset` is how far each boundary had to move for `tsd` to appear, as
`d5,d3`; positive means *into* the element. `0,0` is a duplication sitting
exactly on the called boundary, `0,1` one whose 3' boundary is a base too far
out. It is `NA` when `tsd` is `.`.

A duplication is accepted only if the two k-mers match exactly, contain no `N`,
and hold at least two distinct bases — a homopolymer run matches its own
reflection almost anywhere in a genome and would swamp the signal.

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

395 tests on small synthetic fixtures — no real data required, including a
miniature reference genome for `--genome`. The handful that shell out to
`mmseqs` or import `matplotlib` skip themselves when those are absent.

## Citation

If this is useful in published work, please cite the repository.

## Licence

GNU General Public License v3. See `LICENSE`.

Copyright (c) 2026 Chris Benson.
