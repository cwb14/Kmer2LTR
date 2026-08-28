# ltrk2p — LTR boundary classification and K2P divergence

**Date:** 2026-08-28
**Status:** design approved, pending spec review

## 1. Purpose

Given a multi-sequence FASTA of putative intact LTR-RTs, locate the 5' and 3' LTRs
in each record and report the Kimura 2-parameter distance between them.

The tool assumes exactly one thing: **the 5' and 3' regions of the input sequence are
homologous to each other.** It makes no assumption about element family, LTR length,
internal content, terminal motifs, or strand. (A direct terminal repeat remains a direct
terminal repeat under reverse-complement, so no strand handling is required.)

Secondary purpose, arising for free from the primary one: detect **overextension**, where
a de novo structure-based annotation has included flanking genomic sequence beyond the
element. This is reported through the coordinates themselves, not a separate mechanism.

### Non-goals

- Classifying element family or superfamily.
- Handling solo LTRs, internal-only entries, or non-LTR TEs as anything other than
  a reported non-result.
- Detecting nested elements as such. Where a nested internal element was not excised,
  the tool reports the **outer** pair.

## 2. Contract

**Input:** one FASTA file, plain or gzip, read as a stream. Records may be wrapped or
unwrapped, may contain whitespace inside sequence lines, IUPAC ambiguity codes, mixed
case, and duplicate IDs.

**Output:** one TSV, exactly **one row per input record, in input order, always.**
Line count of the output body equals record count of the input. Nothing is silently
dropped, so the output doubles as a complete accounting of the input and as a clean
intact-LTR-RT filter.

**Coordinates:** 1-based, inclusive, relative to the input sequence *as given*.

A perfectly bounded element therefore reports `ltr5_start = 1` and
`ltr3_end = seq_len`. Any other value is overextension, quantified directly by
`flank5_len = ltr5_start - 1` and `flank3_len = seq_len - ltr3_end`.

Note that input headers sometimes carry genomic coordinates whose span disagrees with the
sequence length (nested-processed outer elements, whose inner insertion was excised while
the header retained the original genomic interval). Headers are treated as opaque
identifiers and are never parsed for coordinates. All output coordinates are offsets into
the supplied sequence.

## 3. Algorithm

### Stage 0 — sanitise (streaming)

Uppercase; strip whitespace within sequence lines; map every character outside `ACGT`
(IUPAC codes, `*`, `-`, digits) to `N`. `N` scores 0 against everything (uninformative)
and is excluded from K2P site counts. Duplicate IDs pass through untouched — the
row-per-record guarantee makes them harmless.

Records are read one at a time. Whole-genome-scale inputs (the poa set is 674 Mbp across
75,180 records) never reside in memory at once.

### Stage 1 — discovery, adaptive window

Let `L = len(S)` and `wstart = L - W`. Set `W = min(floor(L/2), 1500)`. Take prefix
`P = S[0:W]` and suffix `Q = S[wstart:L]`.

Score-only SIMD Smith-Waterman of `P` vs `Q` returns the maximal score and its end
position `(qe, re)` — the two **inner** LTR ends. A second score-only pass on
`reverse(P[0:qe+1])` vs `reverse(Q[0:re+1])` returns an end position that maps back to
`(qb, rb)` — the two **outer** ends. Traceback is never run at this stage: it costs ~7x
more than score-only (1.89 ms vs 0.27 ms at 2 kb) and is not needed until the pair has
been located and its boundaries settled.

**Coordinate convention.** All internal arithmetic is 0-based half-open, matching Python
slicing and the aligner APIs. Conversion to 1-based inclusive happens once, in the output
writer. This is stated explicitly because mixing the two is the most likely source of
off-by-one errors, and the tests assert the convention by slicing the input sequence with
the reported coordinates.

Double `W` and redo — up to `floor(L/2)` — if **either** the hit touches the inner edge of
a window, **or** the hit is not significant.

The second trigger is not redundant, and omitting it is a real failure mode. When
`LTR length >= 2W`, the prefix window covers LTR-relative positions `[0, W)` while the
suffix window covers `[LTR_len - W, LTR_len)` — **disjoint ranges**, so the two windows
hold non-homologous parts of the same LTR and the best available alignment is background
noise. Noise has no reason to land near a window edge, so an edge-only trigger never
fires and the element is silently mis-called with a spurious hit. Measured with
`W0 = 1500`: a hard cliff at exactly `2*W0` — LTRs of 1400/2500/2950 bp resolved 20/20,
while 3000/4000/6000 bp resolved 2/20, 2/20 and 0/20. Maize Gypsy LTRs (1.3-2.5 kb) and
legume Ogre/Tat LTRs (4-5 kb) fall in the failing range, so this is not a corner case.

Growing on insignificance restores 20/20 at every length tested and leaves the common case
unchanged (0.33 ms for a typical 400 bp-LTR element). Non-LTR input grows to the ceiling
and is then rejected outright — an insignificant hit at `W = floor(L/2)` means there is no
terminal repeat, so discovery returns nothing rather than a noise hit a caller could not
distinguish from a real pair. Bounded at ~22 ms for a 35 kb record, and it saves the
downstream calibration and refinement work that a spurious hit would have triggered.

`floor(L/2)` is always sufficient: detection needs `W > LTR_len / 2`, and since an element
contains two LTRs, `L >= 2 * LTR_len`, so `floor(L/2) >= LTR_len > LTR_len / 2`.

**No LTR-length parameter exists anywhere in the tool**; cost tracks the actual LTR size
rather than a worst-case guess.

`W = floor(L/2)` is a hard ceiling. It makes the two LTR intervals disjoint by
construction, so `ltr5_end < ltr3_start` cannot be violated, and it is sufficient: an
element consisting of two abutting LTRs with no internal region is still fully covered.

Discovery is exact dynamic programming, not seed-and-extend. This is deliberate. At 30%
divergence a shared exact 15-mer occurs at only 1.3% of positions, so k-mer seeding
develops a sensitivity cliff precisely in the regime this tool exists to serve. Measured
cost of exact DP is 0.27 ms for a 2 kb window at ~15 GCUPS, so there is no performance
pressure to trade accuracy away.

### Stage 2 — self-calibration

Globally align the core pair and estimate identity `p̂`, transition/transversion ratio
`κ̂`, and indel rate. Build the substitution matrix implied by *this element's own*
divergence and composition:

```
s(x,y) = log2( P_xy(d̂, κ̂) / (4 · f_x · f_y) )
```

where `P` is the K2P transition-probability matrix and `f` the element's base composition.
This is the log-odds of homology against independence: the numerator is the joint
probability of the aligned pair under K2P (whose ancestral distribution is uniform, hence
the factor 4), the denominator the joint probability under independence at the element's
own composition.

**Both `f_x` and `f_y` are required.** The conditional form `P_xy / f_y` is asymmetric
whenever the composition is non-uniform — K2P's `P` assumes equal base frequencies, so
dividing by an observed skewed composition mixes two models. The consequence is not
cosmetic: it makes the alignment score depend on which sequence is the query (measured at
735 vs 755 on one AT-rich pair), which would affect nearly every element in an AT- or
GC-biased genome.

Three consequences, each load-bearing:

1. It is a full 4x4 matrix, so **transitions cost less than transversions** — which is how
   diverged LTRs actually differ, and it is the same substitution model the tool reports.
2. Scores are in **bits**, so "no LTR pair here" becomes a bit-score/E-value decision
   rather than an arbitrary identity floor.
3. Composition adjustment via `f_x · f_y` down-weights matches in AT-rich context, which is
   the principal defence against over-extension into AT-rich flanking DNA.

Re-run Stage 1 once with the calibrated matrix. This self-tuning is what delivers "zero
customization": the scoring adapts per element instead of the user tuning it.

**Integer scaling.** parasail requires integer scores, so bit-valued log-odds are scaled
by `SCALE = 4` (0.25-bit resolution) and rounded. At `p̂ = 0.7` this gives roughly
match `+6`, mismatch `-5`; `T = 5 bits` becomes `20` units. All thresholds are declared in
bits and converted at one place, so nothing downstream reasons in raw score units. Score
magnitudes on a large window can exceed int16, so the `_sat` parasail variants are used
throughout — they escalate 8 -> 16 -> 32 bit automatically on saturation rather than
silently overflowing.

**Fallback.** If the Stage 1 core is too short or too diverged to estimate `p̂` stably
(fewer than 50 ungapped sites), calibration is skipped and the generic +1/-1 matrix is
retained. This keeps a degenerate element from producing a degenerate scoring matrix.

### Stage 2b — significance

The calibrated matrix puts scores in bits, which makes "is this a real LTR pair?" a
statistical question rather than an identity cutoff. The bit score of the final alignment
is `S' = (raw_score / SCALE)`, and significance uses the Karlin-Altschul form

```
E = K * m * n * 2^(-S')
```

with `m`, `n` the two window lengths and `K` a constant for the matrix. A pair is reported
as `status = pass` only if `E` is below threshold; otherwise `status = no_pair`.

Because gapped alignment perturbs the analytic `K` and `lambda`, the threshold is **not
trusted analytically**. It is calibrated empirically against the negative controls of
section 6.6 (non-LTR TEs and shuffled sequence) to hit a target false-positive rate, and
`--min-bitscore` exposes it. The analytic form supplies the correct *shape* — the
dependence on window size, which a bare score threshold would get wrong for long
elements — while the data supplies the constant.

`bitscore` (output column 21) reports `S'`.

### Stage 3 — terminal boundary model selection

**This stage exists because of a measured failure.** Plain Smith-Waterman trims the true
terminus whenever the terminal bases are diverged, because ending one base earlier scores
better. Measured on synthetic perfectly-bounded elements (n=300 per point,
400 bp LTRs):

| p-distance | 5' terminus trimmed | mean bp lost | max | 3' terminus trimmed | mean bp lost | max |
|---|---|---|---|---|---|---|
| 0.05 | 9.3%  | 0.13 | 2  | 8.7%  | 0.16 | 5  |
| 0.15 | 30.3% | 0.90 | 13 | 27.7% | 0.67 | 10 |
| 0.25 | 48.7% | 1.69 | 20 | 43.3% | 1.33 | 20 |
| 0.35 | 65.0% | 3.76 | 44 | 63.0% | 3.40 | 47 |

At ~25% p-distance (K2P ~0.30, the stated target regime) a naive local aligner would
**falsely report overextension on roughly half of perfectly bounded inputs.** That is
disqualifying for the tool's headline feature.

The fix is to ask directly whether homology reaches each terminus, as an exact model
comparison rather than a greedy endpoint:

- **5' test:** align `reverse(S[0 : ltr5_start])` against `reverse(S[wstart : ltr3_start])`,
  begins anchored to the core, ref-end free. Score `s5`.
- **3' test:** align `S[ltr5_end+1 : W]` against `S[ltr3_end+1 : L]`, begins anchored to
  the core, query-end free. Score `s3`.

Snap the boundary to the terminus iff the extension costs less than `T` bits. Under the
penalised objective `score - T * (number of free ends)`, choosing the anchored model is
equivalent to `s_ext > -T`.

`T` is the single meaningful knob and has a clean reading: **how much evidence is required
to claim a flank exists.** A 1-3 bp terminal mismatch run costs ~1-3 bits and will not
trigger a flank call. Unrelated DNA costs ~0.8 bits/base, so a 15 bp flank costs ~12 bits
and will. `T`'s default is set by the benchmark.

Verified: parasail's `sg_de` and `sg_qe` provide exactly these semantics (begins anchored,
one specified end free) — confirmed empirically, score 50 when anchored vs 37 when forced
to skip a lead-in.

Each extension alignment is optimal given a fixed core. The approximation is that the core
alignment itself does not shift; this is validated by ablation.

### Stage 4 — outermost pair

If a flank is still called after Stage 3, search strictly outside it — `S[0:ltr5_start]`
against `S[ltr3_end:L]` — for a pair significant by the Stage 2b criterion. If one exists,
it wins, and Stages 2-3 are re-run on it.

This is what makes a retained (un-excised) nested element report the outer element's LTRs
rather than the nested element's, which are typically younger and would otherwise score
higher. The check is cheap because it runs only when a flank was called.

**This stage is provisional.** The benchmark measures how often it changes the answer on
real data; if the answer is never, it is deleted rather than kept "just in case."

### Stage 5 — refinement and reporting

Globally align the final LTR pair with **WFA** to produce the extended CIGAR, then count
transitions and transversions over ungapped, unambiguous columns.

WFA was chosen on measured evidence, not theory. Head-to-head against parasail
`nw_trace_striped_sat` on the identical task (global traceback, matched affine penalties):

| LTR length | p=0.02 | p=0.05 | p=0.15 | p=0.25 | p=0.35 |
|---|---|---|---|---|---|
| 400 bp  | 168x | 60x  | 9.6x | 3.4x | 1.9x |
| 1000 bp | 168x | 73x  | 12x  | 3.9x | 1.4x |
| 2000 bp | 307x | 132x | 14x  | 3.5x | 1.4x |
| 4000 bp | 505x | 126x | 14x  | 3.1x | 1.1x |

WFA wins at every divergence level, by 100-500x on the high-identity bulk of real data and
still 1.1-1.9x at 35% p-distance. It is exact, not heuristic: WFA and parasail return
identical optimal scores (-704) under matched penalties.

Traceback, not discovery, is the real cost centre (1.89 ms vs 0.27 ms at 2 kb), so this is
where the speed matters.

**CIGAR conversion.** pywfa emits `M` for true matches and `X` for mismatches. Because `X`
is emitted separately, `M` here unambiguously means "equal" despite the SAM spec's
contrary meaning. Conversion to extended CIGAR is `M -> =`. This is pinned down here
because silently mis-reading it would corrupt column 23.

**Open trade-off, resolved by ablation.** WFA2-lib supports only a uniform mismatch
penalty, not the ti/tv-aware calibrated matrix from Stage 2. Boundaries come from
Stages 1-4 (parasail, calibrated matrix), so only the internal alignment path is affected,
and the effect on K2P is expected to be second-order. The benchmark measures whether the
calibrated matrix materially changes K2P at high divergence. If it does, Stage 5 switches
to parasail above a measured crossover; if not, WFA is used throughout and the tool stays
simpler.

### K2P

Over ungapped columns where both bases are unambiguous:

```
P = transitions / n_sites
Q = transversions / n_sites
d = -0.5 * ln(1 - 2P - Q) - 0.25 * ln(1 - 2Q)
```

with the standard Kimura (1980) variance estimator for `k2p_se`. Gaps are excluded from
substitution counts (standard practice) and reported separately as `n_gapcols`.

When `1 - 2P - Q <= 0` or `1 - 2Q <= 0` the estimate is undefined (saturated). It is
reported as `NA` with `status = k2p_undefined`, **never as a silently clamped number.**

#### Terminal sites are not trimmed by default

A tempting safeguard is to exclude the outermost few bases of each LTR from the
substitution counts, on the grounds that boundary placement is uncertain there. This is
**not** done by default, because the bias it introduces does not point in one direction:

- At the two **outer (snapped) ends**, the boundary is fixed by the sequence terminus, not
  by the alignment score. There is no selection on those bases, so they are a fair sample.
  Excluding them costs data and removes no bias.
- At the two **inner ends**, the boundary *is* score-determined, and a local alignment
  always terminates on a match. The retained inner-terminal bases are therefore
  **match-enriched**. Trimming them would bias K2P *upward* — the opposite of the
  contamination the trim was meant to remove.
- Only when Stage 3 has **erred** and retained a few bases of true flank does trimming
  help, and there it biases K2P downward.

A blanket 5 bp margin would discard ~5% of sites (20 bp across four ends of a 400 bp LTR)
from every element, to guard against an error that Stage 3 exists to prevent, at exactly
the divergence where sites are scarcest. That trades a measured, correctable bias for an
unmeasured one.

Two things are done instead:

1. **Ablation rather than assumption** (section 6.4): K2P is measured at
   trim in {0, 3, 5, 10} bp per end for bias *and* RMSE across the divergence range and
   all flank lengths. If a trim reduces RMSE without introducing bias, it becomes the
   default. If it trades bias for variance, it does not. The evidence decides, and if the
   benefit turns out to be confined to the flank-called subset, the trim is applied only
   there.
2. **Uncertain elements are flagged, not silently degraded.** Stage 3 already computes its
   decision margin, so column `flank_margin_bits` reports how decisively each boundary
   call was made. Filtering ambiguous *elements* downstream is statistically cleaner than
   trimming sites from *every* element, and it keeps the tool's behaviour visible rather
   than hidden inside the estimator.

## 4. Output columns

```
1  seq_id        7  ltr3_end     13 n_sites     19 k2p
2  seq_len       8  ltr5_len     14 n_ts        20 k2p_se
3  status        9  ltr3_len     15 n_tv        21 bitscore
4  ltr5_start   10  flank5_len   16 n_gapcols   22 flank_margin_bits
5  ltr5_end     11  flank3_len   17 identity    23 cigar
6  ltr3_start   12  aln_len      18 p_dist
```

The six originally requested fields are `cut -f1,4-7,19,23`. `--cs` swaps column 23 to a
minimap2-style short `cs` string.

`flank_margin_bits` reports the smaller of the two Stage 3 decision margins,
`min(|s5 + T|, |s3 + T|)` in bits: how decisively the boundary call was made. Large values
mean the call was unambiguous; values near zero mark elements whose boundaries sit at the
detection floor and which a cautious downstream analysis may wish to exclude.

`status` values: `pass`, `no_pair`, `too_short`, `all_ambiguous`, `k2p_undefined`.
Non-`pass` rows carry `NA` in the coordinate and distance columns.

## 5. Package

```
ltrk2p/
├── pyproject.toml, environment.yml, README.md, LICENSE
├── src/ltrk2p/
│   ├── cli.py      argparse entry point
│   ├── fasta.py    streaming gzip-aware reader, sanitisation
│   ├── scoring.py  calibrated log-odds matrix, bit scores
│   ├── align.py    Stages 1-4
│   ├── k2p.py      distance and variance
│   ├── cigar.py    extended CIGAR and cs emission
│   └── runner.py   parallel driver, ordered writer, resume
├── tests/          pytest, tiny synthetic fixtures
└── bench/          build_truth · simulate · run_bench · figures
```

**CLI:** `ltrk2p [-o OUT] [--cs] [-t 20] [--resume] [-v] input.fa[.gz]`

Advanced knobs (`--flank-bits`, `--min-bitscore`, `--max-window`) exist and are
documented, but every default is set by the benchmark rather than by hand.

Verbosity follows the project default: minimal by default (start, milestones, done,
errors); `-v` adds per-step progress and sanity checks.

Failure policy: fail fast with a clear message on bad input, missing dependencies, or a
corrupt file; warn and skip for a single malformed record.

**Parallelism:** `ProcessPoolExecutor` over records with chunking; output written in input
order. Records are independent, so this scales linearly.

**Resume:** `--resume` counts data lines already present in the output and skips that many
input records, then appends. This is correct even with duplicate IDs, which an ID-keyed
resume would not be.

## 6. Validation

### 6.1 Ground truth by construction

Repbase, Dfam, MTEC and riceTElib split many families into `X-LTR` and `X-I` entries
(Repbase alone: ~33k `-LTR`, ~30k `-I`). Pairing these by family and concatenating
`LTR + I + LTR` yields full-length elements whose boundaries are known **exactly** — the
truth is a construction, not an annotation, so there is no annotation error to confound
the measurement.

Messy entries (unpairable, duplicate, too short, mostly `N`) are dropped rather than
repaired. There is ample material.

### 6.2 Simulation

From each truth element, evolve the two LTR copies **independently from the ancestral
consensus**, each for t/2, so the pairwise distance equals the target d. Parameters:
d in {0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6};
kappa in {1, 2, 4}; indels at a realistic rate with geometric lengths.

Flanks of length {0, 1, 2, 5, 10, 20, 50, 100, 200, 500} are drawn from
**dinucleotide-shuffled real TE sequence** — composition-matched, so flank detection is
tested against a hard null rather than against uniform random DNA.

Truth is recorded two ways, and the distinction matters:

- **nominal d** — the simulated divergence.
- **realized K2P** — computed from the true alignment.

Comparing against realized K2P isolates *tool* error. Comparing against nominal d measures
the K2P estimator itself, which is a separate and already-known quantity. Reporting both
keeps them from being confounded, which is the usual way this kind of benchmark misleads.

### 6.3 Metrics

- Per-coordinate boundary error distributions vs. divergence.
- Flank-detection ROC vs. flank length.
- K2P bias and RMSE vs. d, out to d = 0.6.
- Failure rate; runtime and scaling.

### 6.4 Ablations

Every design element must earn its place:

| Ablation | Expectation |
|---|---|
| BLASTN-style 1/-3 scoring baseline | fails badly past ~20% divergence (break-even ~75% identity) |
| fixed +1/-1 vs. calibrated log-odds | quantifies what Stage 2 buys |
| Stage 3 on/off | headline: the ~46% false-flank rate at p=0.25 should collapse |
| Stage 4 on/off | delete Stage 4 if it never fires |
| WFA uniform penalty vs. calibrated matrix in Stage 5 | decides the Stage 5 aligner |
| Stage 1 iteration count | confirms one recalibration pass suffices |
| terminal trim of 0/3/5/10 bp per LTR end | decides whether trimming is adopted; measured for bias *and* RMSE, overall and within the flank-called subset |

### 6.5 Real data, three independent checks

1. All seven datasets end-to-end (arab 10,307; poa 75,180; human 16,336; plus the four
   libraries): status breakdown, K2P distributions, runtime.
2. **The TG...CA check.** Because the tool uses no structural priors, the fraction of
   called boundaries landing exactly on `TG`...`CA` is an *unbiased external accuracy
   proxy* on real data where no truth exists. The tool never sees this signal, so it
   cannot game it. Declining the prior is precisely what makes this measurement
   informative.
3. **TSD detection** at called flanks — a second independent signal, same logic.

### 6.6 Negative controls

Dfam and Repbase non-LTR entries (DNA transposons, LINEs, SINEs, Helitrons) plus shuffled
sequence, measuring the false-positive rate for `status = pass`. This sets
`--min-bitscore`.

### 6.7 Cross-check

K2P arithmetic verified against an independent implementation on identical alignments, so
that any tool-vs-truth gap cannot be attributed to the distance calculation.

### 6.8 Unit tests

pytest with tiny synthetic fixtures (not real data):

- K2P against hand-computed values and an independent implementation.
- Extended CIGAR and `cs` round-trip: both sequences reconstructable from the string.
- Identical LTRs -> exact boundaries, K2P = 0.
- Known flanks -> detected at the correct lengths.
- Solo LTR / random sequence -> `status = no_pair`.
- `N`, IUPAC, whitespace, mixed case, blank lines, gzip, empty file, duplicate IDs.
- Window growth: element whose LTR exceeds the initial window.
- Nested element -> outer pair reported.
- Coordinate convention: 1-based inclusive, verified by slicing the input sequence.

### 6.9 Compute

Benchmark jobs: `bio250178` allocation, `shared` partition, 20 threads. Each run produces
a `memo_*.md` recording date, purpose, environment and versions, input provenance, exact
commands, expected outputs, and notes.

## 7. Known limits — to be measured, not assumed

- **Inner-boundary trimming.** Stage 3 anchors the two *outer* boundaries. The inner
  boundaries have the same trimming tendency with no terminus to snap to, and the trimmed
  bases are enriched for mismatches, so K2P may be biased slightly low. This will be
  quantified. It will be corrected only if material, since the trimmed boundary is the
  maximum-likelihood answer under the model.
- **Short flanks (<5 bp) are undetectable in principle** — indistinguishable from terminal
  divergence. This will be reported as a measured detection floor rather than papered over.
- **AT-rich and low-complexity flanks** may invite over-extension. Composition-adjusted
  scoring is the mitigation; the benchmark determines whether it suffices. A
  low-complexity diagnostic column is added only if the data demands one.
- **Composition adjustment can invert `match > transition > transversion` at high
  divergence under strong base-composition skew** (d >~ 0.5 with skew as mild as 65/35;
  and specific cross-pair comparisons such as a common-base transition scoring below a
  rare-to-rare transversion can invert at lower d). This is not a defect: matching a very
  common base genuinely carries less evidence than a substitution involving a rare one, and
  it is the same rare-residue log-odds inflation seen in BLOSUM-style matrices. It is a
  property of *any* composition-adjusted log-odds matrix, present identically before and
  after the symmetry fix, on exactly the same set of (composition, d, kappa) combinations.
  The `calibrated vs fixed +1/-1` ablation in section 6.4 is what determines empirically
  whether it costs accuracy at the divergences this tool targets; if it does, the remedy is
  to drop composition adjustment, not to hand-patch the ordering.
- **Integer scaling at `SCALE = 4` can round a transition and a transversion score to a
  tie** at the extreme edge (measured 20 of 2880 grid points, all at d >= 0.4 with strong
  skew). It never *inverts* a correctly ordered pair (0 of 2880). Pre-existing and
  unaffected by the symmetry fix; raising SCALE is the remedy if the benchmark shows it
  matters.
- **Saturation.** Beyond d ~ 0.6-0.7 K2P becomes unstable and eventually undefined. The
  tool reports `NA` rather than a misleading number, and the benchmark characterises where
  the variance becomes unacceptable.

## 8. Environment

conda env `ltrrt4` at `/anvil/projects/x-bio250178/conda/envs/ltrrt4`:
python 3.11, parasail-python, pywfa (built against the env toolchain; the system OpenMPI
compiler wrapper must not be used), numpy, scipy, pandas, matplotlib, pytest, biopython.

Shipped as `environment.yml`.
