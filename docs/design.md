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

**The composition is the whole element's, not the LTR core's.** Estimating it from
the core alone was measured (2026-08-31) and is worse on every panel of both
benchmark grids: on the homology grid's substitution panel it moves boundary MAE
3.29 → 3.79, 25 bp flank detection 79.2% → 74.1% and the false-flank rate 0.20% →
0.41%. Consequence 3 above is why: the denominator is the *null* a match is being
judged against, and the null for "does homology continue past this boundary?" is
the surrounding sequence, not the repeat itself. The whole element — internal
region included — is the better proxy for that.

Re-run Stage 1 once with the calibrated matrix. This self-tuning is what delivers "zero
customization": the scoring adapts per element instead of the user tuning it.

**Integer scaling.** parasail requires integer scores, so bit-valued log-odds are scaled
by `SCALE = 4` (0.25-bit resolution) and rounded. At `p̂ = 0.7` this gives roughly
match `+6`, mismatch `-5`; `T = 5 bits` becomes `20` units. All thresholds are declared in
bits and converted at one place, so nothing downstream reasons in raw score units. Score
magnitudes on a large window can exceed int16, so the `_sat` parasail variants are used
throughout — they escalate 8 -> 16 -> 32 bit automatically on saturation rather than
silently overflowing.

**Gap penalties are calibrated too (2026-08-31).** Stage 2 has always been specified
as estimating an indel rate, and the shipped code never used one: gaps were fixed at
6 bits to open and 2 bits to extend everywhere. Under a geometric indel-length model
the log-odds cost of a gap of length `k` is `-log2(mu) - (k-1)*log2(P_continue)`, which
is exactly parasail's affine form — so `open = -log2(mu)` with `mu` the per-column
probability that a gap opens, and `extend = -log2(1 - 1/mean_length)`, both read off
the same core alignment the matrix comes from. Clamped to 4-16 bits and 0.25-3 bits
respectively: a core with no gaps at all implies an infinite opening cost, and a core
with one long gap an infinite mean.

The measured effect is not second-order. Against the fixed 6/2 scheme, on the
homology grid's substitution panel, per-element gaps move flank detection at 5 bp
from 7.9% to 23.5% and at 25 bp from 79.2% to 91.1%, boundary MAE from 3.29 to
2.60, and pair loss from 2.16% to 1.95%; on the indel panel they cut K2P RMSE by
35% and K2P bias by 52%. A fixed 10/1 scheme captures part of the same gain but
is beaten by the per-element estimate on every metric of the real-composition
panel.

**Fallback.** If the Stage 1 core is too short or too diverged to estimate `p̂` stably
(fewer than 50 ungapped sites), calibration is skipped and the generic model is
retained **whole** — matrix and gap penalties both. An indel rate estimated from under
50 ungapped columns is noise, and adopting it while rejecting the matrix derived from
the same counts would be incoherent.

### Stage 2b — significance

The calibrated matrix puts scores in bits, which makes "is this a real LTR pair?" a
statistical question rather than an identity cutoff. The bit score of the final alignment
is `S' = (raw_score / SCALE)`, and significance uses the Karlin-Altschul form

```
E = K * m * n * 2^(-S')
```

with `m`, `n` the two window lengths and `K` a constant for the matrix. A pair is reported
as `status = pass` only if `E` is below threshold; otherwise `status = no_pair`.

**Significance is always scored with the generic model, never the calibrated one
— and the generic model now means the generic matrix _and_ a fixed pair of gap
penalties (`SIG_GAPS`, 6 bits open / 2 bits extend).** Pinning the gap penalties
alongside the matrix is what allows the *alignment* gap penalties to change
(section 3, Stage 2) without moving `MAX_EVALUE` underneath them: the final
significance score is `nw(q, r, SIG_GAPS, GENERIC_MATRIX)` regardless of how the
element itself was aligned. The two
are on different score scales — at low divergence the calibrated matrix scores a match at
`+8` in SCALE units where the generic matrix scores `+4` — so a fixed E-value threshold
calibrated against one is invalid against the other. Applying the generic-calibrated
threshold to calibrated scores made a 13 bp chance match report roughly double the bits and
an E-value about 2^13 too small, which leaked spurious "outer pairs" out of unrelated
flanking DNA in Stage 4 (measured 14 of 2500 elements, every instance at low divergence;
zero after gating on the generic matrix). The division of labour is therefore explicit: **the generic matrix determines significance
everywhere.** For the element's own LTR pair, the calibrated matrix still determines the
alignment and the boundaries — it is calibrated to exactly that pair. Stage 4 is the one
exception and uses the generic matrix for its boundary recovery too, because the outer pair
it is searching for is a *different, uncharacterised* pair whose divergence the calibrated
matrix does not describe; scoring it with a matrix tuned to the inner pair both leaked false
positives and degraded the true positives it did find (outer-pair boundaries correct
118/120 with the calibrated matrix, 120/120 with the generic one). Sensitivity of the
generic-scored outer search was verified out to 35% divergence (40/40 detected and correct
at 20-35%, degrading only at an unrealistic 40%).

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

Those are **per-terminus** rates. The user-facing quantity is **per-element** — an element
is wrong if *either* terminus is — and because the two termini fail near-independently that
rate is approximately `1 - (1-a)(1-b)`, which is substantially higher. Measured directly
(n=150 per level):

| p-distance | 5' trimmed | 3' trimmed | **per-element wrong** | 1-(1-a)(1-b) |
|---|---|---|---|---|
| 0.05 | 7.3%  | 10.7% | **18.0%** | 17.2% |
| 0.15 | 26.7% | 27.3% | **48.7%** | 46.7% |
| 0.25 | 50.0% | 46.7% | **75.3%** | 73.3% |
| 0.35 | 59.3% | 64.7% | **88.0%** | 85.6% |

So at ~25% p-distance (K2P ~0.30, the stated target regime) a naive local aligner would
**falsely report overextension on roughly three quarters of perfectly bounded inputs.** That
is disqualifying for the tool's headline feature.

**Measured outcome of Stage 3** (80 perfectly-bounded elements per level, per-element rate):

| p-distance | raw SW | Stage 3 | improvement |
|---|---|---|---|
| 0.05 | 14/80 | 3/80 | 4.7x |
| 0.15 | 37/80 | 2/80 | 18.5x |
| 0.25 | 63/80 (79%) | 1/80 (1.2%) | **63x** |
| 0.35 | 69/80 | 2/80 | 34.5x |
| 0.45 | 75/80 | 3/80 | 25x |

Detection of genuine flanks is retained: at p=0.25, 20 bp flanks are found 95% of the time
and 50 bp or larger essentially always, with called length accurate to under a base (mean
49.8 for a true 50, 199.9 for a true 200). The false-positive rate on zero-flank elements is
2.5%. `T_BITS = 5.0` sits in a plateau, not on a cliff: sweeping it gives 57.5% false-flank
at t=0, 3.8% at t=5, 0% at t=10 but with true-flank detection falling to 83.8%, and at
t=25 the degenerate 0%/0% regime where no flank is ever called.

**Task 16 update — `T_BITS` raised to 10.0, on real data, not the synthetic sweep above.**
The synthetic calibration above under-stated the false-flank problem badly: Task 14 swept
`t_bits in {2,5,8,10,15,20,30}` over 260,876 records built by perturbing 1,694 REAL gold
elements (arabidopsis + human + library-consensus truth.fa) with known added flank length —
`bench/out/gold_perturbed.fa`/`gold_truth.tsv`, scored per-cell in
`bench/out/cells_tbits_*.json`, full sweep table in `bench/out/memo_bench.md` and
`bench/out/memo_bench.md`. On that real-sequence grid, the
shipped default of 5.0 gave a false-flank rate of **26.8% at d=0.3** (42.6% at d=0.4) — a
~20x gap from what synthetic-only calibration suggested. At `t_bits=10`:

| d | false-flank @ t=5 | false-flank @ t=10 |
|---|---|---|
| 0.05 | 4.13% | 0.65% |
| 0.10 | 7.62% | 2.04% |
| 0.20 | 16.12% | 4.52% |
| 0.30 | 26.76% | 8.47% |
| 0.40 | 42.61% | 17.33% |

a 2.5-3.6x reduction across the range, while large-flank detection — the practically
important failure mode, catching genuine overextension — is nearly unchanged: det@50 falls
93.0% -> 91.3% (1.7 points) and det@100 barely moves, 91.7% -> 91.6%. The cost concentrates
entirely on 10-20 bp flanks (det@10 75.5% -> 39.4%, det@20 88.8% -> 70.5%), which section 7
already documents as sitting near the theoretical detection floor. **`T_BITS = 10.0` is now
the shipped default** (`src/ltrk2p/align.py`); `--flank-bits` follows it automatically.

A **divergence-aware `T_BITS`** (varying the threshold per-element with the tool's own
`calibrate()`-estimated `d_hat`) was also measured and is a genuine, non-cherry-picked
Pareto improvement over any single fixed threshold — e.g. a modest per-bin schedule reaches
4.3% pooled false-flank rate at the *same* det@10 as the flat `t=10` recommendation (40.2%
vs 40.3%), and beats the nearest fixed threshold matching its false-flank rate (`t=15`,
3.7%) by +26 points of det@10 and +11 of det@20. **Not implemented.** This task recommends
the value, not a schedule; picking a per-bin schedule needs an explicit detection floor,
not just a false-flank target, or it overcorrects exactly where catching overextension
matters most (see `task-14-report.md`, "Does a divergence-aware threshold help?"). Recorded
here as future work, deliberately deferred.

**Stage 3 takes its extension regions from the internal region, not the discovery
window (2026-08-31).** The partner region for the 5' test used to be
`S[wstart:l3b]` — the part of the *suffix window* lying before the 3' LTR. That
had two degeneracies: it is empty whenever the hit begins exactly at the window's
inner edge, so the snap could never fire; and it overlaps the 5' LTR itself
whenever the window starts before `ltr5_end`. Bounding it by the internal region
instead removes both, and — decisively — removes `w` from this stage, which is
what makes Stage 3 re-runnable after Stage 4, where no window exists. Measured
old-code-vs-new-code at identical settings on a 20,068-record stratified sample:
99.20% of records identical, false-flank rate 7.36% → 6.30%, every detection rate
within half a point, flank-length MAE better at three of four flank sizes.

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

#### Two alternatives to the binary snap, measured and rejected (2026-08-31)

**A scored/graded extension**, replacing the binary snap with "extend to the
furthest endpoint whose running score is still above the noise floor". The
literal version of this idea — take the extension's *optimal* endpoint — is
provably a no-op: discovery's core is the Smith-Waterman optimum over the same
windows, so a positive-scoring outward extension would contradict that
optimality. Confirmed empirically at `H.max() > 0` in 0 of 1,194 records. The
non-trivial version, thresholding at `-T` instead, is measurably worse: on the
homology grid it moves `bnd_mae` 3.29 → 8.81, `flank_mae` 4.70 → 12.63, K2P bias
+0.0006 → +0.0168 and pair loss 2.2% → 3.9%, for a false-flank rate that does not
move at all. The mechanism is that a `T`-bit budget spent on the snap decision
and *again* on the endpoint lets the boundary creep `T / (bits per base)` bases
into non-homologous sequence — on a control element with 8,000 bp of random flank
it calls 7,994.

**Joint inner boundaries.** Section 7 notes that the two inner boundaries are
by-products of the opposite terminus's snap. Re-deriving both from a single
alignment anchored at the settled outer ends (`sg_qe_db`, plus a reversed `sg_de`
for the ref begin) changes 418 of 62,644 located records — so it is a real
alternative, not a no-op — but it is right more often than wrong (276 vs 140) by
*less* than it is wrong by: mean summed inner error on the records it touches
242 → 557, pooled `bnd_mae` 3.287 → 3.286. Rejected.

The second result is the more informative one: the inner-boundary error that
motivated the idea is **not** a fixable artifact of the inner ends being
unexamined. Deriving them jointly and optimally leaves the pooled error unchanged
to four significant figures, which means the residual is estimation error under
the model, not a greedy-trimming bug.

### Stage 4 — outermost pair

If a flank is still called after Stage 3, search strictly outside it — `S[0:ltr5_start]`
against `S[ltr3_end:L]` — for a pair significant by the Stage 2b criterion. If one exists,
it wins, and Stages 2-3 are re-run on it.

**The re-run was specified here from the beginning and was missing from the code
until 2026-08-31.** `outermost` returned the raw Smith-Waterman ends of the outer
pair and `classify` used them directly, so the recovered outer element was the
one pair in the tool whose termini were never tested — precisely the pair most
likely to need it, since it is older and therefore more diverged than the nested
pair that displaced it, and raw SW trims a terminus exactly when the terminal
bases are diverged. The re-run recalibrates on the outer pair (its divergence is
not the inner pair's, so neither the inner pair's matrix nor its gap penalties
describe it), re-runs the outer search under that calibrated matrix, and then
applies Stage 3. The acceptance gate is untouched and still generic-scored, so
no new pair can be admitted by the change.

Making this possible required removing the discovery window from Stage 3 — see
the Stage 3 section — because after Stage 4 there is no window to bound the
extension regions with.

This is what makes a retained (un-excised) nested element report the outer element's LTRs
rather than the nested element's, which are typically younger and would otherwise score
higher. The check is cheap because it runs only when a flank was called.

**Measured on real data, and kept.** Task 14 compared Stage 4 on/off on the raw, unperturbed
real datasets (deliberately not the gold-perturbed grid, whose gold-selection criteria
require clean non-overextended raw boundaries and so structurally exclude the nested-element
scenario this stage exists to recover):

| dataset | n | n_changed | fraction |
|---|---|---|---|
| arabidopsis | 10,307 | 110 | 1.07% |
| human | 16,336 | 299 | 1.83% |
| truth.fa (library consensus) | 31,315 | 4 | 0.01% |
| **total** | **57,958** | **413** | **0.71%** |

Non-zero, and non-trivially so on genuine raw genomic calls (up to 1.83% on human) — the
deletion criterion above is not met, so Stage 4 stays. Every changed record recovers the
outer pair exactly as designed (verified by inspection, e.g. arabidopsis
`LR999451.1:10387972-10390268#LTR/unknown/unknown`: `[6-166]`/`[2144-2297]` with Stage 4 vs.
a much shorter, more central `[850-1135]`/`[1150-1438]` without it). Full detail:
`bench/out/stage4_diff_*.json`, `bench/out/memo_bench.md`.

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

**Resolved: keep WFA everywhere, no crossover to parasail.** WFA2-lib supports only a
uniform mismatch penalty, not the ti/tv-aware calibrated matrix from Stage 2. Boundaries
come from Stages 1-4 (parasail, calibrated matrix), so only the internal alignment path is
affected. Task 14's `wfa_vs_matrix` ablation (`refine="matrix"` vs the default `"wfa"`,
identical boundaries by construction, measured on the 260,876-record real gold-perturbed
grid) found the opposite of "second order": switching to the calibrated matrix roughly
**doubles the K2P bias magnitude** in both the correctly-bounded population (-0.0087 ->
-0.0156) and the flank-called population (-0.0066 -> -0.0134), with no reliable RMSE gain
(marginally better in one population, worse in the other). No crossover divergence exists
where the calibrated matrix wins outright, so WFA stays the aligner at every divergence
level. Full table: `bench/out/memo_bench.md`, `task-14-report.md`.

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

   **Resolved: `trim=0` stays the default, in both populations.** Measured on the same
   260,876-record real gold-perturbed grid (`d_nominal` as ground truth, matching
   `gold_robustness.score_gold`'s convention), comparing trim in {0,3,5,10}:

   | population (n) | trim | bias | RMSE |
   |---|---|---|---|
   | correctly-bounded, no flank called (49,259) | 0 | -0.0087 | 0.0761 |
   | | 10 | -0.0102 | 0.0761 |
   | flank-called subset (195,404) | 0 | -0.0066 | 0.0508 |
   | | 10 | -0.0045 | 0.0508 |

   In the correctly-bounded population RMSE is flat to 3 significant figures at every trim
   depth and bias does not improve — exactly the "costs data, removes no bias" prediction
   above. In the flank-called subset, bias genuinely moves toward zero as trim grows (a real
   ~32% relative reduction, -0.0066 -> -0.0045) — the predicted correction for retained
   flank contamination — but RMSE does not move (0.0508 -> 0.0508): the variance cost of
   losing sites at these trim depths cancels the bias gain. Adopting a trim requires
   reducing RMSE, not just bias; neither population clears that bar, so no trim is applied
   anywhere. Full per-trim-depth table: `bench/out/memo_bench.md`, `task-14-report.md`.
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

`status` values: `pass`, `weak_pair`, `no_pair`, `too_short`, `all_ambiguous`,
`k2p_undefined`.

`no_pair`, `too_short` and `all_ambiguous` rows carry `NA` in every data column — no LTR
pair was located, so there is nothing to report.

**`weak_pair` is the significance gate declining to destroy a measurement.** A pair was
located, its boundaries settled and its divergence measured, and only then did the
E-value (or an explicit `--min-bitscore`) fall short. Nulling all twenty data columns at
that point discards work that is correct and useful — coordinates, substitution counts,
identity, K2P — to record a single bit of information. `weak_pair` records that bit in
`status` and reports everything else, exactly as `k2p_undefined` already does for
saturation.

`status == "pass"` keeps its meaning unchanged, so the TSV remains a clean
intact-LTR-RT filter and no downstream `pass` filter shifts. What changes is that
`grep -v pass` is no longer a synonym for "nothing was found here". The cost is zero
runtime: the gate already sat *after* Stages 3-5, so the work was being done and then
thrown away.

Verified on the 70,000-record homology grid: of 62,648 records the old gate accepted,
**zero** have a single differing field under the new one, while 1,386 records it had
nulled are now reported as `weak_pair`. The change adds rows; it moves no boundary.

Where a pair is both insignificant and saturated, `weak_pair` wins — the significance
failure is the more fundamental statement, and `k2p`/`k2p_se` are `NA` either way.

**`k2p_undefined` is deliberately different**: it means the pair *was* located and its
boundaries are valid, but the divergence is saturated so the K2P correction has no defined
value. Only `k2p` and `k2p_se` are `NA`; the coordinates, lengths, flank calls, substitution
counts, identity, p-distance, bit score and CIGAR are all reported, because they are
correct and useful. Nulling them would discard good measurements to satisfy a blanket rule.

In practice this status is close to unreachable through `classify`: the significance gate
structurally requires roughly >50% identity to accept a pair, while K2P stays defined below
p-distance 0.5, so the two conditions are nearly mutually exclusive. Systematic search over
~2500 trials (pairwise distance to 3.0, `max_evalue` relaxed to 1e6) got no closer than
p-distance 0.5017 with both K2P denominators still well clear of zero. The status is kept
because saturation is a real possibility on adversarial input and reporting `NA` is the
honest response, not because it is expected to fire.

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

**Parallelism:** `ProcessPoolExecutor` over records with **bounded submission** — at most
`threads * 4` futures in flight, popped in submission order. `Executor.map` cannot be used:
it drains its input generator completely before dispatching, which materialises the entire
FASTA in the parent and breaks the streaming guarantee on exactly the multi-threaded path
that exists for whole-genome input. Measured with the alignment stubbed out, driver peak RSS
was flat at ~39 MB across 48/192/479 MB inputs single-threaded, but grew to 52/69/96 MB at
8 threads under `map`. Output is byte-identical regardless of thread count.

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

### 6.2b Homology-only grid (added 2026-08-31)

Sections 6.1-6.2 build elements from library consensus; section 6.5's gold-subset
benchmark perturbs *real* elements but selects them with `select_gold`, which filters on
the tool's own output **and** requires canonical `TG`..`CA` termini. Both filters are
deliberate there and both are disqualifying for one specific question: *how well are
boundaries recovered from homology alone?* The first makes the answer circular; the
second smuggles in a structural prior the tool itself declines to use.

`bench/homology_grid.py` answers that question with neither. Its elements come straight
from `truth.fa`'s `X-LTR` + `X-I` + `X-LTR` constructions, which are **perfect by
construction** — the two LTR copies are literally the same string, so divergence is
exactly zero and boundaries exactly known — and are selected on length alone: no motif,
no TSD, no call by ltrk2p.

Three perturbation axes on those perfect elements:

- **substitutions**, parameterised by the *target pairwise p-distance*
  {0, 5, 10, 15, 20, 25, 30, 35}% rather than by a branch length, so "add 25% mutations"
  means 25% of sites differ between the copies. `d_for_p` inverts
  `1 − p_same(d, κ) = p_target` by bisection and records the implied K2P distance.
- **flanks** of {0, 5, 25, 45, 65} bp, drawn from a composition-matched shuffle of the
  element's own sequence, so no detector can succeed by noticing a compositional step.
- **indels**, swept independently of substitutions at
  {0.001, 0.002, 0.005, 0.01, 0.02, 0.05} events per site per branch, geometric lengths.

**The true alignment is tracked through evolution.** `evolve_tracked` returns a per-base
ancestral trace alongside each descendant; `true_alignment` merges two traces into the
exact pairwise alignment. `realized_p` and `realized_k2p` are therefore measurements, not
approximations, even under indels. This is the one thing the older simulator cannot do:
`bench/simulate.py` computes `realized_k2p` by truncate-and-compare on unaligned strings,
which Task 14 found disagrees materially with a proper alignment even at `indel_frac=0`
and is simply wrong once indels are switched on. `simulate.py` is unchanged and still
carries that flaw; nothing measured on this grid depends on it.

Two panels are reported separately: `lib` (library constructions — the non-circular one,
and the deciding evidence) and `real` (motif-free Arabidopsis elements — real composition
and real internal structure, but still selected from a prediction TSV).

### 6.3 Metrics

- Per-coordinate boundary error distributions vs. divergence.
- Flank-detection ROC vs. flank length.
- K2P bias and RMSE vs. d, out to d = 0.6.
- Failure rate; runtime and scaling.

### 6.4 Ablations

Every design element must earn its place:

| Ablation | Expectation | Result (Task 14/16, measured on the 260,876-record real gold-perturbed grid) |
|---|---|---|
| BLASTN-style 1/-3 scoring baseline | fails badly past ~20% divergence (break-even ~75% identity) | confirmed: pooled false-flank 46.4% vs calibrated's 20.2%; crosses 50% between d=0.2 (46.1%) and d=0.3 (74.5%), matching the predicted ~75% identity break-even |
| fixed +1/-1 vs. calibrated log-odds | quantifies what Stage 2 buys | fixed_1_1 sits between calibrated and blastn_1_3 throughout (pooled false-flank 28.6%); calibration roughly halves the false-flank rate fixed_1_1 would otherwise show at moderate-high d |
| Stage 3 on/off | headline: the ~46% false-flank rate at p=0.25 should collapse | confirmed: no_stage3 pooled false-flank 56.2% vs calibrated (Stage 3 on) 20.2%; on synthetic perfectly-bounded elements at p=0.25, raw SW is wrong 63/80 (79%) vs Stage 3's 1/80 (1.2%), a 63x improvement |
| Stage 4 on/off | delete Stage 4 if it never fires | **kept** — fires on 0.71% of raw real predictions (413/57,958: arabidopsis 1.07%, human 1.83%, truth.fa 0.01%), always recovering the outer pair; see Stage 4 section above |
| WFA uniform penalty vs. calibrated matrix in Stage 5 | decides the Stage 5 aligner | **WFA kept everywhere** — the calibrated matrix roughly doubles K2P bias magnitude with no reliable RMSE gain; see Stage 5 section above |
| Stage 1 iteration count | confirms one recalibration pass suffices | not measured in Task 14/16 (out of scope); `classify()` still performs exactly one recalibration pass (`discover` -> `calibrate` -> `discover` again), unchanged |
| terminal trim of 0/3/5/10 bp per LTR end | decides whether trimming is adopted; measured for bias *and* RMSE, overall and within the flank-called subset | **not adopted, `trim=0` stays default** — RMSE flat at every trim depth in both populations; trim_10 cuts flank-called-subset bias ~32% but RMSE is unchanged (variance from fewer sites cancels it); see "Terminal sites are not trimmed by default" above |

### 6.5 Real data, three independent checks

1. All seven datasets end-to-end (arab 10,307; poa 75,180; human 16,336; plus the four
   libraries): status breakdown, K2P distributions, runtime.
2. **The TG...CA check.** Because the tool uses no structural priors, the fraction of
   called boundaries landing exactly on `TG`...`CA` is an *unbiased external accuracy
   proxy* on real data where no truth exists. The tool never sees this signal, so it
   cannot game it. Declining the prior is precisely what makes this measurement
   informative.
3. **TSD detection** at called flanks — a second independent signal, same logic.

**Task 16 ran all three** against the final defaults (T_BITS=10.0, MAX_EVALUE=1e-10;
`bench/out/real_*.tsv`, job `20201250`, 256,776 records across the seven datasets in
00:05:06). Headline: on every genomic/putative-intact dataset (arabidopsis, poa, human,
plus MTEC's maize library), boundaries the tool calls perfectly-bounded (`flank_len==0`)
land on `TG`...`CA` 4-10x more often than boundaries it calls flank-corrected, and at or
above the raw-input-ends baseline — an unbiased, tool-blind confirmation the boundary logic
is doing its job. TSD presence at flank-called boundaries is 3-7x enriched over a shifted-
10bp local control at k=5 on every dataset with enough records to be meaningful. Full
tables, per-dataset breakdown, and the before/after comparison against the pre-Task-16
defaults: `bench/out/memo_real.md`.

### 6.6 Negative controls

Dfam and Repbase non-LTR entries (DNA transposons, LINEs, SINEs, Helitrons) plus shuffled
sequence, measuring the false-positive rate for `status = pass`. This sets `MAX_EVALUE`
(and, through it, `--min-bitscore`'s effective floor).

**Task 16 measured this against `bench/out/negatives.fa`** (46,823 real non-LTR TEs from
two independent libraries -- Dfam-RepeatMasker and Repbase, plus MTEC and riceTElib --
with LTR-class and ambiguous tyrosine-recombinase classes excluded) and a mononucleotide-
shuffled null drawn from the library-consensus truth set (`bench/simulate.shuffle_dinuc`,
the same "shuffle" null section 6.2 already uses). Sweeping `MAX_EVALUE`:

| MAX_EVALUE | negatives.fa `pass` rate | shuffled-null `pass` rate | arab_ltr_all_clean `pass` rate | hardest real-gold cells (d=0.4-0.5, no flank) `pass` rate |
|---|---|---|---|---|
| 1e-2 | 10.53% | 0.019% | -- | -- |
| 1e-3 (old default) | 9.09% | 0.003% | 99.80% | 72.80% |
| 1e-4 | 7.87% | 0.000% | -- | -- |
| 1e-6 | 6.60% | 0.000% | -- | -- |
| 1e-8 | 5.79% | 0.000% | -- | -- |
| **1e-10 (new default)** | **5.23%** | **0.000%** | **99.67%** | **53.78%** |
| 1e-30 | 2.86% | -- | 90.19% | 20.96% |
| 1e-50 | 2.14% | -- | 74.90% | 11.36% |
| 1e-100 | 1.36% | -- | 55.35% | 3.99% |
| 1e-200 | 0.41% | -- | 31.14% | 2.58% |

**The shuffled null is clean at every threshold tested** (<=0.02% throughout, 0% at the
new default) -- confirming the residual `negatives.fa` false-"pass" rate is not chance
alignment noise near a significance boundary. It is real, strongly-significant direct
terminal-repeat structure inside specific TE subclasses, confirmed by inspecting the
records that still pass at extreme thresholds (e.g. `EnSpm-N1a_CR`, bitscore 757,
`k2p=0.0`, flank5_len=flank3_len=0 -- the entire 1,514 bp record is an exact 757 bp direct
repeat with nothing else in it; `TART-A`, a Drosophila telomeric retrotransposon with a
known, genuine terminal-repeat-like structure despite its non-LTR classification). Breaking
the residual down by class at an extreme threshold (1e-200) confirms it is concentrated,
not diffuse: Satellite entries pass at 6.7-37-100% depending on source library (tandem
repeats are close to tautologically "terminal-repeat-like"), tRNA/rRNA/Simple-Repeat
(small n, inherently repetitive) similarly high, Helitron/Neptune 1.6-3.6% (library-
consensus internal duplication), and a long tail of individual entries across many other
classes each under ~1.2% -- while the numerically dominant classes (hAT, Mariner/Tc1, L1,
MuDR, DNA, EnSpm/CACTA) are at or below ~0.1-0.8%. No `MAX_EVALUE` threshold can separate
these without the family classification this tool's own "Non-goals" (section 1) explicitly
declines to perform.

**Reaching <1% false-`pass` on `negatives.fa` is achievable only at unacceptable cost to
true positives**, and this is the operative finding: false-`pass` rate falls only slowly
as `MAX_EVALUE` tightens (9.09% -> 5.23% -> 2.86% -> ... -> 0.41% at 1e-200), while real-
data sensitivity collapses much faster and much earlier. On `arab_ltr_all_clean.fa.gz`
(10,307 real, structurally-supported LTR-RT calls -- exactly the tool's stated input
contract, not an edge case), `pass` rate is still 99.67% at `MAX_EVALUE=1e-10` but falls to
90.2% by `1e-30` (2.86% FP -- still nowhere near the 1% target) and to 55.3% by `1e-100`
(1.36% FP -- just barely above target). The cost concentrates earliest and worst in exactly
the population the tool is supposed to serve at the divergent end of its range: the
hardest real gold-perturbed cells (d_nominal in {0.4, 0.5}, no added flank -- correctly-
bounded, high-divergence, the lowest-bitscore-per-base population by construction) already
lose over a quarter of their `pass` calls between the old and new default (72.8% -> 53.8%),
and that population is reduced to noise (<4%) by `1e-100`. **`MAX_EVALUE = 1e-10` is the
new default**: it is the last point before the steep part of this cost curve (arab
`pass` rate moves only 13/10,307 records) while nearly halving the false-`pass` rate on
real non-LTR TEs (9.09% -> 5.23%, a genuine, substantial improvement). It does **not**
reach <1%, and no threshold does at acceptable cost; this is reported as a measured,
inherent limitation rather than forced past the point the data supports. Full sweep,
per-class breakdown and example records: `bench/out/memo_real.md`,
`bench/out/negatives_diagnosis.json`, `bench/out/maxevalue_sweep.json`,
`bench/out/maxevalue_truepos_sweep.json`.

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
