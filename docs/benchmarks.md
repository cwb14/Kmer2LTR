# Benchmarks

Every default in `ltrk2p` is set from measurement rather than by hand. This
document records what was measured, on what, and what it decided. `docs/design.md`
gives the reasoning; this gives the numbers.

Three properties are asked of every benchmark here:

- **The truth must not come from the tool.** Where it does — one benchmark
  deliberately selects elements the tool already resolves cleanly, in order to
  ask a different question — that is stated and its numbers are not used for any
  claim about accuracy on arbitrary input.
- **A configuration that answers less often must not look more accurate.**
  Boundary error is reported over located records with the loss rate printed
  beside it; called-flank error counts a missed flank as a called length of zero.
- **The evidence must be reproducible from the repository.** Every table below
  comes from a script in `bench/`.

---

## The datasets

| | records | what it is |
|---|---|---|
| `truth.fa` | 31,315 | `X-LTR` + `X-I` + `X-LTR` concatenations built from Repbase, Dfam, MTEC and riceTElib entries that split a family into LTR and internal parts. Boundaries are known **exactly** because the element is a construction, not an annotation. |
| homology grid | 70,000 | perfect elements from `truth.fa` under known substitutions, flanks and indels |
| gold-perturbed grid | 260,876 | 1,694 real elements under known perturbation |
| `negatives.fa` | 46,823 | real non-LTR transposons: DNA transposons, LINEs, SINEs, Helitrons, satellites |
| shuffled null | 31,315 | composition-matched mononucleotide shuffles of `truth.fa` |
| real call sets | 256,776 | arabidopsis, poa, human LTR-RT calls plus four repeat libraries |

### Homology grid — `bench/homology_grid.py`

The primary benchmark, and the only one whose element selection involves no
motif, no target-site duplication and no call by `ltrk2p`.

Its elements come from `truth.fa` and are **perfect by construction**: the two
LTR copies are literally the same string, so their divergence is exactly zero and
their boundaries exactly known. Everything measured is a deviation the harness
itself introduced.

Three axes, 100 cells per element:

- **substitutions**, parameterised by target *pairwise p-distance*
  {0, 5, 10, 15, 20, 25, 30, 35}% rather than by branch length, so "25%
  mutations" means 25% of sites differ. The K2P distance delivering it is solved
  for by bisection.
- **flanks** of {0, 5, 25, 45, 65} bp, drawn from a composition-matched shuffle
  of the element's own sequence — a detector cannot succeed by noticing a
  compositional step.
- **indels**, swept independently at {0.001, 0.002, 0.005, 0.01, 0.02, 0.05}
  events per site per branch, geometric lengths.

**The true alignment is tracked through evolution**, so realized p-distance and
realized K2P are measurements rather than approximations even under indels. This
is what `bench/simulate.py` cannot do: it computes its realized K2P by
truncate-and-compare on unaligned strings, which disagrees with a proper
alignment even at zero indel rate. Nothing here depends on it.

Two panels are reported separately: `lib` (library constructions — the
non-circular one, and the deciding evidence) and `real` (motif-free arabidopsis
elements — real composition and internal structure, but selected from a
prediction file, so not independent).

### Gold-perturbed grid — `bench/gold_subset.py`, `bench/gold_robustness.py`

Real elements the tool already resolves cleanly, perturbed by a known amount.
Selection uses the tool's own output *and* requires canonical `TG`…`CA` termini,
so **this benchmark cannot show that boundaries are called correctly on arbitrary
input** — that premise is what selection assumes. It answers a narrower question
it is valid for: given an element already resolved cleanly, does a known
perturbation break it, and does the known-correct answer come back? The truth
comes entirely from the perturbation, never from a second call to the tool.

Half its flanks are lifted from a different element's internal region — real
sequence with real repeat structure — which makes it a harder false-flank test
than the homology grid's shuffles, and explains most of the disagreement between
the two.

---

## Where the defaults come from

### Boundary model selection is the load-bearing stage

Plain Smith-Waterman trims the true terminus whenever the terminal bases are
diverged, because ending one base earlier scores better. Per-element error rates
on perfectly bounded synthetic elements:

| p-distance | naive local aligner wrong | with Stage 3 |
|---|---|---|
| 0.05 | 14/80 | 3/80 |
| 0.15 | 37/80 | 2/80 |
| 0.25 | **63/80 (79%)** | **1/80 (1.2%)** |
| 0.35 | 69/80 | 2/80 |
| 0.45 | 75/80 | 3/80 |

Without this stage the tool would falsely report overextension on roughly three
quarters of perfectly bounded input at its target divergence. It is the reason
the tool exists in this form.

### Scoring: calibrated log-odds beats fixed penalties

Pooled false-flank rate on the gold grid: BLASTN-style 1/−3 gives 46.4%, a flat
+1/−1 gives 28.6%, the per-element calibrated matrix gives 20.2%. The 1/−3 scheme
crosses 50% between d=0.2 and d=0.3, matching its predicted ~75%-identity
break-even. The raw match/mismatch magnitude matters more than the
transition/transversion distinction, but both help.

### Gap penalties are calibrated per element

Deriving affine penalties from the core alignment's own indel statistics, rather
than fixing them at 6 bits to open and 2 to extend. At **matched false-flank
rate** — the only fair comparison, since the change moves the whole ROC:

| | homology grid | gold grid |
|---|---|---|
| smallest flank (5 / 10 bp) | **+10.5 pts** | **+14.0 pts** |
| mid flank (25 / 20 bp) | **+7.5 pts** | **+4.0 pts** |
| large flank (45 / 50 bp) | **+2.1 pts** | −1.7 pts |
| largest (65 / 100 bp) | **+0.7 pts** | −1.2 pts |
| pairs lost | **−0.4 pts** | +0.5 pts |

It dominates outright on the non-circular grid and cuts K2P RMSE 35% under
indels. On the gold grid it trades 1.7 points of 50 bp detection for 14 points at
10 bp. Runtime cost: +2.6%.

### The flank threshold is divergence-aware

A single constant fits badly: at fixed `t_bits` the false-flank rate swings
40–80× across the observed divergence range while the threshold does not move.
The schedule is derived by `bench/calibrate_flank_threshold.py` under a rule
fixed before the sweep it consumes is run, with an explicit detection floor —
a threshold that reaches a false-flank target by abandoning large-flank detection
is inadmissible.

```
T_BITS_SCHEDULE = ((0.025, 2), (0.15, 8), (float("inf"), 10))
```

The shape is the opposite of the obvious guess. Tightening at high divergence is
what one reaches for, but the floor rules it out: `t=15` and `t=20` cost more
large-flank detection than the floor permits, in *every* divergence bin. What
survives is a relaxation at low divergence, where false flanks are cheap and
boundaries unambiguous — in the lowest bin, `t=2` lifts 5 bp flank detection from
83.7% to 98.9% for 1.5 points of false-flank rate.

Against a flat threshold re-tuned to the same pooled false-flank rate the
schedule is worth +3.7 points of det@10 on the gold grid and +2.6 of det@5 on the
homology grid, costing 0.0–1.4 points at every larger flank length. A real but
modest gain.

**Reference curves, measured at the shipped `MAX_EVALUE`.** Note the homology
curve is not monotone: `t=30` and `t=20` reach the same false-flank rate, but
`t=30` detects a third as many 25 bp flanks and loses 47% more pairs outright, so
the top of the range is simply dominated.

| t | gold ff | det10 | det50 | | homology ff | det5 | det45 |
|---|---|---|---|---|---|---|---|
| 2 | 0.3879 | 0.8378 | 0.8719 | | 0.2609 | 0.8135 | 0.9808 |
| 5 | 0.1722 | 0.7129 | 0.8699 | | 0.0368 | 0.4580 | 0.9762 |
| 8 | 0.0910 | 0.5128 | 0.8642 | | 0.0048 | 0.1730 | 0.9637 |
| 10 | 0.0616 | 0.3793 | 0.8562 | | 0.0020 | 0.0788 | 0.9455 |
| 15 | 0.0272 | 0.1300 | 0.8090 | | 0.0010 | 0.0035 | 0.8508 |
| 20 | 0.0169 | 0.0226 | 0.7218 | | 0.0008 | 0.0003 | 0.7090 |
| 30 | 0.0117 | 0.0075 | 0.5157 | | 0.0008 | 0.0003 | 0.4377 |

### The significance threshold, and a limit that cannot be tuned away

`MAX_EVALUE = 1e-10`, calibrated against `negatives.fa` and the shuffled null.

The shuffled null is clean at every threshold tested (≤0.02%, and 0.00% at the
shipped value), which proves the residual false-positive rate on `negatives.fa`
is not chance alignment noise. It is real terminal-repeat structure inside
specific subclasses: satellites are close to tautologically terminal-repeat-like,
and a handful of library-consensus Helitron and CACTA entries are near-exact
tandem constructions. Inspecting survivors at extreme thresholds confirms it —
one entry is a 1,514 bp record consisting of an exact 757 bp direct repeat and
nothing else.

**Reaching <1% false positives is achievable only at unacceptable cost**, and
that is the operative finding rather than a shortfall. False positives fall
slowly as the threshold tightens (9.09% → 5.23% → 2.86% → 0.41% at 1e-200) while
sensitivity on real input collapses far faster: arabidopsis pass rate is 99.67%
at 1e-10, 90.2% at 1e-30, 55.3% at 1e-100. No threshold separates these without
the family classification the tool explicitly declines to perform.

### Terminal trimming is not applied

Trimming the outermost few bases before counting substitutions was measured at
0, 3, 5 and 10 bp per end. In the correctly-bounded population RMSE is flat to
three significant figures at every depth and bias does not improve — exactly the
"costs data, removes no bias" prediction. In the flank-called subset bias does
move toward zero (a real ~32% reduction) but RMSE does not, because the variance
cost of losing sites cancels it. Adopting a trim requires reducing RMSE, not just
bias; neither population clears that bar.

---

## Alternatives that were measured and rejected

Recorded because a design is only as trustworthy as the things it declined.

**A graded flank extension**, replacing the binary snap with "extend to the
furthest endpoint whose running score is still above the noise floor". The
literal version is *provably* a no-op: discovery's core is the Smith-Waterman
optimum over the same windows, so a positive-scoring outward extension would
contradict that optimality. Confirmed at 0 of 1,194 records. The thresholded
version is measurably worse — boundary MAE 3.29 → 8.81, flank MAE 4.70 → 12.63,
K2P bias +0.0006 → +0.0168 — because a budget spent on the snap decision and
again on the endpoint lets the boundary creep into non-homologous sequence. On a
control with 8,000 bp of random flank it calls 7,994.

**Joint inner boundaries.** The two inner boundaries are by-products of the
opposite terminus's snap; re-deriving both from one alignment anchored at the
settled outer ends changes 418 of 62,644 records, so it is a real alternative.
It is right more often than wrong (276 vs 140) by *less* than it is wrong by:
mean summed inner error on the records it touches 242 → 557, pooled boundary MAE
unchanged at 3.29 → 3.29.

That null result is the informative one. The inner-boundary error is real — the
inner ends are called *too long*, mean +6.3 bp at d=0.1 rising to +14.6 at d=0.5,
which is the opposite sign from the greedy-trimming story one would expect — but
deriving them jointly and optimally leaves the pooled error unchanged to four
significant figures. It is estimation error under the model, not a fixable
artifact of how the ends are obtained.

**Composition estimated from the LTR core** rather than the whole element. Worse
on every panel of both grids: boundary MAE 3.29 → 3.79, 25 bp detection 79.2% →
74.1%, false-flank 0.20% → 0.41%. The log-odds denominator is the *null* a match
is judged against, and the null for "does homology continue past this boundary?"
is the surrounding sequence, not the repeat itself.

**A uniform-penalty aligner for the final alignment** was compared against the
calibrated matrix. The calibrated matrix roughly doubles K2P bias magnitude in
both populations with no reliable RMSE gain, so wavefront alignment stays
everywhere. It is also 100–500× faster on the high-identity bulk of real data and
still 1.1–1.9× faster at 35% divergence, and returns identical optimal scores
under matched penalties.

---

## The short-flank detection floor

A flank of length *k* can supply at most `k × α` bits of evidence that it is a
flank, where α is the expected per-base cost of aligning unrelated sequence under
the element's own matrix. Whenever `k × α < t_bits` the snap test cannot fire
whatever the sequence says: the boundary is decided before it is looked at.

α is read off the calibrated matrix and falls with divergence:

| d | t_bits demanded | α (bits/bp) | detection floor |
|---|---|---|---|
| 0.00 | 2.0 | 4.24 | 0.5 bp |
| 0.05 | 8.0 | 2.54 | 3.1 bp |
| 0.10 | 8.0 | 1.84 | 4.3 bp |
| 0.20 | 10.0 | 1.19 | 8.4 bp |
| 0.35 | 10.0 | 0.73 | 13.7 bp |

Robust to composition — at a 65/35 AT skew the floors move to 0.5 / 4.3 / 12.4
bp. This floor reproduces the shape of the measured detection curve, and turns
"short flanks are undetectable in principle" from an observation into an
equation.

`--flank-sensitivity` caps the demand at `β × α × k`. Measured on the homology
grid, with false-flank rate re-weighted to the **real** divergence distribution
rather than the grid's:

| setting | β | real-weighted ff | det5 | det25 | flank MAE | mean abs K2P bias | K2P RMSE |
|---|---|---|---|---|---|---|---|
| `strict` | off | 0.0057 | 0.303 | 0.911 | 3.69 | 0.00374 | 0.01547 |
| `balanced` | 1.0 | 0.0782 | 0.562 | 0.918 | 3.38 | **0.00237** | 0.01400 |
| `sensitive` | 0.5 | 0.1065 | 0.793 | 0.949 | 3.18 | 0.00292 | 0.01360 |

The two β values are not fitted. `β = 1.0` tests the observed extension cost
against its own expectation under the flank hypothesis, which is the unbiased
test and lands at roughly even odds on a genuine short flank. `β = 0.5` puts the
boundary midway between the two hypotheses — a homologous continuation scores
about 0, a non-homologous one about `−k·α` — which is the maximum-likelihood
split for a segment of that length. Values below 0.5 have no such derivation and
are simply permissive; they also turn K2P bias back upward.

Below the floor no threshold recovers certainty: a 1 bp flank carries one base of
evidence and is near-chance however it is tested. The cap raises detection to the
information-theoretic ceiling, not past it.

**Re-weighting matters.** The grids spread elements nearly evenly across
divergence; real input does not, so a grid-pooled false-flank rate understates
what a user sees in the low-divergence bins where most real elements live:

| | d̂ 0 | 0.05 | 0.1 | 0.2 | 0.3 | 0.4 |
|---|---|---|---|---|---|---|
| real (pooled, n=102,300) | 0.509 | 0.281 | 0.113 | 0.075 | 0.019 | 0.003 |
| gold grid | 0.173 | 0.163 | 0.156 | 0.160 | 0.161 | 0.138 |

The distribution is strongly dataset-dependent — arabidopsis and poa put 59–64%
of elements in the lowest bin, human only 4.7%, because most human ERV activity
is ancient. Any pooled rate should be read with that in mind.

---

## End-to-end validation

Run with no flags, so every value comes from the shipped source.

**Row-count invariant.** Output body lines equal input records on all seven real
datasets — 10,307 / 75,180 / 16,336 / 124,517 / 26,292 / 1,517 / 2,627, all
exact. A mismatch would be a bug, not a rounding detail.

**Status breakdown.**

| dataset | n | pass | weak_pair | no_pair | pass % |
|---|---|---|---|---|---|
| arabidopsis | 10,307 | 10,260 | 36 | 11 | 99.54% |
| poa | 75,180 | 75,034 | 90 | 56 | 99.81% |
| human | 16,336 | 16,164 | 97 | 75 | 98.95% |
| repbase | 124,517 | 3,986 | 226 | 119,486 | 3.20% |
| Dfam | 26,292 | 1,825 | 73 | 23,927 | 6.94% |
| MTEC maize | 1,517 | 618 | 1 | 836 | 40.74% |
| rice | 2,627 | 89 | 12 | 2,520 | 3.39% |
| **negatives.fa** | 46,823 | 2,406 | 200 | 44,217 | **5.14%** |
| **shuffled null** | 31,315 | 0 | 0 | 31,315 | **0.00%** |

The three putative-intact call sets pass at 99–100%, as they should. The four
general repeat libraries pass at 3–41%, correctly low, since they are not
filtered to LTR-class entries.

**Independent confirmation from signals the tool never uses.** Neither `TG`…`CA`
termini nor target-site duplications appear anywhere in `classify()`, so they
cannot be gamed and both are usable as external accuracy proxies on real data
where no truth exists.

| dataset | 5′TG at flank=0 | raw-ends baseline | TSD k=5 | shifted control | enrichment |
|---|---|---|---|---|---|
| arabidopsis | 0.666 | 0.605 | 0.0443 | 0.0053 | **8.4×** |
| human | 0.327 | 0.297 | 0.0221 | 0.0019 | **11.6×** |
| poa | 0.791 | 0.724 | 0.0277 | 0.0034 | **8.1×** |
| MTEC maize | 0.902 | 0.353 | — | — | — |

Boundaries called perfectly bounded land on a canonical terminus more often than
boundaries called flank-corrected, and above the raw-input-ends baseline, on
every genomic dataset. `repbase` reverses the pattern, which connects back to the
negative-control finding: it is not LTR-filtered, so much of its small passing
population is exactly the non-LTR-but-genuinely-repetitive material significance
alone cannot exclude, which has no reason to carry a retroviral terminus.

---

## Reproducing

```bash
python bench/build_truth.py --help          # ground truth from split library entries
python bench/homology_grid.py build --help  # the homology grid
python bench/gold_robustness.py --help      # the gold-perturbed grid
python bench/run_bench.py --help            # ablation grid and scoring
python bench/run_configs.py --help          # named pipeline configurations
python bench/calibrate_flank_threshold.py --help
python bench/report.py --help               # comparison tables
```

SLURM launchers for each are in `bench/`. They assume 20 cores and write to
`bench/out/`, which is gitignored.
