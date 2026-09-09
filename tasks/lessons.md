# Lessons

Patterns to apply before starting work in this repo. Reviewed at the top of any
non-trivial change.

## 2026-09-09 — Never design a repeat-finding mechanism that consumes the sequence

**What happened.** I proposed enumerating alternative repeat pairs by iterated
masked Smith-Waterman: find a pair, overwrite both copies with `N`, re-run. The
user replaced it with diagonal suppression.

**Why the correction is right.** Masking is per-base, so it removes the sequence
from *every* register at once, including the true register that needs those same
bases to score. In the shifted-diagonal case the wrong pair sits on the exact
sequence the right pair must use. The mechanism was self-defeating.

**Rule for next time.** When the competing hypotheses are *registers* over the
same sequence, suppress in offset space, never in sequence space. Ask "what does
this suppression cost the alternative I am trying to find?" before proposing it.

## 2026-09-09 — Reproduce a supplied baseline before building on it

**What happened.** The user supplied a fixture set with a `kmer2ltr_current_call`
column. I ran the tool against it before writing code. Zero of forty rows
reproduced, and several were impossible outputs: negative starts, ends past the
sequence length, `ltr5_end == ltr3_start`. Five of twenty "controls it currently
gets right" were also wrong. Building against that column would have measured a
fiction and hidden five real failures.

**Rule for next time.** A supplied truth column and a supplied *baseline* column
are different things. Verify the truth is internally consistent, and recompute
the baseline from the code every time. Any benchmark driver I write must
recompute rather than read a recorded baseline.

## 2026-09-09 — Validate on real fixtures before writing synthetics

**What happened.** My first synthetic reproduced the wrong failure. It placed the
competing pair strictly *inside* the LTR pair, which Stage 4 already recovers, so
the test passed for the wrong reason and proved nothing about the new rule. The
real records showed the competing pair *overlaps* both LTRs, which is exactly why
Stage 4 misses them. Rebuilding the synthetic to the measured geometry made it a
real test.

**Rule for next time.** Derive the synthetic's geometry from the measured real
case, not from my model of the case. Then assert the synthetic actually fails
under the old behaviour, so a test that stops reproducing the bug fails loudly
instead of passing silently.

## 2026-09-09 — Measure the diagnostic before writing it in the README

**What happened.** I wrote that the failure is recognisable because `motif` is not
`tg...ca`. Measured: 0 of 20 correct calls show `tg...ca` either. The claim was
worthless. The symptom that does hold is that the pair stops short of the record's
ends, 23 of the 25 wrong records against 1 false alarm in 15.

**Rule for next time.** Every user-facing diagnostic, symptom or rule of thumb
gets measured on the data before it is written down. Plausible-sounding biology
is not evidence.
