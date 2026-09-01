"""Gold-subset selection and known-perturbation of real LTR-RT elements.

`simulate.py` builds elements from library CONSENSUS sequence. This
module instead starts from REAL elements the tool already resolves cleanly,
applies a KNOWN perturbation (mutation and/or flank), and hands the result to
`gold_robustness.py` to check whether the tool recovers the known-correct
answer. Real sequence carries real composition bias and internal structure a
simulator does not reproduce.

Circularity note (binding -- see docs/benchmarks.md): `select_gold` filters on
the tool's OWN output, so a subset built here must never be used to argue
"the tool calls boundaries correctly" -- that premise is exactly what is
under test. It is valid only for the before/after question: given an element
already called clean, does a KNOWN perturbation break it, and does the tool
recover the KNOWN-correct answer? The truth comes from the perturbation
(`perturb`'s returned truth dict), never from a second call to the tool.
Selection additionally requires signals the tool never uses when it makes its
own boundary call (canonical TG..CA termini, equal LTR lengths), which is
what keeps selection from simply rubber-stamping the tool's own opinion.
"""
from __future__ import annotations

import csv
import sys
from functools import lru_cache
from pathlib import Path
from typing import Callable, Iterator

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from kmer2ltr.align import classify      # noqa: E402
from kmer2ltr.fasta import read_fasta    # noqa: E402
from bench.simulate import evolve, shuffle_dinuc  # noqa: E402

# Indel rate used while perturbing gold LTRs. Matches `bench.simulate.simulate_element`'s own
# default (`indel_frac=0.10`) rather than introducing a second, untested free parameter for this benchmark alone.
INDEL_FRAC = 0.10


def select_gold(pred_tsv, fasta, min_bitscore: float = 200.0, max_k2p: float = 0.05,
                 require_motif: bool = True) -> Iterator[tuple[str, str]]:
    """Yield `(seq_id, sequence)` for elements the tool already resolves cleanly.

    ALL of the following must hold:
      - `status == "pass"`
      - `ltr5_start == 1` and `ltr3_end == seq_len` -- no overextension called,
        i.e. the tool did not report any 5' or 3' flank on the raw input.
      - `ltr5_len == ltr3_len` -- a signal the tool's boundary call never
        uses, so agreement with it is not automatic.
      - `k2p <= max_k2p` -- the two copies are still close enough that a
        perturbation grid reaching d=0.5 is a real stress test, not noise
        added on top of noise.
      - `bitscore >= min_bitscore` -- the pair is unambiguously significant.
      - if `require_motif`: the 5' LTR begins `TG` and the 3' LTR ends `CA`
        -- canonical termini, a second signal the tool's own call never uses.

    `pred_tsv` and `fasta` are joined by RECORD ORDER, not by id (ids may
    repeat in either file) -- the same contract `runner.py` documents for its
    own output ("streams records, preserves input order, one row per
    record"), which is what makes a positional zip here correct rather than
    merely convenient.

    Streaming throughout: one TSV row and one FASTA record in flight at a
    time, never the whole file.
    """
    with open(pred_tsv, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        records = read_fasta(fasta)
        for row, (seq_id, seq) in zip(reader, records):
            if row["status"] != "pass":
                continue
            seq_len = int(row["seq_len"])
            ltr5_start = int(row["ltr5_start"])
            ltr3_end = int(row["ltr3_end"])
            ltr5_len = int(row["ltr5_len"])
            ltr3_len = int(row["ltr3_len"])
            k2p = float(row["k2p"])
            bitscore = float(row["bitscore"])

            if ltr5_start != 1 or ltr3_end != seq_len:
                continue
            if ltr5_len != ltr3_len:
                continue
            if k2p > max_k2p:
                continue
            if bitscore < min_bitscore:
                continue
            if require_motif and (seq[:2] != "TG" or seq[-2:] != "CA"):
                continue
            yield seq_id, seq


@lru_cache(maxsize=8192)
def gold_ltr_len(seq: str) -> int:
    """Re-derive an already-gold element's LTR length from the tool's own
    (deterministic) boundary call.

    `perturb`'s signature (below) carries no boundary argument -- by design,
    so it must locate the LTR/
    internal split itself. Re-running `classify` on the untouched input is
    the only source of that split that cannot drift from what `select_gold`
    already required of it: `classify` is a pure function of `seq` (no
    randomness anywhere in discover/calibrate/terminal_snap/outermost), so
    calling it again on the SAME sequence that was already gold-selected is
    guaranteed to reproduce the identical boundary, not merely a similar one.

    Public (not `_`-prefixed) and memoized on `seq`: `perturb` uses it to
    split the element, and `gold_robustness.py` reuses the SAME cached call
    to carve a neighbour's internal region for `flank_source="other"`, so
    the cost is paid once per unique gold element rather than once per grid
    cell (7 x 11 x 2 = 154 cells) or twice per element (once per module).
    """
    r = classify("_gold_probe", seq)
    if (r.status != "pass" or r.ltr5_start != 1 or r.ltr3_end != len(seq)
            or r.ltr5_len != r.ltr3_len):
        raise ValueError(
            "perturb() requires an already gold-selected element (status=pass, "
            "no overextension, equal LTR lengths); got a sequence that does not "
            "satisfy that on re-classification -- was it really gold-selected?"
        )
    return r.ltr5_len


def perturb(seq: str, d: float, kappa: float, flank5: int, flank3: int, rng,
            flank_source: Callable[[int, object], str]) -> tuple[str, dict]:
    """Apply a KNOWN perturbation to a gold element and return the truth.

    A single ancestral copy of the element's own 5' LTR is evolved TWICE,
    independently, for a branch length of d/2 each (via `bench.simulate.
    evolve`, reused rather than reimplemented) -- exactly the construction
    `bench.simulate.simulate_element` already uses which was verified to land
    within 10% of nominal out to d=0.7. Using one ancestor rather than
    evolving the gold element's existing 5' and 3' copies independently from
    where they already sit matters: those two copies already differ by up to
    `max_k2p` (0.05 by default) from gold selection, so "evolve each one
    further by d/2" would not land at pairwise distance d -- it would land at
    d plus that leftover baseline, contaminating the one number (K2P bias vs
    KNOWN d) this benchmark exists to measure, worst at the most important
    grid cell, d=0.

    The internal region is sliced out and spliced back UNCHANGED. Flanks are
    obtained by calling `flank_source(n, rng)` once per side (only when
    n > 0), so this function is agnostic to how flanks are actually sourced
    -- `flank_from` below is the production source; the unit tests pass a
    trivial callable directly.

    Returns `(seq, truth)`; `truth` carries `ltr5_start`, `ltr5_end`,
    `ltr3_start`, `ltr3_end`, `d_nominal`, `kappa`, `flank5`, `flank3` --
    all 1-based inclusive coordinates into the returned `seq`, computed from
    the ACTUAL lengths of the evolved pieces (indels can change LTR length),
    never assumed.
    """
    ltr_len = gold_ltr_len(seq)
    total_len = len(seq)
    ltr_anc = seq[:ltr_len]
    internal = seq[ltr_len: total_len - ltr_len]

    l5 = evolve(ltr_anc, d / 2.0, kappa, INDEL_FRAC, rng)
    l3 = evolve(ltr_anc, d / 2.0, kappa, INDEL_FRAC, rng)

    f5 = flank_source(flank5, rng) if flank5 else ""
    f3 = flank_source(flank3, rng) if flank3 else ""

    out = f5 + l5 + internal + l3 + f3
    truth = {
        "ltr5_start": len(f5) + 1,
        "ltr5_end": len(f5) + len(l5),
        "ltr3_start": len(f5) + len(l5) + len(internal) + 1,
        "ltr3_end": len(out) - len(f3),
        "d_nominal": d,
        "kappa": kappa,
        "flank5": flank5,
        "flank3": flank3,
    }
    return out, truth


def flank_from(pool_seq: str, n: int, rng, mode: str) -> str:
    """Draw an n-bp composition-realistic flank from `pool_seq`.

    mode="other": a contiguous, randomly-placed window of real sequence
    (meant to be a DIFFERENT element's internal region) -- realistic in both
    composition and local order, since it is unmodified genomic sequence,
    just not this element's own flank.

    mode="shuffle": a mononucleotide shuffle of `pool_seq` (reuses
    `bench.simulate.shuffle_dinuc` rather than a second implementation) --
    the harder null. Matched to the element's OWN base composition exactly,
    so a detector that merely notices a compositional shift at the boundary
    cannot succeed against it; only a detector that notices the loss of
    homology can.

    `pool_seq` is tiled up to at least length n before drawing, so a flank
    request longer than the supplied pool still returns exactly n bases
    rather than raising or silently truncating.
    """
    if n <= 0:
        return ""
    if not pool_seq:
        raise ValueError("flank_from: pool_seq is empty")
    reps = -(-n // len(pool_seq))          # ceil division
    tiled = pool_seq * reps
    if mode == "shuffle":
        return shuffle_dinuc(tiled, rng)[:n]
    if mode == "other":
        start = rng.randrange(0, len(tiled) - n + 1)
        return tiled[start:start + n]
    raise ValueError(f"flank_from: unknown mode {mode!r}, expected 'other' or 'shuffle'")
