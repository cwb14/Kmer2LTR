import random
import pytest
from ltrk2p.align import discover, Hit, GAP_OPEN, GAP_EXTEND
from ltrk2p.scoring import GENERIC_MATRIX

def _rnd(n, seed=None):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))

def test_finds_exact_ltr_pair_with_no_flank():
    ltr = _rnd(300, 1); internal = _rnd(1200, 2)
    S = ltr + internal + ltr
    h = discover(S, GENERIC_MATRIX)
    assert h.qb == 0 and h.qe == 299
    # 3' LTR occupies the last 300 bases
    assert len(S) - h.w + h.rb == len(S) - 300
    assert len(S) - h.w + h.re == len(S) - 1

def test_finds_pair_when_flanked_on_both_sides():
    ltr = _rnd(300, 3); internal = _rnd(1200, 4)
    f5, f3 = _rnd(80, 5), _rnd(120, 6)
    S = f5 + ltr + internal + ltr + f3
    h = discover(S, GENERIC_MATRIX)
    assert h.qb == 80 and h.qe == 379
    assert len(S) - h.w + h.re == len(S) - 121
    # 3' LTR starts at absolute position 80+300+1200=1580
    assert len(S) - h.w + h.rb == 1580

def test_window_grows_when_ltr_exceeds_initial_window():
    """LTR of 2500 bp exceeds w0=1500, so the window must double to find it."""
    ltr = _rnd(2500, 7); internal = _rnd(3000, 8)
    S = ltr + internal + ltr
    h = discover(S, GENERIC_MATRIX)
    assert h.w > 1500
    assert h.qb == 0 and h.qe == 2499

def test_window_never_exceeds_half_length():
    ltr = _rnd(400, 9)
    S = ltr + ltr                     # no internal region at all
    h = discover(S, GENERIC_MATRIX)
    assert h.w == len(S) // 2
    assert h.qb == 0 and h.qe == 399

def test_ltr5_end_always_before_ltr3_start():
    for seed in range(10):
        ltr = _rnd(200, seed); internal = _rnd(400, seed + 100)
        S = ltr + internal + ltr
        h = discover(S, GENERIC_MATRIX)
        ltr3_start_abs = len(S) - h.w + h.rb
        assert h.qe < ltr3_start_abs

def test_returns_none_for_sequence_too_short():
    assert discover("ACGT", GENERIC_MATRIX) is None

def test_random_sequence_yields_low_score():
    """No terminal repeat -> whatever is found must score far below a real pair."""
    S = _rnd(3000, 42)
    h = discover(S, GENERIC_MATRIX)
    ltr = _rnd(300, 1)
    real = discover(ltr + _rnd(1200, 2) + ltr, GENERIC_MATRIX)
    assert h is None or h.score < real.score / 4

def test_window_grows_for_ltr_at_and_above_twice_initial_window():
    """LTRs >= 2*w0 leave the two windows covering disjoint parts of the LTR,
    so the hit is pure noise and never touches an edge. Growth must still fire."""
    for ltr_len in (3000, 4000, 6000):
        ltr = _rnd(ltr_len, ltr_len)
        S = ltr + _rnd(5000, ltr_len + 1) + ltr
        h = discover(S, GENERIC_MATRIX)
        assert h is not None, f"no hit at ltr_len={ltr_len}"
        assert h.qb == 0 and h.qe == ltr_len - 1, \
            f"ltr_len={ltr_len}: got qb={h.qb} qe={h.qe}, want 0..{ltr_len-1}"

def test_degenerate_inputs_do_not_raise():
    for S in ("", "ACGT", "N" * 500, "A" * 500, _rnd(150, 1)):
        discover(S, GENERIC_MATRIX)      # must not raise

def test_non_ltr_sequence_returns_none_not_a_noise_hit():
    """Random sequence has no terminal repeat. Returning a Hit anyway would hand
    callers noise they cannot distinguish from a real pair -- Hit has no
    significance field. Measured E-values on such hits were 58-467 against a
    1e-3 threshold."""
    for length in (3000, 9000):
        for seed in range(5):
            assert discover(_rnd(length, seed), GENERIC_MATRIX) is None, \
                f"spurious hit on random {length}bp seed={seed}"
