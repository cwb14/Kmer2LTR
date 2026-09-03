import random
import pytest
from kmer2ltr.align import discover, calibrate, terminal_snap, T_BITS
from kmer2ltr.scoring import GENERIC_MATRIX

def _rnd(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))

def _evolve(s, p, seed):
    r = random.Random(seed)
    ti = {"A": "G", "G": "A", "C": "T", "T": "C"}
    tv = {"A": "CT", "G": "CT", "C": "AG", "T": "AG"}
    out = []
    for c in s:
        if r.random() < p:
            out.append(ti[c] if r.random() < 2/3 else r.choice(tv[c]))
        else:
            out.append(c)
    return "".join(out)

def _bounds(S):
    h = discover(S, GENERIC_MATRIX)
    m, _, _ = calibrate(S, h)
    h2 = discover(S, m) or h
    return terminal_snap(S, h2, m)

def test_perfect_element_snaps_to_both_termini():
    ltr = _rnd(400, 1)
    S = ltr + _rnd(1200, 2) + ltr
    b = _bounds(S)
    assert b.l5b == 0
    assert b.l3e == len(S) - 1

def test_perfect_element_snaps_even_at_high_divergence():
    """The headline requirement. Plain SW trims the terminus on ~46% of these
    at p=0.25; Stage 3 must not."""
    hits = 0
    trials = 40
    for seed in range(trials):
        anc = _rnd(400, seed)
        l5 = _evolve(anc, 0.125, seed + 500)
        l3 = _evolve(anc, 0.125, seed + 900)
        S = l5 + _rnd(1200, seed + 300) + l3
        b = _bounds(S)
        if b.l5b == 0 and b.l3e == len(S) - 1:
            hits += 1
    assert hits >= int(0.9 * trials), f"only {hits}/{trials} snapped to both termini"

def test_real_flank_is_detected():
    ltr = _rnd(400, 3)
    f5, f3 = _rnd(60, 4), _rnd(90, 5)
    S = f5 + ltr + _rnd(1200, 6) + ltr + f3
    b = _bounds(S)
    assert b.l5b == pytest.approx(60, abs=5)
    assert (len(S) - 1 - b.l3e) == pytest.approx(90, abs=5)

def test_asymmetric_flank_only_five_prime():
    ltr = _rnd(400, 7)
    S = _rnd(70, 8) + ltr + _rnd(1200, 9) + ltr
    b = _bounds(S)
    assert b.l5b == pytest.approx(70, abs=5)
    assert b.l3e == len(S) - 1

def test_tiny_flank_below_floor_is_absorbed():
    """<5 bp flanks are information-theoretically undetectable; the spec
    commits to snapping rather than guessing."""
    ltr = _rnd(400, 10)
    S = "ACG" + ltr + _rnd(1200, 11) + ltr + "TT"
    b = _bounds(S)
    assert b.l5b <= 3 and (len(S) - 1 - b.l3e) <= 2

def test_margin_bits_is_larger_for_unambiguous_calls():
    ltr = _rnd(400, 12)
    clear = _rnd(300, 13) + ltr + _rnd(1200, 14) + ltr + _rnd(300, 15)
    b_clear = _bounds(clear)
    tiny = "AC" + ltr + _rnd(1200, 16) + ltr
    b_tiny = _bounds(tiny)
    assert b_clear.margin_bits is not None
    if b_tiny.margin_bits is not None:
        assert b_clear.margin_bits > b_tiny.margin_bits

def test_bounds_stay_ordered_and_in_range():
    ltr = _rnd(300, 17)
    S = _rnd(40, 18) + ltr + _rnd(800, 19) + ltr + _rnd(40, 20)
    b = _bounds(S)
    assert 0 <= b.l5b <= b.l5e < b.l3b <= b.l3e <= len(S) - 1


def test_external_credit_reaches_the_terminal_snap_entry_point():
    """`terminal_snap` is the Stage 1 entry point the benchmarks call; its
    `credit` has to arrive at `_extend` the same way `snap_bounds`' does."""
    ltr = _rnd(400, 71)
    S = _rnd(40, 72) + ltr + _rnd(1200, 73) + _evolve(ltr, 0.1, 74) + _rnd(40, 75)
    h = discover(S, GENERIC_MATRIX)
    m, _, _ = calibrate(S, h)
    h2 = discover(S, m) or h
    tight = terminal_snap(S, h2, m)
    loose = terminal_snap(S, h2, m, credit=1e6)
    assert tight.l5b > 0 and tight.l3e < len(S) - 1
    assert loose.l5b == 0 and loose.l3e == len(S) - 1
