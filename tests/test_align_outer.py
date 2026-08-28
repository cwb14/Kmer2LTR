import random
from ltrk2p.align import discover, calibrate, terminal_snap, outermost
from ltrk2p.scoring import GENERIC_MATRIX

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

def test_retained_nested_element_reports_outer_pair():
    """Old outer element (15% divergent LTRs) with a young nested element
    (1% divergent LTRs) left inside. The nested pair scores higher, so plain
    discovery picks it; Stage 4 must recover the outer pair."""
    outer = _rnd(500, 1)
    inner = _rnd(400, 2)
    nested = inner + _rnd(600, 3) + _evolve(inner, 0.01, 4)
    S = outer + _rnd(400, 5) + nested + _rnd(400, 6) + _evolve(outer, 0.15, 7)
    h = discover(S, GENERIC_MATRIX)
    m, _, _ = calibrate(S, h)
    b = outermost(S, terminal_snap(S, h, m), m)
    assert b.l5b < 50, f"expected outer LTR at ~0, got {b.l5b}"
    assert len(S) - 1 - b.l3e < 50

def test_does_not_fire_when_boundaries_already_terminal():
    ltr = _rnd(400, 8)
    S = ltr + _rnd(1000, 9) + ltr
    h = discover(S, GENERIC_MATRIX)
    m, _, _ = calibrate(S, h)
    b0 = terminal_snap(S, h, m)
    assert outermost(S, b0, m) == b0

def test_does_not_fire_on_genuine_flank():
    """Real flanking DNA has no terminal repeat, so nothing outer exists."""
    ltr = _rnd(400, 10)
    S = _rnd(200, 11) + ltr + _rnd(1000, 12) + ltr + _rnd(200, 13)
    h = discover(S, GENERIC_MATRIX)
    m, _, _ = calibrate(S, h)
    b0 = terminal_snap(S, h, m)
    assert outermost(S, b0, m) == b0

def test_does_not_invent_an_outer_pair_from_unrelated_flanks():
    """Random flanks share no homology. Adopting a chance match as an 'outer
    pair' corrupts a correctly-called element -- worse than this stage not
    existing. Leaked 14/2500 when the gate was scored with the calibrated
    matrix instead of GENERIC_MATRIX."""
    spurious = 0
    trials = 0
    for div in (0.01, 0.05):
        for flank in (300, 600):
            for seed in range(25):
                trials += 1
                anc = _rnd(400, seed * 7)
                S = (_rnd(flank, seed + 1) + _evolve(anc, div / 2, seed + 2)
                     + _rnd(1200, seed + 3) + _evolve(anc, div / 2, seed + 4)
                     + _rnd(flank, seed + 5))
                h = discover(S, GENERIC_MATRIX)
                if h is None:
                    continue
                m, _, _ = calibrate(S, h)
                h2 = discover(S, m) or h
                b0 = terminal_snap(S, h2, m)
                b = outermost(S, b0, m)
                if (b0.l5b, b0.l5e, b0.l3b, b0.l3e) != (b.l5b, b.l5e, b.l3b, b.l3e):
                    spurious += 1
    assert spurious == 0, f"invented an outer pair on {spurious}/{trials} unrelated-flank elements"
