import random
from kmer2ltr.align import discover, calibrate, terminal_snap
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

def test_stage3_beats_raw_sw_on_terminal_accuracy():
    """Raw SW trims the terminus on ~46% of perfect elements at p=0.25
    (measured). Stage 3 must cut that by at least 5x."""
    raw_bad = snap_bad = 0
    n = 60
    for seed in range(n):
        anc = _rnd(400, seed)
        S = _evolve(anc, 0.125, seed+11) + _rnd(1000, seed+22) + _evolve(anc, 0.125, seed+33)
        h = discover(S, GENERIC_MATRIX)
        if h is None:
            continue
        if h.qb != 0 or (len(S) - h.w + h.re) != len(S) - 1:
            raw_bad += 1
        m, _, _ = calibrate(S, h)
        b = terminal_snap(S, discover(S, m) or h, m)
        if b.l5b != 0 or b.l3e != len(S) - 1:
            snap_bad += 1
    assert snap_bad * 5 <= raw_bad or snap_bad == 0, \
        f"raw SW wrong on {raw_bad}/{n}, Stage 3 wrong on {snap_bad}/{n}"
