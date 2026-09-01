import math, random
import pytest
from kmer2ltr.k2p import count_substitutions, k2p_distance

def _reference_k80(a, b):
    """Independent K80 implementation, written from the 1980 paper formulas
    without reference to k2p.py, over the same site-inclusion rule."""
    ts = tv = n = 0
    ring = {"A": 0, "G": 0, "C": 1, "T": 1}
    for x, y in zip(a, b):
        if x in "ACGT" and y in "ACGT":
            n += 1
            if x != y:
                ts += 1 if ring[x] == ring[y] else 0
                tv += 0 if ring[x] == ring[y] else 1
    p, q = ts / n, tv / n
    return 0.5 * math.log(1 / (1 - 2 * p - q)) + 0.25 * math.log(1 / (1 - 2 * q))

def test_matches_independent_implementation_on_random_alignments():
    random.seed(0)
    for _ in range(200):
        n = random.randint(50, 400)
        a = "".join(random.choice("ACGT") for _ in range(n))
        b = "".join(random.choice("ACGT") if random.random() < 0.25 else c for c in a)
        d, _ = k2p_distance(count_substitutions(a, b))
        if d is not None:
            assert d == pytest.approx(_reference_k80(a, b), rel=1e-10)
