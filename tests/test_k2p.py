import math
import pytest
from kmer2ltr.k2p import SubstCounts, count_substitutions, k2p_distance, p_distance

def test_counts_basic():
    #        match  ts(A>G)  tv(A>C)  gap      N
    a = "AAAA" "A" "A" "A" "A"
    b = "AAAA" "G" "C" "-" "N"
    c = count_substitutions(a, b)
    assert (c.n_match, c.n_ts, c.n_tv, c.n_gapcols) == (4, 1, 1, 1)
    assert c.n_sites == 6
    assert c.aln_len == 8

def test_counts_all_four_transitions():
    a, b = "ACGT", "GTAC"       # A>G ts, C>T ts, G>A ts, T>C ts
    c = count_substitutions(a, b)
    assert (c.n_ts, c.n_tv) == (4, 0)

def test_counts_transversions():
    a, b = "AACC", "CCAA"       # A<->C are all transversions
    c = count_substitutions(a, b)
    assert (c.n_ts, c.n_tv) == (0, 4)

def test_identical_sequences_zero_distance():
    c = count_substitutions("ACGTACGT", "ACGTACGT")
    d, se = k2p_distance(c)
    assert d == 0.0 and se == 0.0
    assert p_distance(c) == 0.0

def test_k2p_matches_hand_computed_value():
    # P=0.10, Q=0.05, n=100 -> d = -0.5*ln(1-0.25) - 0.25*ln(1-0.10)
    c = SubstCounts(n_sites=100, n_match=85, n_ts=10, n_tv=5, n_gapcols=0, aln_len=100)
    expected = -0.5 * math.log(1 - 2*0.10 - 0.05) - 0.25 * math.log(1 - 2*0.05)
    d, se = k2p_distance(c)
    assert d == pytest.approx(expected, rel=1e-12)
    # exact SE value is checked in test_k2p_standard_error_matches_delta_method_closed_form
    assert se > 0

def test_k2p_undefined_when_saturated():
    # 1 - 2P - Q <= 0
    c = SubstCounts(n_sites=100, n_match=0, n_ts=50, n_tv=50, n_gapcols=0, aln_len=100)
    assert k2p_distance(c) == (None, None)

def test_k2p_undefined_with_no_sites():
    c = SubstCounts(n_sites=0, n_match=0, n_ts=0, n_tv=0, n_gapcols=10, aln_len=10)
    assert k2p_distance(c) == (None, None)
    assert p_distance(c) is None

def test_k2p_monotonic_in_divergence():
    ds = []
    for mism in (5, 10, 20, 30):
        c = SubstCounts(n_sites=100, n_match=100-mism, n_ts=mism*2//3,
                        n_tv=mism - mism*2//3, n_gapcols=0, aln_len=100)
        ds.append(k2p_distance(c)[0])
    assert ds == sorted(ds)

def test_k2p_exceeds_p_distance():
    # multiple-hit correction must inflate the raw proportion
    c = SubstCounts(n_sites=100, n_match=75, n_ts=17, n_tv=8, n_gapcols=0, aln_len=100)
    assert k2p_distance(c)[0] > p_distance(c)

def test_k2p_standard_error_matches_delta_method_closed_form():
    """Exact-value check on the SE, derived independently via the delta method.

    d = -0.5*ln(w1) - 0.25*ln(w2), so dd/dP = 1/w1 = a and
    dd/dQ = 0.5*(1/w1 + 1/w2) = b. Under multinomial sampling
    Var(P)=P(1-P)/n, Var(Q)=Q(1-Q)/n, Cov(P,Q)=-PQ/n, giving
    Var(d) = a^2*P + b^2*Q - (a*P + b*Q)^2, all over n.
    """
    n, ts, tv = 100, 10, 5
    c = SubstCounts(n_sites=n, n_match=n - ts - tv, n_ts=ts, n_tv=tv,
                    n_gapcols=0, aln_len=n)
    P, Q = ts / n, tv / n
    w1, w2 = 1 - 2 * P - Q, 1 - 2 * Q
    a = 1 / w1
    b = 0.5 * (1 / w1 + 1 / w2)
    expected_se = math.sqrt((a * a * P + b * b * Q - (a * P + b * Q) ** 2) / n)
    _, se = k2p_distance(c)
    assert se == pytest.approx(expected_se, rel=1e-12)
    assert se == pytest.approx(0.04633146812126295, rel=1e-12)   # pinned value

def test_k2p_standard_error_matches_monte_carlo_sampling_sd():
    """The reported SE must match the actual spread of d across replicates."""
    import random
    rng = random.Random(0)
    P_true, Q_true, n, reps = 0.10, 0.05, 500, 3000
    ds, ses = [], []
    for _ in range(reps):
        ts = tv = 0
        for _ in range(n):
            u = rng.random()
            if u < P_true:
                ts += 1
            elif u < P_true + Q_true:
                tv += 1
        d, se = k2p_distance(SubstCounts(n_sites=n, n_match=n - ts - tv, n_ts=ts,
                                         n_tv=tv, n_gapcols=0, aln_len=n))
        if d is not None:
            ds.append(d)
            ses.append(se)
    mean = sum(ds) / len(ds)
    empirical_sd = (sum((x - mean) ** 2 for x in ds) / (len(ds) - 1)) ** 0.5
    reported_se = sum(ses) / len(ses)
    assert reported_se == pytest.approx(empirical_sd, rel=0.05)


def test_insertion_time_divides_by_two_branches():
    """The two LTRs are identical at insertion and diverge independently, so
    the divergence spans TWO branches: t = d / (2*mu), not d / mu. Getting this
    wrong is the classic factor-of-two error in LTR dating."""
    from kmer2ltr.k2p import insertion_time
    assert insertion_time(0.14, 7e-9) == 10_000_000
    assert insertion_time(0.0, 7e-9) == 0


def test_insertion_time_halves_when_the_rate_doubles():
    from kmer2ltr.k2p import insertion_time
    assert insertion_time(0.1, 2e-8) * 2 == insertion_time(0.1, 1e-8)


def test_insertion_time_is_undefined_without_both_a_distance_and_a_rate():
    """A saturated pair has no distance and an unspecified species has no rate;
    inventing a default rate would report an age nobody asked for."""
    from kmer2ltr.k2p import insertion_time
    assert insertion_time(None, 7e-9) is None
    assert insertion_time(0.1, None) is None
    assert insertion_time(0.1, 0.0) is None
    assert insertion_time(0.1, -1e-9) is None
