import random
import pytest
from ltrk2p.align import classify
from ltrk2p.fasta import sanitize

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

def test_identical_ltrs_give_zero_k2p_and_exact_bounds():
    ltr = _rnd(400, 1)
    S = ltr + _rnd(1200, 2) + ltr
    r = classify("x", S)
    assert r.status == "pass"
    assert r.ltr5_start == 1                 # 1-based
    assert r.ltr3_end == len(S)              # 1-based inclusive
    assert r.k2p == 0.0
    assert r.ltr5_len == 400 and r.ltr3_len == 400
    assert r.cigar == "400="

def test_coordinates_are_one_based_inclusive_and_slice_correctly():
    ltr = _rnd(300, 3)
    S = _rnd(50, 4) + ltr + _rnd(900, 5) + ltr + _rnd(60, 6)
    r = classify("x", S)
    assert S[r.ltr5_start - 1: r.ltr5_end] == ltr
    assert S[r.ltr3_start - 1: r.ltr3_end] == ltr
    assert r.flank5_len == r.ltr5_start - 1
    assert r.flank3_len == r.seq_len - r.ltr3_end

def test_k2p_recovers_simulated_divergence():
    anc = _rnd(2000, 7)
    S = _evolve(anc, 0.10, 8) + _rnd(1500, 9) + _evolve(anc, 0.10, 10)
    r = classify("x", S)
    assert r.status == "pass"
    assert r.k2p == pytest.approx(0.20, abs=0.04)   # two branches of 0.10 each

def test_random_sequence_reports_no_pair():
    r = classify("x", _rnd(3000, 11))
    assert r.status == "no_pair"
    assert r.k2p is None and r.ltr5_start is None

def test_solo_ltr_reports_no_pair():
    r = classify("x", _rnd(400, 12))
    assert r.status == "no_pair"

def test_too_short_reports_too_short():
    r = classify("x", "ACGT" * 10)
    assert r.status == "too_short"

def test_all_N_reports_all_ambiguous():
    r = classify("x", "N" * 2000)
    assert r.status == "all_ambiguous"

def test_sanitised_iupac_does_not_crash_and_excludes_N_sites():
    ltr = _rnd(400, 13)
    dirty = ltr[:100] + "RYKM" + ltr[104:]
    S = sanitize(dirty + _rnd(1000, 14) + ltr)
    r = classify("x", S)
    assert r.status == "pass"
    assert r.n_sites <= r.aln_len

def test_cs_flag_switches_alignment_string():
    ltr = _rnd(400, 15)
    S = ltr + _rnd(1000, 16) + _evolve(ltr, 0.05, 17)
    r_cig = classify("x", S)
    r_cs = classify("x", S, cs=True)
    assert "=" in r_cig.cigar
    assert r_cs.cigar.startswith(":")

def test_every_result_has_all_23_fields():
    from dataclasses import fields
    r = classify("x", _rnd(3000, 18))
    assert len(fields(r)) == 23

def test_identity_and_p_dist_are_consistent():
    ltr = _rnd(600, 19)
    S = ltr + _rnd(900, 20) + _evolve(ltr, 0.12, 21)
    r = classify("x", S)
    assert r.identity == pytest.approx(1.0 - r.p_dist, abs=1e-9)
