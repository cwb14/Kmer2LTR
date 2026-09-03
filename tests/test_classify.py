import random
from dataclasses import astuple, fields

import pytest
from kmer2ltr.align import Result, classify
from kmer2ltr.fasta import sanitize

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

def test_every_result_has_all_29_fields():
    r = classify("x", _rnd(3000, 18))
    assert len(fields(r)) == 29

def test_identity_and_p_dist_are_consistent():
    ltr = _rnd(600, 19)
    S = ltr + _rnd(900, 20) + _evolve(ltr, 0.12, 21)
    r = classify("x", S)
    assert r.identity == pytest.approx(1.0 - r.p_dist, abs=1e-9)


# --------------------------------------------------------------------------- #
# Terminal motif and insertion time
# --------------------------------------------------------------------------- #

def test_motif_reports_the_dinucleotides_at_the_called_boundaries():
    ltr = "TG" + _rnd(296, 30) + "CA"
    S = ltr + _rnd(1000, 31) + ltr
    r = classify("x", S)
    assert r.status == "pass" and (r.ltr5_start, r.ltr3_end) == (1, len(S))
    assert r.motif == "tg...ca"


def test_motif_is_read_off_the_boundary_and_never_searched_for():
    """The tool uses no terminal-motif prior anywhere, which is what makes the
    `TG`..`CA` rate on real data an INDEPENDENT accuracy check rather than a
    restatement of its own call. An element with non-canonical termini must be
    bounded just as well, and simply report what is there."""
    ltr = "AA" + _rnd(296, 32) + "GG"
    S = ltr + _rnd(1000, 33) + ltr
    r = classify("x", S)
    assert (r.ltr5_start, r.ltr3_end) == (1, len(S))
    assert r.motif == "aa...gg"


def test_motif_follows_a_flank_correction():
    """The motif must come from the CORRECTED boundary, not the raw input ends
    -- otherwise it could not measure whether the correction was right."""
    ltr = "TG" + _rnd(396, 34) + "CA"
    S = _rnd(120, 35) + ltr + _rnd(900, 36) + ltr + _rnd(90, 37)
    r = classify("x", S)
    assert r.flank5_len > 0 and r.flank3_len > 0
    assert r.motif == "tg...ca", f"motif {r.motif} with flanks {r.flank5_len}/{r.flank3_len}"


def test_k2p_time_needs_a_mutation_rate_and_matches_the_distance():
    ltr = _rnd(400, 38)
    S = ltr + _rnd(1000, 39) + _evolve(ltr, 0.08, 40)
    assert classify("x", S).k2p_time is None
    r = classify("x", S, mutation_rate=7e-9)
    assert r.k2p_time == round(r.k2p / (2 * 7e-9))


def test_a_mutation_rate_changes_nothing_else_in_the_row():
    """`-u` is a unit conversion applied after the measurement; it must not be
    able to move a boundary or a distance."""
    from dataclasses import replace
    ltr = _rnd(350, 41)
    S = _rnd(60, 42) + ltr + _rnd(900, 43) + _evolve(ltr, 0.15, 44)
    plain = classify("x", S)
    dated = classify("x", S, mutation_rate=1.3e-8)
    assert dated.k2p_time is not None
    assert replace(dated, k2p_time=None) == plain


def test_rows_without_a_pair_report_no_motif_and_no_time():
    for r in (classify("x", "ACGT" * 10, mutation_rate=7e-9),
              classify("x", "N" * 500, mutation_rate=7e-9),
              classify("x", _rnd(3000, 45), mutation_rate=7e-9)):
        assert r.motif is None and r.k2p_time is None, r.status


# --------------------------------------------------------------------------- #
# Genome columns and the external-evidence credit
# --------------------------------------------------------------------------- #

def test_result_appends_the_genome_columns_after_k2p_time():
    """They are appended, never inserted, so every documented `cut -f` recipe
    still selects the same fields."""
    names = [f.name for f in fields(Result)]
    assert names[-5:] == ["k2p_time", "orientation", "tsd", "tsd_offset", "tsd_input"]
    assert names.index("cigar") == 22
    assert names.index("motif") == 23
    assert names.index("k2p_time") == 24


def test_the_genome_columns_start_empty():
    ltr = _rnd(400, 21)
    r = classify("x", ltr + _rnd(1200, 22) + ltr)
    assert (r.orientation, r.tsd, r.tsd_offset, r.tsd_input) == (None, None, None, None)


def test_zero_credit_reproduces_the_shipped_call_exactly():
    ltr = _rnd(400, 23)
    S = _rnd(40, 24) + ltr + _rnd(1200, 25) + _evolve(ltr, 0.1, 26) + _rnd(17, 27)
    assert astuple(classify("x", S)) == astuple(classify("x", S, tsd_credit=0.0))


def test_a_credit_pushes_a_boundary_out_to_the_terminus_and_nowhere_else():
    """The binary snap returns the whole candidate flank or none of it, so
    external evidence can move a boundary to the sequence end and to no other
    position -- it cannot invent an interior boundary."""
    ltr = _rnd(400, 28)
    S = _rnd(40, 29) + ltr + _rnd(1200, 30) + _evolve(ltr, 0.1, 31) + _rnd(40, 32)
    tight = classify("x", S)
    loose = classify("x", S, tsd_credit=1e6)
    assert tight.flank5_len > 0 and tight.flank3_len > 0
    assert (loose.flank5_len, loose.flank3_len) == (0, 0)
    assert loose.ltr5_start == 1 and loose.ltr3_end == loose.seq_len


def test_a_credit_leaves_an_already_unflanked_element_alone():
    ltr = _rnd(400, 33)
    S = ltr + _rnd(1200, 34) + _evolve(ltr, 0.1, 35)
    assert astuple(classify("x", S)) == astuple(classify("x", S, tsd_credit=1e6))
