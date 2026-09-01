"""Tests for the Stage 3 rewrite (window-free extension regions, graded mode,
joint inner boundaries) and the Stage 4 re-run of Stages 2-3."""
import random

import pytest

from kmer2ltr.align import (SIG_GAPS, Bounds, calibrate, classify, discover,
                          ltr_spans, outermost, snap_bounds, terminal_snap)
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
            out.append(ti[c] if r.random() < 2 / 3 else r.choice(tv[c]))
        else:
            out.append(c)
    return "".join(out)


def _hit_and_matrix(S):
    h = discover(S, GENERIC_MATRIX, SIG_GAPS)
    m, _, _ = calibrate(S, h)
    return (discover(S, m, SIG_GAPS, w0=h.w) or h), m


# --------------------------------------------------------------------------- #
# Stage 3 is now window-free
# --------------------------------------------------------------------------- #

def test_snap_bounds_needs_no_window_and_agrees_with_the_hit_form():
    """snap_bounds takes spans, not a Hit -- that is what makes it re-runnable
    after Stage 4, where no discovery window exists."""
    ltr = _rnd(400, 1)
    S = _rnd(40, 2) + ltr + _rnd(1200, 3) + ltr + _rnd(40, 4)
    hit, m = _hit_and_matrix(S)
    a = terminal_snap(S, hit, m)
    b = snap_bounds(S, ltr_spans(S, hit), m, SIG_GAPS, 10.0)
    assert a.spans == b.spans


def test_extension_partner_survives_a_hit_flush_against_the_window_edge():
    """The window-bounded form used S[wstart:l3b] as the 5' partner region,
    which is EMPTY when the hit begins exactly at the suffix window's start --
    the snap could then never fire. Feeding spans that reproduce that geometry
    must still produce a decision."""
    ltr = _rnd(300, 5)
    S = _rnd(30, 6) + ltr + _rnd(900, 7) + ltr
    hit, m = _hit_and_matrix(S)
    l5b, l5e, l3b, l3e = ltr_spans(S, hit)
    # spans identical to the real hit, but scored with no window at all
    b = snap_bounds(S, (l5b, l5e, l3b, l3e), m, SIG_GAPS, 10.0)
    assert 0 <= b.l5b <= b.l5e < b.l3b <= b.l3e <= len(S) - 1
    assert b.margin_bits is not None


def test_bounds_stay_ordered_in_every_stage3_mode():
    ltr = _rnd(300, 8)
    S = _rnd(60, 9) + ltr + _rnd(800, 10) + _evolve(ltr, 0.2, 11) + _rnd(60, 12)
    hit, m = _hit_and_matrix(S)
    for mode in ("binary", "graded"):
        for inner in ("none", "joint"):
            b = snap_bounds(S, ltr_spans(S, hit), m, SIG_GAPS, 10.0,
                            mode=mode, inner=inner)
            assert 0 <= b.l5b <= b.l5e < b.l3b <= b.l3e <= len(S) - 1, (mode, inner)


def test_unknown_stage3_modes_fail_loudly():
    ltr = _rnd(300, 13)
    S = ltr + _rnd(800, 14) + ltr
    hit, m = _hit_and_matrix(S)
    with pytest.raises(ValueError):
        snap_bounds(S, ltr_spans(S, hit), m, SIG_GAPS, 10.0, mode="nope")
    with pytest.raises(ValueError):
        snap_bounds(S, ltr_spans(S, hit), m, SIG_GAPS, 10.0, inner="nope")


# --------------------------------------------------------------------------- #
# Graded mode
# --------------------------------------------------------------------------- #

def test_graded_mode_still_snaps_a_perfectly_bounded_element():
    """Grading changes where a flank STARTS, never whether a clean element is
    left alone."""
    for seed in range(8):
        anc = _rnd(400, seed)
        S = _evolve(anc, 0.06, seed + 50) + _rnd(1000, seed + 90) + _evolve(anc, 0.06, seed + 130)
        hit, m = _hit_and_matrix(S)
        b = snap_bounds(S, ltr_spans(S, hit), m, SIG_GAPS, 10.0, mode="graded")
        assert b.l5b == 0 and b.l3e == len(S) - 1


def test_graded_mode_never_calls_a_longer_flank_than_binary():
    """Grading can only move a boundary OUTWARD from the core, so the flank it
    reports is at most the one the binary rule reports."""
    for seed in range(10):
        ltr = _rnd(400, seed + 200)
        S = (_rnd(80, seed + 210) + ltr + _rnd(1000, seed + 220)
             + _evolve(ltr, 0.25, seed + 230) + _rnd(80, seed + 240))
        hit, m = _hit_and_matrix(S)
        spans = ltr_spans(S, hit)
        bb = snap_bounds(S, spans, m, SIG_GAPS, 10.0, mode="binary")
        bg = snap_bounds(S, spans, m, SIG_GAPS, 10.0, mode="graded")
        assert bg.l5b <= bb.l5b
        assert (len(S) - 1 - bg.l3e) <= (len(S) - 1 - bb.l3e)


def test_graded_and_binary_agree_whenever_the_snap_fires():
    """When the whole flank candidate is homologous both rules reach the
    terminus, so they can only differ on elements that were going to be
    flank-called anyway."""
    ltr = _rnd(400, 300)
    S = ltr + _rnd(1200, 301) + _evolve(ltr, 0.05, 302)
    hit, m = _hit_and_matrix(S)
    spans = ltr_spans(S, hit)
    assert (snap_bounds(S, spans, m, SIG_GAPS, 10.0, mode="binary").spans
            == snap_bounds(S, spans, m, SIG_GAPS, 10.0, mode="graded").spans)


# --------------------------------------------------------------------------- #
# Joint inner boundaries
# --------------------------------------------------------------------------- #

def test_joint_inner_leaves_the_outer_boundaries_alone():
    """The inner refinement is anchored on the outer calls; it must never move
    them, or Stage 3's whole decision would be silently re-opened."""
    ltr = _rnd(500, 400)
    S = _rnd(70, 401) + ltr + _rnd(1200, 402) + _evolve(ltr, 0.15, 403)
    hit, m = _hit_and_matrix(S)
    spans = ltr_spans(S, hit)
    base = snap_bounds(S, spans, m, SIG_GAPS, 10.0, inner="none")
    joint = snap_bounds(S, spans, m, SIG_GAPS, 10.0, inner="joint")
    assert (joint.l5b, joint.l3e) == (base.l5b, base.l3e)


def test_joint_inner_recovers_a_clean_ltr_length():
    ltr = _rnd(500, 410)
    S = ltr + _rnd(1500, 411) + ltr
    r = classify("x", S, inner="joint")
    assert r.status == "pass"
    assert r.ltr5_len == pytest.approx(500, abs=5)
    assert r.ltr3_len == pytest.approx(500, abs=5)


# --------------------------------------------------------------------------- #
# Stage 4 re-runs Stages 2-3
# --------------------------------------------------------------------------- #

def test_stage4_rerun_snaps_the_recovered_outer_pair_to_the_termini():
    """The spec says the outer pair wins and Stages 2-3 are re-run on it. The
    shipped code stopped at raw Smith-Waterman ends, so the recovered outer pair
    was the one part of the element whose termini were never tested -- and raw
    SW trims a terminus precisely when the outer pair is diverged."""
    reran = plain = n = 0
    for seed in range(12):
        outer = _rnd(500, 500 + seed)
        inner = _rnd(400, 600 + seed)
        nested = inner + _rnd(600, 700 + seed) + _evolve(inner, 0.01, 800 + seed)
        S = (_evolve(outer, 0.15, 900 + seed) + _rnd(400, 1000 + seed) + nested
             + _rnd(400, 1100 + seed) + _evolve(outer, 0.15, 1200 + seed))
        a = classify("x", S, stage4_recal=True)
        b = classify("x", S, stage4_recal=False)
        if a.ltr5_start is None or b.ltr5_start is None:
            continue
        n += 1
        # both must find the OUTER pair rather than the younger nested one
        assert a.ltr5_start < 120 and b.ltr5_start < 120, seed
        reran += (a.ltr5_start == 1 and a.ltr3_end == len(S))
        plain += (b.ltr5_start == 1 and b.ltr3_end == len(S))
    assert n >= 10
    assert reran > plain, (
        f"re-running Stages 2-3 after Stage 4 reached both termini on {reran}/{n}, "
        f"raw Stage 4 ends on {plain}/{n} -- the re-run bought nothing")


def test_stage4_rerun_does_not_disturb_an_element_stage4_never_fires_on():
    ltr = _rnd(400, 510)
    S = ltr + _rnd(1000, 511) + ltr
    assert (classify("x", S, stage4_recal=True).ltr5_start
            == classify("x", S, stage4_recal=False).ltr5_start == 1)


def test_stage4_still_refuses_to_invent_a_pair_from_unrelated_flanks():
    """The re-run must not weaken the acceptance gate: it recalibrates only
    AFTER a pair has already cleared the generic-matrix significance test."""
    spurious = 0
    trials = 0
    for div in (0.01, 0.05):
        for flank in (300, 600):
            for seed in range(15):
                trials += 1
                anc = _rnd(400, seed * 7)
                S = (_rnd(flank, seed + 1) + _evolve(anc, div / 2, seed + 2)
                     + _rnd(1200, seed + 3) + _evolve(anc, div / 2, seed + 4)
                     + _rnd(flank, seed + 5))
                hit, m = _hit_and_matrix(S)
                b0 = snap_bounds(S, ltr_spans(S, hit), m, SIG_GAPS, 10.0)
                if outermost(S, b0, m).spans != b0.spans:
                    spurious += 1
    assert spurious == 0, f"invented an outer pair on {spurious}/{trials}"


# --------------------------------------------------------------------------- #
# The parasail semantics the inner refinement depends on
# --------------------------------------------------------------------------- #

def test_sg_qe_db_anchors_the_outer_ends_and_frees_the_inner_ones():
    """`_joint_inner` is only correct if sg_qe_db means "query begin and ref end
    anchored, query end and ref begin free", and if a reversed sg_de pass
    recovers the ref begin. Both are verified here rather than assumed, because
    parasail's flag names describe free GAPS, not free ends, and the two read
    the same until they disagree."""
    import parasail
    from kmer2ltr.scoring import GENERIC_MATRIX, SCALE
    go, ge = 6 * SCALE, 2 * SCALE
    core = "ACGTACGTACGTACGTACGT"          # 20 bp shared repeat
    q = core + "T" * 15                     # repeat at the query's START
    r = "G" * 15 + core                     # repeat at the ref's END
    res = parasail.sg_qe_db_striped_sat(q, r, go, ge, GENERIC_MATRIX)
    assert res.score == 20 * SCALE          # all 20 matches, nothing else paid for
    assert res.end_query == 19               # query end free: stops after the repeat
    assert res.end_ref == len(r) - 1         # ref end anchored
    rev = parasail.sg_de_striped_sat(q[:res.end_query + 1][::-1], r[::-1],
                                     go, ge, GENERIC_MATRIX)
    assert len(r) - 1 - rev.end_ref == 15    # ref begin recovered exactly


# --------------------------------------------------------------------------- #
# The parasail table lifetime footgun
# --------------------------------------------------------------------------- #

def test_score_table_read_matches_the_reported_score():
    """`score_table` is a numpy VIEW onto memory the parasail result owns, and it
    does not keep that result alive. Reading it off a temporary
    (`np.asarray(f(...).score_table)`) is a use-after-free: silent garbage while
    the freed pages happen to stay mapped, a segfault once the table is large
    enough for the allocator to return them.

    For `sg_de` the query is fully consumed, so the maximum of the table's LAST
    ROW is exactly the reported score. That identity holds for a correctly-owned
    table and fails for a dangling one, which is what makes it a usable guard.
    """
    import numpy as np
    import parasail
    from kmer2ltr.scoring import GENERIC_MATRIX, SCALE
    go, ge = 6 * SCALE, 2 * SCALE
    q = _rnd(600, 900)
    r = q[:300] + _rnd(1200, 901)
    res = parasail.sg_de_table_striped_sat(q, r, go, ge, GENERIC_MATRIX)
    tab = np.asarray(res.score_table)
    assert tab.shape == (len(q), len(r))
    assert int(tab[-1].max()) == res.score


def test_graded_path_handles_a_large_table_without_crashing_or_garbage():
    """The record that first exposed the lifetime bug allocated a 2000 x 4200
    table. Garbage in that table makes the chosen endpoint effectively random, so
    a sane flank call on a big, heavily diverged element is the observable
    signature of a correctly-owned one."""
    ltr = _rnd(700, 910)
    S = (_rnd(500, 911) + _evolve(ltr, 0.25, 912) + _rnd(9000, 913)
         + _evolve(ltr, 0.25, 914) + _rnd(500, 915))
    r = classify("x", S, snap_mode="graded", t_bits=10.0)
    assert r.status in ("pass", "weak_pair")
    assert 0 <= r.flank5_len <= 500 and 0 <= r.flank3_len <= 500
    assert r.ltr5_start <= r.ltr5_end < r.ltr3_start <= r.ltr3_end == r.seq_len - r.flank3_len
