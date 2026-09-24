"""External credit moves an OUTER boundary only; the partner's INNER boundary follows
homology, and credited bases with no partner are gap columns, never mismatches.

`classify(..., spans=...)` measures a pair the caller already knows instead of
re-discovering one.
"""
import random

import pytest

from kmer2ltr.align import SIG_GAPS, Pairing, _extend, classify
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


def test_a_credited_flank_with_no_partner_leaves_the_inner_ends_alone():
    """Breaks if the credit drags the partner's inner end over the flank (forced pairing)."""
    ltr = _rnd(400, 100)
    S = _rnd(40, 101) + ltr + _rnd(1200, 102) + _evolve(ltr, 0.1, 103) + _rnd(40, 104)
    tight = classify("x", S)
    loose = classify("x", S, tsd_credit=1e6)
    assert tight.flank5_len > 0 and tight.flank3_len > 0
    assert (loose.ltr5_start, loose.ltr3_end) == (1, len(S))
    assert (loose.ltr5_end, loose.ltr3_start) == (tight.ltr5_end, tight.ltr3_start)
    assert (loose.k2p, loose.n_sites) == (tight.k2p, tight.n_sites)
    assert loose.n_gapcols == tight.n_gapcols + tight.flank5_len + tight.flank3_len
    assert loose.cigar.startswith(f"{tight.flank5_len}I")
    assert loose.cigar.endswith(f"{tight.flank3_len}D")
    assert loose.status == "pass"


def test_unpaired_flanks_at_both_ends_are_never_slid_into_mismatches():
    """Breaks if the unpaired runs go through end-to-end alignment: with a run at each
    end, sliding the two LTRs past each other is cheaper than two long gaps."""
    ltr = _rnd(400, 110)
    S = _rnd(400, 111) + ltr + _rnd(1200, 112) + _evolve(ltr, 0.05, 113) + _rnd(600, 114)
    tight = classify("x", S)
    loose = classify("x", S, tsd_credit=1e6)
    assert (loose.ltr5_start, loose.ltr3_end) == (1, len(S))
    assert loose.k2p == tight.k2p
    assert loose.n_sites == tight.n_sites


def test_a_credit_pairs_as_far_as_significant_homology_reaches():
    """Breaks if the partner's inner end ignores real homology beyond a deletion in the
    partner (under-pairing) or stops short of it."""
    a, b, c = _rnd(100, 120), _rnd(60, 121), _rnd(300, 122)
    a2, c2 = _evolve(a, 0.02, 123), _evolve(c, 0.02, 124)
    flank = _rnd(40, 125)
    S = flank + a + b + c + _rnd(1200, 126) + a2 + c2
    loose = classify("x", S, tsd_credit=1e6)
    a2_start = len(flank) + len(a) + len(b) + len(c) + 1200 + 1      # 1-based
    assert loose.ltr5_start == 1 and loose.ltr3_end == len(S)
    assert abs(loose.ltr3_start - a2_start) <= 3
    assert loose.cigar.startswith(f"{len(flank)}I")
    assert loose.k2p < 0.05                   # A, C pair; B and the flank are gaps


def test_homology_deep_in_the_partner_side_does_not_drag_the_inner_end():
    """Breaks if any significant hit in the partner window moves the inner end: a copy of
    the credited segment 500 bp inside the internal region is not its partner."""
    x = _rnd(400, 130)
    ltr = _rnd(400, 131)
    # 15% diverged: too far off-diagonal for the uncredited snap, still a strong local hit
    internal = _rnd(700, 132) + _evolve(x, 0.15, 133) + _rnd(500, 134)
    S = x + ltr + internal + _evolve(ltr, 0.05, 135)
    tight = classify("x", S)
    loose = classify("x", S, tsd_credit=1e6)
    assert tight.ltr5_start > 1
    assert loose.ltr5_start == 1
    assert loose.ltr3_start == tight.ltr3_start
    assert loose.k2p == tight.k2p


def test_a_pairing_beyond_an_insertion_next_to_the_core_keeps_the_insertion_as_gaps():
    """Breaks if a credited flank's pairing far out (beyond an insertion next to the
    called LTR) counts the insertion as paired, so that it is aligned end to end onto
    the partner, or if the real pairing beyond it is dropped."""
    p, ins, ltr = _rnd(300, 200), _rnd(1500, 201), _rnd(400, 202)
    S = p + ins + ltr + _rnd(1200, 203) + _evolve(p, 0.02, 204) + _evolve(ltr, 0.05, 205)
    l5b = len(p) + len(ins)
    l3b = l5b + 400 + 1200 + 300
    spans = (l5b, l5b + 399, l3b, l3b + 399)
    tight = classify("x", S, spans=spans)
    loose = classify("x", S, tsd_credit=1e6, spans=spans)
    assert tight.ltr5_start == l5b + 1                  # homology alone does not snap it
    assert loose.ltr5_start == 1
    assert abs(loose.ltr3_start - (l3b + 1 - len(p))) <= 3
    assert loose.n_gapcols >= len(ins)
    assert loose.k2p < 0.08


def test_a_credited_flank_pairs_only_as_far_as_homology_reaches_even_if_homology_carries_it():
    """Breaks if a credited flank whose inner part pairs well enough to carry the whole
    flank by homology alone gets every base paired: its outer part has no partner and
    would be aligned end to end, as a run of mismatches."""
    x, h, ltr = _rnd(400, 210), _rnd(1000, 211), _rnd(400, 212)
    S = x + h + ltr + _rnd(1200, 213) + _evolve(h, 0.05, 214) + _evolve(ltr, 0.05, 215)
    l5b = len(x) + len(h)
    l3b = l5b + 400 + 1200 + 1000
    spans = (l5b, l5b + 399, l3b, l3b + 399)
    tight = classify("x", S, spans=spans)
    loose = classify("x", S, tsd_credit=1e6, spans=spans)
    assert tight.ltr5_start == 1                        # homology alone carries all of it
    assert loose.ltr5_start == 1
    assert abs(loose.ltr3_start - (l3b + 1 - len(h))) <= 3
    assert loose.cigar.startswith(f"{len(x)}I")
    assert loose.k2p < 0.08 < tight.k2p


def test_the_stretch_between_the_core_and_a_credited_flanks_pairing_is_gap_columns():
    """Breaks if a block that sits between the core and the paired part of a credited
    flank on BOTH sides (the credited LTR's 250 bp and the partner's 150 bp, unrelated)
    is aligned end to end: two facing blocks are cheaper as mismatches than as gaps."""
    h, bl, br, ltr = _rnd(300, 220), _rnd(250, 221), _rnd(150, 222), _rnd(400, 223)
    S = h + bl + ltr + _rnd(1200, 224) + _evolve(h, 0.02, 225) + br + _evolve(ltr, 0.02, 226)
    l5b = len(h) + len(bl)
    l3b = l5b + 400 + 1200 + len(h) + len(br)
    spans = (l5b, l5b + 399, l3b, l3b + 399)
    loose = classify("x", S, tsd_credit=1e6, spans=spans)
    assert loose.ltr5_start == 1
    assert abs(loose.ltr3_start - (l3b + 1 - len(br) - len(h))) <= 3
    assert loose.k2p < 0.04
    assert loose.n_gapcols >= len(bl) + len(br)


def test_a_long_credited_flank_without_a_collinear_pairing_stays_unpaired():
    """Breaks if a credited flank that homology alone would carry is paired with a
    homolog that starts more than MAX_PARTNER_SKIP partner bases from the core: the
    inner boundary must not be dragged there on the strength of that hit."""
    h, br, ltr = _rnd(600, 230), _rnd(300, 231), _rnd(400, 232)
    S = h + ltr + _rnd(1200, 233) + _evolve(h, 0.02, 234) + br + _evolve(ltr, 0.02, 235)
    l5b = len(h)
    l3b = l5b + 400 + 1200 + len(h) + len(br)
    spans = (l5b, l5b + 399, l3b, l3b + 399)
    tight = classify("x", S, spans=spans)
    loose = classify("x", S, tsd_credit=1e6, spans=spans)
    assert tight.ltr5_start == 1                        # homology alone carries it
    assert loose.ltr5_start == 1
    assert loose.ltr3_start == l3b + 1
    assert loose.cigar.startswith(f"{len(h)}I")


@pytest.mark.parametrize("flank, div, seed", [
    (50, 0.05, 2), (50, 0.10, 1), (60, 0.15, 0),       # too short to show significance
    (60, 0.15, 6), (80, 0.10, 5),                       # a few mismatched end bases
])
def test_a_homologous_credited_flank_pairs_as_it_does_without_the_credit(flank, div, seed):
    """Breaks if a credit makes a flank that homology alone carries pair worse than no
    credit does: a 50-60 bp homolog can miss the local significance test, and a local
    alignment drops a few mismatched end bases -- neither is evidence the flank lacks
    a partner, and marking it unpaired moves the partner's inner end and lowers K2P."""
    ltr = _rnd(500, 1000 + seed)
    S = ltr + _rnd(1500, 3000 + seed) + _evolve(ltr, div, 2000 + seed)
    spans = (flank, 499, 2000 + flank, len(S) - 1)
    tight = classify("x", S, spans=spans)
    loose = classify("x", S, tsd_credit=1e6, spans=spans)
    assert tight.ltr5_start == loose.ltr5_start == 1
    assert (loose.ltr3_start, loose.cigar, loose.k2p) == (tight.ltr3_start, tight.cigar,
                                                          tight.k2p)


def test_a_long_credited_unpaired_flank_does_not_make_the_pair_weak():
    """Breaks if significance is scored over credited bases that have no partner."""
    ltr = _rnd(400, 140)
    S = _rnd(1500, 141) + ltr + _rnd(1200, 142) + _evolve(ltr, 0.05, 143)
    loose = classify("x", S, tsd_credit=1e6)
    assert loose.ltr5_start == 1
    assert loose.status == "pass"


def test_extend_reports_the_pairing_only_when_the_credit_decided():
    """Breaks if `_extend` stops telling its caller how a credited flank pairs."""
    core_side = _rnd(300, 150)
    homologous = _evolve(core_side[:100], 0.01, 151)
    junk = _rnd(100, 152)
    k_out, k_in, _, unpaired = _extend(homologous, core_side, GENERIC_MATRIX, SIG_GAPS,
                                       10.0, False)
    assert (k_out, unpaired) == (100, None) and 95 <= k_in <= 105
    assert _extend(junk, core_side, GENERIC_MATRIX, SIG_GAPS, 10.0, False)[::3] == (0, None)
    k_out, k_in, margin, pairing = _extend(junk, core_side, GENERIC_MATRIX, SIG_GAPS,
                                           10.0, False, credit=1e6)
    assert (k_out, k_in, margin, pairing) == (100, 0, None, Pairing(100, 0, 0, 0, 0))


def test_the_graded_scan_ignores_the_credit():
    """Breaks if a credit too small to carry a flank still lowers the graded noise floor."""
    ltr = _rnd(400, 160)
    S = _rnd(300, 161) + ltr + _rnd(1200, 162) + _evolve(ltr, 0.1, 163) + _rnd(300, 164)
    free = classify("x", S, snap_mode="graded")
    small = classify("x", S, snap_mode="graded", tsd_credit=2.0)
    assert (small.ltr5_start, small.ltr5_end, small.ltr3_start, small.ltr3_end) == \
        (free.ltr5_start, free.ltr5_end, free.ltr3_start, free.ltr3_end)


def test_known_spans_reproduce_the_call():
    """Breaks if measuring a known pair differs from measuring the pair discovery found."""
    ltr = _rnd(400, 170)
    S = _rnd(50, 171) + ltr + _rnd(900, 172) + _evolve(ltr, 0.05, 173) + _rnd(60, 174)
    ref = classify("x", S)
    spans = (ref.ltr5_start - 1, ref.ltr5_end - 1, ref.ltr3_start - 1, ref.ltr3_end - 1)
    again = classify("x", S, spans=spans)
    for f in ("status", "ltr5_start", "ltr5_end", "ltr3_start", "ltr3_end", "k2p",
              "n_sites", "cigar", "motif"):
        assert getattr(again, f) == getattr(ref, f), f


def test_known_spans_are_measured_instead_of_the_pair_discovery_prefers():
    """Breaks if `spans=` is ignored, or if Stage 4 swaps the given pair for an outer one:
    discovery (with Stage 4) reports the outer pair; the caller asks for the nested one."""
    outer, inner = _rnd(300, 180), _rnd(300, 181)
    nested = inner + _rnd(800, 182) + inner
    pre = outer + _rnd(400, 183)
    S = pre + nested + _rnd(400, 184) + _evolve(outer, 0.08, 185)
    found = classify("x", S)
    assert found.ltr5_start == 1                              # discovery took the outer pair
    b = len(pre)
    spans = (b, b + 299, b + 300 + 800, b + 300 + 800 + 299)
    given = classify("x", S, spans=spans)
    assert (given.ltr5_start, given.ltr5_end) == (b + 1, b + 300)
    assert (given.ltr3_start, given.ltr3_end) == (b + 1101, b + 1400)
    assert given.k2p == 0.0


def test_invalid_spans_are_reported_not_measured():
    """Breaks if an impossible pair is silently measured."""
    S = _rnd(2000, 190)
    assert classify("x", S, spans=(10, 5, 20, 30)).status == "bad_spans"
    assert classify("x", S, spans=(0, 99, 50, 199)).status == "bad_spans"     # overlapping
    assert classify("x", S, spans=(0, 99, 1900, 2000)).status == "bad_spans"  # past the end
