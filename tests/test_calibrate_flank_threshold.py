"""The schedule rule decides a shipped default, so its logic is tested on
synthetic stats where the right answer is known by construction."""
import math

import pytest

from bench.calibrate_flank_threshold import (DHAT_LABELS, LARGE_TOL, MID_TOL, MIN_BIN_N,
                                   REFERENCE_T, as_python, derive)

TS = (2.0, 5.0, 10.0, 20.0)


def _stats(ff, det, lost=0.1, n_flank=10_000, large=(50, 100), mid=20):
    row = {"ff": ff, "lost": lost, "n_flank": n_flank}
    for f in large:
        row[f"det{f}"] = det
    row[f"det{mid}"] = det
    return row


def _grid(per_t, large=(50, 100), mid=20, n_flank=10_000):
    """per_t: {t: {bin: (ff, det)}} -> the shape derive() consumes."""
    return {t: {b: _stats(ff, det, n_flank=n_flank, large=large, mid=mid)
                for b, (ff, det) in bins.items()}
            for t, bins in per_t.items()}


def _flat(ff_by_t, det_by_t=None):
    det_by_t = det_by_t or {t: 0.9 for t in ff_by_t}
    return {t: {b: (ff_by_t[t], det_by_t[t]) for b in DHAT_LABELS} for t in ff_by_t}


def test_an_easy_bin_drops_below_the_reference_threshold():
    """Rule step 3 is smallest-that-suffices: in a bin whose false-flank rate is
    already under the pooled target even at a permissive threshold, `t` must come
    DOWN from the reference, because every increment above what the target needs
    is detection given away for nothing."""
    per_t = {2.0: {}, 5.0: {}, 10.0: {}, 20.0: {}}
    for b in DHAT_LABELS:
        easy = b <= 0.05
        per_t[2.0][b] = (0.001 if easy else 0.40, 0.9)
        per_t[5.0][b] = (0.001 if easy else 0.30, 0.9)
        per_t[10.0][b] = (0.001 if easy else 0.20, 0.9)
        per_t[20.0][b] = (0.001 if easy else 0.05, 0.9)
    g = _grid(per_t)
    h = _grid(per_t, large=(45, 65), mid=25)
    schedule, ff_target, _ = derive(g, h, verbose=False)
    got = dict(schedule)
    assert 0.001 < ff_target < 0.20
    assert got[0.0] == 2.0 and got[0.05] == 2.0      # easy bins relax
    assert got[0.5] > 2.0                             # hard bins do not


def test_a_bin_above_target_climbs_until_it_meets_it():
    per_t = {2.0: {}, 5.0: {}, 10.0: {}, 20.0: {}}
    for b in DHAT_LABELS:
        hard = b >= 0.3
        per_t[2.0][b] = (0.50 if hard else 0.01, 0.9)
        per_t[5.0][b] = (0.30 if hard else 0.008, 0.9)
        per_t[10.0][b] = (0.20 if hard else 0.005, 0.9)
        per_t[20.0][b] = (0.05 if hard else 0.002, 0.9)
    g = _grid(per_t)
    h = _grid(per_t, large=(45, 65), mid=25)
    schedule, ff_target, _ = derive(g, h, verbose=False)
    got = dict(schedule)
    assert got[0.0] == 2.0                      # easy bins sit at the bottom
    assert got[0.4] == 20.0                     # hard bins climb
    assert got[0.5] == 20.0


def test_the_detection_floor_blocks_a_t_that_would_meet_the_target():
    """The whole point of the floor: a threshold that reaches the false-flank
    target by abandoning large-flank detection is not admissible."""
    per_t = {}
    for t, ff, det in ((2.0, 0.50, 0.90), (5.0, 0.40, 0.90),
                       (10.0, 0.30, 0.90), (20.0, 0.01, 0.50)):
        per_t[t] = {b: (ff, det) for b in DHAT_LABELS}
    g = _grid(per_t)
    h = _grid(per_t, large=(45, 65), mid=25)
    schedule, _, rows = derive(g, h, verbose=False)
    # t=20 reaches the target but drops detection 40 points below the reference,
    # so it is inadmissible and the schedule cannot use it anywhere
    assert all(t == 10.0 for _, t in schedule)
    assert all(20.0 not in r["admissible"] for r in rows if r["t"] is not None)


def test_a_detection_drop_inside_tolerance_is_still_admissible():
    """The floor is a tolerance, not equality: a threshold that gives up slightly
    less than LARGE_TOL of detection stays available."""
    per_t = {}
    for t, ff, det in ((2.0, 0.50, 0.90), (5.0, 0.40, 0.90),
                       (10.0, 0.30, 0.90), (20.0, 0.01, 0.90 - LARGE_TOL + 0.001)):
        per_t[t] = {b: (ff, det) for b in DHAT_LABELS}
    _, _, rows = derive(_grid(per_t), _grid(per_t, large=(45, 65), mid=25),
                        verbose=False)
    assert all(20.0 in r["admissible"] for r in rows if r["t"] is not None)


def test_the_schedule_is_forced_monotone():
    """Rule step 4: a more diverged element may never need LESS evidence."""
    per_t = {2.0: {}, 5.0: {}, 10.0: {}, 20.0: {}}
    for b in DHAT_LABELS:
        # deliberately anti-monotone: the MIDDLE bin is the only hard one
        hard = b == 0.2
        per_t[2.0][b] = (0.50 if hard else 0.001, 0.9)
        per_t[5.0][b] = (0.001, 0.9)
        per_t[10.0][b] = (0.001, 0.9)
        per_t[20.0][b] = (0.001, 0.9)
    schedule, _, _ = derive(_grid(per_t), _grid(per_t, large=(45, 65), mid=25),
                            verbose=False)
    values = [t for _, t in schedule]
    assert values == sorted(values), values
    assert values[-1] >= values[DHAT_LABELS.index(0.2)]


def test_an_undersized_gold_bin_is_not_estimated():
    per_t = _flat({2.0: 0.01, 5.0: 0.008, 10.0: 0.005, 20.0: 0.002})
    g = _grid(per_t)
    for t in g:
        g[t][0.5]["n_flank"] = MIN_BIN_N - 1
    h = _grid(per_t, large=(45, 65), mid=25)
    _, _, rows = derive(g, h, verbose=False)
    assert [r for r in rows if r["bin"] == 0.5][0]["why"] == "gold bin too small"


def test_an_undersized_homology_bin_abstains_instead_of_vetoing():
    """The homology grid stops at 35% mutated, so its top bins are empty by
    construction; an empty bin must not veto the gold grid's evidence."""
    per_t = {2.0: {}, 5.0: {}, 10.0: {}, 20.0: {}}
    for b in DHAT_LABELS:
        hard = b == 0.5                            # only the top bin misses the target
        per_t[2.0][b] = (0.60 if hard else 0.01, 0.9)
        per_t[5.0][b] = (0.50 if hard else 0.01, 0.9)
        per_t[10.0][b] = (0.40 if hard else 0.01, 0.9)
        per_t[20.0][b] = (0.02 if hard else 0.01, 0.9)
    g = _grid(per_t)
    h = _grid(per_t, large=(45, 65), mid=25)
    for t in h:                                   # homology says t=20 is terrible...
        h[t][0.5]["det45"] = 0.9 if t != 20.0 else 0.1
        h[t][0.5]["n_flank"] = MIN_BIN_N - 1      # ...but has no evidence there
    schedule, _, _ = derive(g, h, verbose=False)
    assert dict(schedule)[0.5] == 20.0


def test_missing_reference_is_a_hard_error():
    per_t = _flat({2.0: 0.01, 5.0: 0.008})
    with pytest.raises(SystemExit):
        derive(_grid(per_t), _grid(per_t, large=(45, 65), mid=25), verbose=False)


def test_as_python_merges_equal_bins_and_ends_at_infinity():
    lit = as_python([(0.0, 5.0), (0.05, 5.0), (0.1, 10.0), (0.2, 10.0),
                     (0.3, 15.0), (0.4, 15.0), (0.5, 20.0)])
    # 4 merged entries + the outer tuple + float("inf")'s own call
    assert lit.count("(") == 6
    assert 'float("inf")' in lit
    ns = {}
    exec(lit, ns)
    sched = ns["T_BITS_SCHEDULE"]
    assert [t for _, t in sched] == [5, 10, 15, 20]
    assert sched[-1][0] == math.inf
    edges = [e for e, _ in sched]
    assert edges == sorted(edges)
