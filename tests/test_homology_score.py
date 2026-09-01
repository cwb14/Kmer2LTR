"""`score_homology` decides every verdict in the campaign, so its accounting is
tested against hand-built rows rather than trusted."""
import csv

import pytest

from bench.homology_grid import (TRUTH_COLUMNS, cells_from_json, cells_to_json,
                                 dhat_bin, score_homology)
from kmer2ltr.runner import COLUMNS


def _truth_row(**kw):
    row = dict.fromkeys(TRUTH_COLUMNS, "0")
    row.update(seq_id="x", source="lib", panel="subs", p_target="0.1",
               indel_rate="0", flank5="0", flank3="0", ltr5_start="1",
               ltr5_end="100", ltr3_start="201", ltr3_end="300", seq_len="300",
               realized_p="0.1", realized_k2p="0.11")
    row.update(kw)
    return row


def _pred_row(**kw):
    row = dict.fromkeys(COLUMNS, "NA")
    row.update(seq_id="x", seq_len="300", status="pass", ltr5_start="1",
               ltr5_end="100", ltr3_start="201", ltr3_end="300",
               flank5_len="0", flank3_len="0", k2p="0.11")
    row.update(kw)
    return row


def _write(tmp_path, truths, preds):
    t = tmp_path / "truth.tsv"
    p = tmp_path / "pred.tsv"
    with open(t, "w", newline="") as fh:
        w = csv.DictWriter(fh, TRUTH_COLUMNS, delimiter="\t")
        w.writeheader()
        for r in truths:
            w.writerow(r)
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, COLUMNS, delimiter="\t")
        w.writeheader()
        for r in preds:
            w.writerow(r)
    return t, p


KEY = ("lib", "subs", 0.1, 0.0, 0)


def test_a_perfect_call_scores_perfectly(tmp_path):
    t, p = _write(tmp_path, [_truth_row()], [_pred_row()])
    c = score_homology(t, p)["by_p"][KEY]
    assert c["n"] == 1 and c["n_pass"] == 1 and c["n_located"] == 1
    assert c.get("n_lost", 0) == 0 and c.get("n_false_flank", 0) == 0
    for f in ("ltr5_start", "ltr5_end", "ltr3_start", "ltr3_end"):
        assert c[f"abs_{f}"] == 0 and c[f"exact_{f}"] == 1
    assert c["k2p_err"] == pytest.approx(0.0)


def test_boundary_error_is_absolute_and_per_coordinate(tmp_path):
    t, p = _write(tmp_path, [_truth_row()],
                  [_pred_row(ltr5_end="97", ltr3_start="205")])
    c = score_homology(t, p)["by_p"][KEY]
    assert c["abs_ltr5_end"] == 3 and c["exact_ltr5_end"] == 0
    assert c["abs_ltr3_start"] == 4 and c["exact_ltr3_start"] == 0
    assert c["abs_ltr5_start"] == 0 and c["exact_ltr5_start"] == 1


def test_a_false_flank_is_counted_only_where_no_flank_exists(tmp_path):
    t, p = _write(tmp_path, [_truth_row()],
                  [_pred_row(ltr5_start="6", flank5_len="5")])
    assert score_homology(t, p)["by_p"][KEY]["n_false_flank"] == 1
    t2, p2 = _write(tmp_path, [_truth_row(flank5="5", flank3="5", ltr5_start="6",
                                          ltr5_end="105", ltr3_start="206",
                                          ltr3_end="305", seq_len="310")],
                    [_pred_row(seq_len="310", ltr5_start="6", ltr5_end="105",
                               ltr3_start="206", ltr3_end="305",
                               flank5_len="5", flank3_len="5")])
    c2 = score_homology(t2, p2)["by_p"][("lib", "subs", 0.1, 0.0, 5)]
    assert c2.get("n_false_flank", 0) == 0
    assert c2["n_detected"] == 1 and c2["called_abs5"] == 0


def test_detection_requires_both_sides(tmp_path):
    tr = _truth_row(flank5="5", flank3="5", ltr5_start="6", ltr5_end="105",
                    ltr3_start="206", ltr3_end="305", seq_len="310")
    t, p = _write(tmp_path, [tr],
                  [_pred_row(seq_len="310", ltr5_start="6", ltr5_end="105",
                             ltr3_start="206", ltr3_end="310",
                             flank5_len="5", flank3_len="0")])
    c = score_homology(t, p)["by_p"][("lib", "subs", 0.1, 0.0, 5)]
    assert c.get("n_detected", 0) == 0
    # the 3' side was missed entirely, so its called length is 0 vs a true 5
    assert c["called_abs3"] == 5


def test_a_lost_record_still_pays_the_full_flank_error(tmp_path):
    """A configuration that answers less often must not look more accurate:
    a record with no call counts as having called a zero-length flank."""
    tr = _truth_row(flank5="25", flank3="25", ltr5_start="26", ltr5_end="125",
                    ltr3_start="226", ltr3_end="325", seq_len="350")
    t, p = _write(tmp_path, [tr], [_pred_row(status="no_pair", seq_len="350",
                                             ltr5_start="NA", ltr5_end="NA",
                                             ltr3_start="NA", ltr3_end="NA",
                                             flank5_len="NA", flank3_len="NA",
                                             k2p="NA")])
    c = score_homology(t, p)["by_p"][("lib", "subs", 0.1, 0.0, 25)]
    assert c["n_lost"] == 1 and c.get("n_located", 0) == 0
    assert c["called_abs5"] == 25 and c["called_abs3"] == 25


def test_weak_pair_counts_as_located_not_lost(tmp_path):
    """That is the whole point of the non-destructive gate: the measurement
    survives even though the pair did not clear significance."""
    t, p = _write(tmp_path, [_truth_row()], [_pred_row(status="weak_pair")])
    c = score_homology(t, p)["by_p"][KEY]
    assert c["n_weak"] == 1 and c.get("n_lost", 0) == 0 and c["n_located"] == 1
    assert c.get("n_pass", 0) == 0


def test_k2p_error_is_measured_against_the_exact_realized_value(tmp_path):
    t, p = _write(tmp_path, [_truth_row(realized_k2p="0.11", d_nominal="0.13")],
                  [_pred_row(k2p="0.15")])
    c = score_homology(t, p)["by_p"][KEY]
    assert c["k2p_err"] == pytest.approx(0.04)
    assert c["k2p_err_nom"] == pytest.approx(0.02)


def test_dhat_view_reuses_the_same_accounting(tmp_path):
    """The two views differ only in their key; a divergence between them would
    silently invalidate the schedule derivation."""
    t, p = _write(tmp_path, [_truth_row()], [_pred_row(ltr5_end="90")])
    h = tmp_path / "dhat.tsv"
    h.write_text("seq_id\td_hat\nx\t0.22\n")
    got = score_homology(t, p, h)
    assert got["by_dhat"] is not None
    a = got["by_p"][KEY]
    b = got["by_dhat"][("lib", "subs", dhat_bin(0.22), 0.0, 0)]
    assert a == b


def test_dhat_view_is_absent_without_a_dhat_file(tmp_path):
    t, p = _write(tmp_path, [_truth_row()], [_pred_row()])
    assert score_homology(t, p)["by_dhat"] is None


def test_cells_round_trip_through_json(tmp_path):
    t, p = _write(tmp_path, [_truth_row()], [_pred_row()])
    scored = score_homology(t, p)
    path = tmp_path / "cells.json"
    cells_to_json(scored, path)
    back = cells_from_json(path)
    assert back["by_p"] == scored["by_p"]
    assert back["by_dhat"] is None


def test_dhat_bin_edges_are_the_grid_midpoints():
    assert dhat_bin(0.0) == 0.0
    assert dhat_bin(0.024) == 0.0 and dhat_bin(0.026) == 0.05
    assert dhat_bin(0.074) == 0.05 and dhat_bin(0.076) == 0.1
    assert dhat_bin(9.0) == 0.5
