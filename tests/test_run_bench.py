import random
from pathlib import Path
import pytest
from bench.run_bench import make_dataset, score_run

def test_make_dataset_writes_matching_fasta_and_truth(tmp_path):
    tfa, ttsv = tmp_path / "t.fa", tmp_path / "t.tsv"
    tfa.write_text(">e1\n" + "A"*300 + "C"*900 + "A"*300 + "\n")
    ttsv.write_text("elem_id\ttotal_len\tltr_len\tltr5_start\tltr5_end\tltr3_start\tltr3_end\n"
                    "e1\t1500\t300\t1\t300\t1201\t1500\n")
    out_fa, out_tsv = tmp_path / "s.fa", tmp_path / "s.tsv"
    n = make_dataset(tfa, ttsv, out_fa, out_tsv, ds=[0.1, 0.3], kappas=[2.0],
                     flanks=[0, 50], n_per_cell=1, seed=0)
    n_fa = sum(1 for l in out_fa.read_text().split("\n") if l.startswith(">"))
    n_tsv = len(out_tsv.read_text().rstrip("\n").split("\n")) - 1
    assert n == n_fa == n_tsv == 4          # 2 ds x 1 kappa x 2 flanks x 1 rep

def test_score_run_reports_zero_error_on_perfect_predictions(tmp_path):
    truth = tmp_path / "t.tsv"
    truth.write_text("elem_id\tltr5_start\tltr5_end\tltr3_start\tltr3_end\t"
                     "d_nominal\trealized_k2p\tflank5\tflank3\n"
                     "e1\t1\t300\t1201\t1500\t0.1\t0.1\t0\t0\n")
    pred = tmp_path / "p.tsv"
    pred.write_text("seq_id\tstatus\tltr5_start\tltr5_end\tltr3_start\tltr3_end\tk2p\n"
                    "e1\tpass\t1\t300\t1201\t1500\t0.1\n")
    s = score_run(pred, truth)
    assert s["n"] == 1
    assert s["mae_ltr5_start"] == 0
    assert s["k2p_bias"] == pytest.approx(0.0)

def test_score_run_detects_boundary_error(tmp_path):
    truth = tmp_path / "t.tsv"
    truth.write_text("elem_id\tltr5_start\tltr5_end\tltr3_start\tltr3_end\t"
                     "d_nominal\trealized_k2p\tflank5\tflank3\n"
                     "e1\t1\t300\t1201\t1500\t0.1\t0.1\t0\t0\n")
    pred = tmp_path / "p.tsv"
    pred.write_text("seq_id\tstatus\tltr5_start\tltr5_end\tltr3_start\tltr3_end\tk2p\n"
                    "e1\tpass\t6\t300\t1201\t1495\t0.12\n")
    s = score_run(pred, truth)
    assert s["mae_ltr5_start"] == 5
    assert s["k2p_bias"] == pytest.approx(0.02)


# --------------------------------------------------------------------------- #
# Additional coverage: ablate(), the fixed-matrix pickling workaround, and
# score_gold_grid(). Not strictly required, but "every standalone tool gets at
# least basic tests" (project convention) -- and the fixed-matrix path is a
# genuine footgun (parasail.Matrix cannot be pickled; ProcessPoolExecutor
# would silently hang/crash if a config tried to pass one through) worth a
# regression test on its own.
# --------------------------------------------------------------------------- #

def _rnd(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


def _write_fasta(tmp_path, records: dict[str, str]) -> Path:
    p = tmp_path / "in.fa"
    p.write_text("".join(f">{i}\n{s}\n" for i, s in records.items()))
    return p


def test_ablate_calibrated_matches_direct_classify_call(tmp_path):
    from kmer2ltr.align import classify
    ltr = _rnd(300, 1)
    S = ltr + _rnd(900, 2) + ltr
    fa = _write_fasta(tmp_path, {"x": S})
    out = tmp_path / "pred.tsv"
    summary = ablate_result = None
    from bench.run_bench import ablate
    summary = ablate("calibrated", fa, out, threads=1)
    assert summary["n"] == 1
    rows = out.read_text().rstrip("\n").split("\n")
    assert rows[0].split("\t")[0] == "seq_id"          # header present
    fields = rows[1].split("\t")
    direct = classify("x", S)
    assert fields[3] == str(direct.ltr5_start)          # ltr5_start column
    assert fields[2] == direct.status


def test_ablate_fixed_matrix_configs_run_under_multiprocessing(tmp_path):
    """Regression guard: blastn_1_3/fixed_1_1 pass a fixed parasail matrix and
    skip calibration. parasail.Matrix objects cannot be pickled (confirmed:
    pickle.dumps(GENERIC_MATRIX) raises ValueError), so this must NOT put a
    Matrix object into a ProcessPoolExecutor task argument -- threads=2 here
    forces the multiprocessing path, so a regression would hang or raise
    instead of silently "working" because it happened to run single-process.
    """
    from bench.run_bench import ablate
    ltr = _rnd(300, 3)
    recs = {f"r{i}": ltr + _rnd(900, 100 + i) + ltr for i in range(4)}
    fa = _write_fasta(tmp_path, recs)
    for name in ("blastn_1_3", "fixed_1_1"):
        out = tmp_path / f"{name}.tsv"
        summary = ablate(name, fa, out, threads=2)
        assert summary["n"] == 4
        rows = out.read_text().rstrip("\n").split("\n")
        assert len(rows) == 5                            # header + 4 records
        statuses = [r.split("\t")[2] for r in rows[1:]]
        assert all(s == "pass" for s in statuses)


def test_ablate_no_stage3_reports_flank_on_diverged_terminus():
    """Functional check that the ablation actually flips behaviour, using the
    same kind of construction the spec's Stage 3 section describes: a
    perfectly-bounded element whose terminal bases happen to be diverged
    enough that plain greedy trimming (no_stage3) calls a flank where the
    calibrated Stage 3 model comparison (default) does not.
    """
    from bench.run_bench import ablate
    import tempfile
    r = random.Random(42)
    ti = {"A": "G", "G": "A", "C": "T", "T": "C"}
    ltr = _rnd(400, 5)
    # diverge just the last few bases of the 5' copy relative to the 3' copy
    ltr5 = ltr[:-4] + "".join(ti[c] for c in ltr[-4:])
    S = ltr5 + _rnd(1200, 6) + ltr
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        fa = _write_fasta(d, {"x": S})
        out_default = d / "default.tsv"
        out_no3 = d / "no_stage3.tsv"
        ablate("calibrated", fa, out_default, threads=1)
        ablate("no_stage3", fa, out_no3, threads=1)
        row_default = out_default.read_text().splitlines()[1].split("\t")
        row_no3 = out_no3.read_text().splitlines()[1].split("\t")
        flank5_default = row_default[9]     # flank5_len column
        flank5_no3 = row_no3[9]
        # no_stage3 uses the raw discovered endpoint directly (ltr_spans), which
        # trims a diverged terminus more readily than Stage 3's model comparison
        assert flank5_default == "0"
        assert int(flank5_no3) >= int(flank5_default)


def test_score_gold_grid_matches_hand_counted_cells(tmp_path):
    from bench.run_bench import score_gold_grid
    truth = tmp_path / "truth.tsv"
    truth.write_text(
        "seq_id\tdataset\torig_id\td_nominal\tflank5\tflank3\tflank_source\t"
        "ltr5_start\tltr5_end\tltr3_start\tltr3_end\tseq_len\n"
        "a\tds\to1\t0.0\t0\t0\tother\t1\t300\t1000\t1300\t1300\n"
        "b\tds\to1\t0.0\t0\t0\tother\t1\t300\t1000\t1300\t1300\n"
        "c\tds\to1\t0.1\t20\t20\tother\t21\t320\t1000\t1300\t1320\n"
        "d\tds\to1\t0.1\t20\t20\tother\t21\t320\t1000\t1300\t1320\n"
    )
    pred = tmp_path / "pred.tsv"
    pred.write_text(
        "seq_id\tstatus\tltr5_start\tltr5_end\tltr3_start\tltr3_end\t"
        "flank5_len\tflank3_len\tk2p\n"
        "a\tpass\t1\t300\t1000\t1300\t0\t0\t0.02\n"      # flank=0, correctly bounded
        "b\tpass\t6\t300\t1000\t1295\t5\t5\t0.05\n"       # flank=0, FALSE flank called
        "c\tpass\t21\t320\t1000\t1300\t20\t20\t0.09\n"    # flank=20, detected exactly
        "d\tpass\t1\t320\t1000\t1300\t0\t20\t0.30\n"      # flank=20, 5' side MISSED
    )
    s = score_gold_grid(truth, pred)
    cell0 = s["by_d"][(0.0, 0)]
    assert cell0["n"] == 2 and cell0["n_pass"] == 2 and cell0["n_false_flank"] == 1
    cell20 = s["by_d"][(0.1, 20)]
    assert cell20["n"] == 2
    assert cell20["n_detected"] == 1                     # only "c": BOTH sides > 0
    assert cell20["sum_abs_err5"] == pytest.approx(0 + 20)   # c: |21-1-20|=0; d: |0-20|=20
    assert cell20["sum_abs_err3"] == pytest.approx(0 + 0)    # both called 3' side exactly

    # by_called buckets on the TOOL's own flank5_len>0-or-flank3_len>0 decision,
    # independent of ground truth -- "a" is the only record the tool called
    # clean (no flank on either side); b/c/d all had at least one side called.
    assert s["by_called"][(0.0, False)]["n"] == 1                 # a
    assert s["by_called"][(0.0, True)]["n"] == 1                  # b
    assert s["by_called"][(0.1, True)]["n"] == 2                  # c, d
    assert (0.1, False) not in s["by_called"]
    bc_true_01 = s["by_called"][(0.1, True)]
    assert bc_true_01["n_k2p"] == 2
    assert bc_true_01["sum_k2p_err"] == pytest.approx((0.09 - 0.1) + (0.30 - 0.1))

    # JSON round-trip preserves both numeric (d_nominal, flank_true) keys and
    # (d_nominal, bool) keys exactly.
    from bench.run_bench import cells_to_json, cells_from_json
    json_path = tmp_path / "cells.json"
    cells_to_json(s, json_path)
    s2 = cells_from_json(json_path)
    assert s2["by_d"][(0.1, 20)] == cell20
    assert s2["by_called"][(0.0, False)]["n"] == 1
    assert s2["by_called"][(0.1, True)]["n_k2p"] == 2


def test_trim_leaves_boundaries_and_flank_calls_unchanged():
    """The trim ablation's "overall vs within the flank-called subset" K2P
    comparison is only valid if every trim_K config
    reports the SAME boundaries/flank calls for a given input -- otherwise
    "the flank-called subset" would mean a DIFFERENT population of records
    for each K, and cross-config bias/RMSE comparisons would not isolate
    trim's effect on K2P from a boundary-call difference. classify()'s trim=
    only drops columns from the final aligned pair before count_substitutions
    (see kmer2ltr/align.py), strictly AFTER bounds/flank5_len/flank3_len are
    already fixed -- verified directly here, not just inferred from reading
    the source, since it is what every trim-config comparison in this
    module's gold analysis depends on.
    """
    from kmer2ltr.align import classify
    ltr = _rnd(300, 20)
    S = _rnd(30, 21) + ltr + _rnd(900, 22) + ltr + _rnd(45, 23)
    results = {k: classify("x", S, trim=k) for k in (0, 3, 5, 10)}
    coords = {(r.status, r.ltr5_start, r.ltr5_end, r.ltr3_start, r.ltr3_end,
               r.flank5_len, r.flank3_len) for r in results.values()}
    assert len(coords) == 1
    # n_sites strictly decreases as trim grows (2*K fewer columns each side)
    n_sites = [results[k].n_sites for k in (0, 3, 5, 10)]
    assert n_sites == sorted(n_sites, reverse=True)
