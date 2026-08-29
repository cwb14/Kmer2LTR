# tests/test_gold_subset.py
import random
import pytest
from bench.gold_subset import select_gold, perturb, flank_from

def _rnd(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))

def _write(tmp_path, rows, seqs):
    p = tmp_path / "pred.tsv"
    hdr = ("seq_id\tseq_len\tstatus\tltr5_start\tltr5_end\tltr3_start\tltr3_end\t"
           "ltr5_len\tltr3_len\tk2p\tbitscore\n")
    p.write_text(hdr + "".join("\t".join(map(str, r)) + "\n" for r in rows))
    f = tmp_path / "in.fa"
    f.write_text("".join(f">{i}\n{s}\n" for i, s in seqs.items()))
    return p, f

def test_select_gold_requires_terminal_boundaries(tmp_path):
    ltr = "TG" + _rnd(296, 1) + "CA"
    elem = ltr + _rnd(500, 2) + ltr
    L = len(elem)
    good = ("g", L, "pass", 1, 300, L - 299, L, 300, 300, 0.01, 500.0)
    # overextended -> excluded even though everything else is clean
    bad = ("b", L, "pass", 5, 304, L - 299, L, 300, 300, 0.01, 500.0)
    pred, fa = _write(tmp_path, [good, bad], {"g": elem, "b": elem})
    assert [i for i, _ in select_gold(pred, fa)] == ["g"]

def test_select_gold_excludes_high_divergence_and_low_score(tmp_path):
    ltr = "TG" + _rnd(296, 3) + "CA"
    elem = ltr + _rnd(500, 4) + ltr
    L = len(elem)
    rows = [("hi_k2p", L, "pass", 1, 300, L - 299, L, 300, 300, 0.40, 500.0),
            ("lo_bits", L, "pass", 1, 300, L - 299, L, 300, 300, 0.01, 10.0),
            ("ok", L, "pass", 1, 300, L - 299, L, 300, 300, 0.01, 500.0)]
    pred, fa = _write(tmp_path, rows, {r[0]: elem for r in rows})
    assert [i for i, _ in select_gold(pred, fa)] == ["ok"]

def test_select_gold_motif_requirement_is_optional(tmp_path):
    ltr = _rnd(300, 5)                      # no TG..CA
    elem = ltr + _rnd(500, 6) + ltr
    L = len(elem)
    row = ("nomotif", L, "pass", 1, 300, L - 299, L, 300, 300, 0.01, 500.0)
    pred, fa = _write(tmp_path, [row], {"nomotif": elem})
    assert [i for i, _ in select_gold(pred, fa, require_motif=True)] == []
    assert [i for i, _ in select_gold(pred, fa, require_motif=False)] == ["nomotif"]

def test_select_gold_excludes_non_pass_status(tmp_path):
    elem = _rnd(1000, 7)
    row = ("np", 1000, "no_pair", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
    pred, fa = _write(tmp_path, [row], {"np": elem})
    assert list(select_gold(pred, fa)) == []

def test_perturb_truth_coordinates_are_exact():
    ltr = "TG" + _rnd(296, 8) + "CA"
    elem = ltr + _rnd(600, 9) + ltr
    seq, truth = perturb(elem, 0.0, 2.0, 40, 60, random.Random(10), lambda n, r: "A" * n)
    assert truth["ltr5_start"] == 41
    assert truth["ltr3_end"] == len(seq) - 60
    assert seq[40:40 + 300] == ltr          # unmutated at d=0
    assert len(seq) == len(elem) + 100

def test_perturb_leaves_internal_region_untouched():
    ltr = _rnd(200, 11)
    internal = _rnd(700, 12)
    elem = ltr + internal + ltr
    seq, truth = perturb(elem, 0.3, 2.0, 0, 0, random.Random(13), lambda n, r: "")
    # the internal region must be byte-identical; only the LTRs are evolved
    assert internal in seq

def test_flank_from_shuffle_preserves_composition():
    from bench.gold_subset import flank_from
    pool = _rnd(4000, 14)
    f = flank_from(pool, 300, random.Random(15), mode="shuffle")
    assert len(f) == 300
    assert set(f) <= set("ACGT")

def test_perturb_zero_flank_is_perfectly_bounded():
    ltr = _rnd(250, 16)
    elem = ltr + _rnd(600, 17) + ltr
    seq, truth = perturb(elem, 0.1, 2.0, 0, 0, random.Random(18), lambda n, r: "")
    assert truth["ltr5_start"] == 1
    assert truth["ltr3_end"] == len(seq)
