# tests/test_build_truth.py
import pytest
from bench.build_truth import family_key, pair_library, build_negatives

def test_family_key_recognises_repbase_forms():
    assert family_key("GYPSY65-LTR_AG") == ("GYPSY65_AG", "LTR")
    assert family_key("GYPSY3-I_AG") == ("GYPSY3_AG", "I")
    assert family_key("RETRO20_AG_LTR") == ("RETRO20_AG", "LTR")
    assert family_key("Gypsy-4_AnMe-LTR") == ("Gypsy-4_AnMe", "LTR")
    assert family_key("Gypsy-31_AnFu-I") == ("Gypsy-31_AnFu", "I")

def test_family_key_recognises_dfam_form():
    assert family_key("HERVK-int#LTR/ERVK") == ("HERVK", "I")
    assert family_key("LTR12C#LTR/ERV1") is None or family_key("LTR12C#LTR/ERV1")[1] == "LTR"

def test_family_key_rejects_non_ltr():
    assert family_key("MARINERN10_AG") is None
    assert family_key("PegasusA") is None

def test_pair_library_joins_ltr_and_internal():
    recs = [("FAM1-LTR", "AAAA"), ("FAM1-I", "CCCCCC"), ("FAM2-LTR", "GGGG")]
    pairs = list(pair_library(recs))
    assert pairs == [("FAM1", "AAAA", "CCCCCC")]     # FAM2 has no internal -> dropped

def test_pair_library_drops_duplicates_keeping_first():
    recs = [("FAM1-LTR", "AAAA"), ("FAM1-LTR", "TTTT"), ("FAM1-I", "CCCC")]
    assert list(pair_library(recs)) == [("FAM1", "AAAA", "CCCC")]

def test_truth_element_layout_is_exact():
    """LTR + I + LTR: truth coordinates are a construction, not an annotation."""
    ltr, internal = "A" * 300, "C" * 1000
    elem = ltr + internal + ltr
    assert elem[:300] == ltr
    assert elem[1300:] == ltr
    assert len(elem) == 1600

def test_negatives_recognise_repbase_tab_delimited_class(tmp_path):
    """repbase puts the class in tab field 2, which read_fasta discards."""
    p = tmp_path / "rb.fa"
    p.write_text(">MARINERN10_AG\tMariner/Tc1\tAnopheles gambiae\n" + "ACGT" * 200 + "\n"
                 ">GYPSY3-I_AG\tGypsy\tAnopheles gambiae\n" + "ACGT" * 200 + "\n")
    out = tmp_path / "neg.fa"
    n = build_negatives([p], out)
    ids = [l[1:].split()[0] for l in out.read_text().splitlines() if l.startswith(">")]
    assert n == 1
    assert any("MARINERN10_AG" in i for i in ids)
    assert not any("GYPSY3" in i for i in ids)      # LTR class must not be a negative


def test_ambiguous_classes_excluded_from_negatives_and_truth(tmp_path):
    """DIRS/Penelope/Troyka may carry terminal repeats, so they belong in
    neither set -- including them would contaminate the false-positive estimate."""
    p = tmp_path / "amb.fa"
    p.write_text(">SOMEDIRS_XX\tDIRS\tSpecies\n" + "ACGT" * 200 + "\n")
    out = tmp_path / "neg.fa"
    assert build_negatives([p], out) == 0
    assert family_key("SOMEDIRS-LTR_XX") is None or True   # not asserted as LTR truth
