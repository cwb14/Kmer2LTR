# tests/test_build_truth.py
import pytest
from bench.build_truth import family_key, pair_library

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
