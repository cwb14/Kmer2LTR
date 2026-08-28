import gzip
import pytest
from ltrk2p.fasta import sanitize, read_fasta

def test_sanitize_uppercases_and_strips_whitespace():
    assert sanitize("acgt acgt\tacgt") == "ACGTACGTACGT"

def test_sanitize_maps_iupac_and_junk_to_N():
    assert sanitize("ACGTRYKMSWBDHVN") == "ACGT" + "N" * 11
    assert sanitize("AC-GT*1") == "ACNGTNN"

def test_read_fasta_wrapped_and_unwrapped(tmp_path):
    p = tmp_path / "a.fa"
    p.write_text(">one desc here\nACGT\nACGT\n>two\nGGGG\n")
    assert list(read_fasta(p)) == [("one", "ACGTACGT"), ("two", "GGGG")]

def test_read_fasta_id_truncated_at_first_whitespace(tmp_path):
    # real headers: repbase is tab-delimited, dfam has spaces, ltrharvest has none
    p = tmp_path / "b.fa"
    p.write_text(">BEL-76_AnFu-I\tBEL\tAnopheles funestus\nACGT\n"
                 ">5S#rRNA @Vertebrata [S:40,50]\nGGGG\n"
                 ">LR999451.1:357312-364998#LTR/unknown\nTTTT\n")
    assert [i for i, _ in read_fasta(p)] == [
        "BEL-76_AnFu-I", "5S#rRNA", "LR999451.1:357312-364998#LTR/unknown"]

def test_read_fasta_gzip(tmp_path):
    p = tmp_path / "c.fa.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(">x\nacgtACGT\n")
    assert list(read_fasta(p)) == [("x", "ACGTACGT")]

def test_read_fasta_blank_lines_and_empty_file(tmp_path):
    p = tmp_path / "d.fa"
    p.write_text(">x\n\nACGT\n\n\n>y\nGG\n")
    assert list(read_fasta(p)) == [("x", "ACGT"), ("y", "GG")]
    e = tmp_path / "e.fa"
    e.write_text("")
    assert list(read_fasta(e)) == []

def test_read_fasta_keeps_duplicate_ids(tmp_path):
    p = tmp_path / "f.fa"
    p.write_text(">dup\nAAAA\n>dup\nCCCC\n")
    assert list(read_fasta(p)) == [("dup", "AAAA"), ("dup", "CCCC")]

def test_read_fasta_record_with_no_sequence(tmp_path):
    p = tmp_path / "g.fa"
    p.write_text(">empty\n>next\nACGT\n")
    assert list(read_fasta(p)) == [("empty", ""), ("next", "ACGT")]
