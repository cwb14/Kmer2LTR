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

def test_sanitize_covers_whitespace_classified_control_codes():
    # 0x1c-0x1f, 0x85, 0xa0 are isspace()==True and were previously unmapped
    for ch in "\x1c\x1d\x1e\x1f\x85\xa0":
        assert sanitize(f"AC{ch}GT") == "ACGT", f"{ch!r} not handled"


def test_sanitize_maps_non_ascii_to_N():
    assert sanitize("AC中GT") == "ACNGT"       # CJK
    assert sanitize("AC\U0001F600GT") == "ACNGT"   # emoji (astral plane)
    assert sanitize("ACGTé") == "ACGTN"       # accented Latin


def test_read_fasta_raw_indexes_identically_to_the_sanitized_sequence():
    """Coordinates are measured on the sanitized copy and used to slice the raw
    one, so the two must line up character for character."""
    import tempfile, os
    from ltrk2p.fasta import read_fasta_raw
    fd, path = tempfile.mkstemp(suffix=".fa")
    os.close(fd)
    try:
        open(path, "w").write(">a desc\nacgtRY nn\nAC\tGT\n>b\n\n")
        recs = list(read_fasta_raw(path))
    finally:
        os.unlink(path)
    (sid, clean, raw), (bid, bclean, braw) = recs
    assert sid == "a"
    assert raw == "acgtRYnnACGT"
    assert clean == "ACGTNNNNACGT"
    assert len(raw) == len(clean)
    assert (bid, bclean, braw) == ("b", "", "")


def test_read_fasta_raw_preserves_case_and_ambiguity_codes():
    """Soft-masking and IUPAC codes are information the user supplied. The
    sliced FASTA outputs are written from this string precisely so that a
    boundary-corrected element is still the user's own sequence."""
    import tempfile, os
    from ltrk2p.fasta import read_fasta_raw
    fd, path = tempfile.mkstemp(suffix=".fa")
    os.close(fd)
    try:
        open(path, "w").write(">x\nacgtACGTrynN-*\n")
        (_id, clean, raw), = list(read_fasta_raw(path))
    finally:
        os.unlink(path)
    assert raw == "acgtACGTrynN-*"
    assert clean == "ACGTACGTNNNNNN"


def test_read_fasta_is_the_raw_reader_without_the_raw_column():
    import tempfile, os
    from ltrk2p.fasta import read_fasta, read_fasta_raw
    fd, path = tempfile.mkstemp(suffix=".fa")
    os.close(fd)
    try:
        open(path, "w").write(">a\nacgt\n>b\nTTTT\n")
        assert list(read_fasta(path)) == [(i, s) for i, s, _ in read_fasta_raw(path)]
    finally:
        os.unlink(path)


def test_despace_and_sanitize_agree_on_whitespace_for_every_codepoint():
    """The load-bearing invariant behind every sliced FASTA output.

    A coordinate is measured on the sanitized sequence and used to slice the raw
    one, so the two must have the same length and the same character-by-character
    correspondence. That holds only if `_despace` removes exactly the characters
    `sanitize` drops -- one uses `str.split()`, the other `str.isspace()`, and a
    single codepoint where they disagree would silently shift every downstream
    slice by one base with no error anywhere.
    """
    from ltrk2p.fasta import _despace, sanitize
    disagree = [cp for cp in range(0x110000)
                if not 0xD800 <= cp < 0xE000        # lone surrogates: not text
                and len(sanitize(_despace(chr(cp)))) != len(_despace(chr(cp)))]
    assert disagree == [], f"{len(disagree)} codepoints disagree, first: {disagree[:5]}"


def test_raw_and_sanitized_line_up_through_the_reader():
    import random, tempfile, os
    from ltrk2p.fasta import read_fasta_raw
    r = random.Random(0)
    pool = "ACGTacgtRYNn-*0 \t\v\f 　x"
    body = "".join(r.choice(pool) for _ in range(4000))
    fd, path = tempfile.mkstemp(suffix=".fa")
    os.close(fd)
    try:
        wrapped = "\n".join(body[i:i + 61] for i in range(0, len(body), 61))
        open(path, "w").write(f">x\n{wrapped}\n")
        (_id, clean, raw), = list(read_fasta_raw(path))
    finally:
        os.unlink(path)
    assert len(clean) == len(raw)
    assert all((c == "N") or (c == rc.upper()) for c, rc in zip(clean, raw))
