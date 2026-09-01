import random
from dataclasses import replace

import pytest

from kmer2ltr.align import _classify, classify
from kmer2ltr.extras import (ExtraSpec, ExtraWriter, build, fasta_record,
                           iupac_consensus, shift_locus)


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


def _element(seed=1, flank5=0, flank3=0, d=0.10):
    """A synthetic element plus the raw string the extras are cut from."""
    ltr = "TG" + _rnd(296, seed) + "CA"
    core = ltr + _rnd(1000, 100 + seed) + _evolve(ltr, d, 200 + seed)
    return _rnd(flank5, 300 + seed) + core + _rnd(flank3, 400 + seed)


# --------------------------------------------------------------------------- #
# IUPAC consensus
# --------------------------------------------------------------------------- #

def test_consensus_keeps_agreed_bases_and_codes_disagreements():
    assert iupac_consensus("ACGT", "ACGT") == "ACGT"
    assert iupac_consensus("AC", "GT") == "RY"
    assert iupac_consensus("GAGA", "CTTC") == "SWKM"


def test_consensus_takes_the_observed_base_across_a_gap():
    """An indel in one copy is not evidence about the other copy's base.

    Dropping the column instead would silently shorten every consensus by the
    element's indel load, which is exactly the elements a family clustering
    most needs to keep whole.
    """
    assert iupac_consensus("A-CT", "AGCT") == "AGCT"
    assert iupac_consensus("AGCT", "A-CT") == "AGCT"


def test_consensus_maps_anything_unreadable_to_N():
    assert iupac_consensus("ANCN", "ATCG") == "ANCN"
    assert iupac_consensus("N-", "-N") == "NN"


def test_consensus_rejects_unequal_lengths():
    """A length mismatch means the caller passed something that is not an
    alignment; silently zipping to the shorter one would emit a plausible
    consensus for a pair that was never aligned."""
    with pytest.raises(ValueError):
        iupac_consensus("ACGT", "ACG")


def test_consensus_is_as_long_as_the_alignment():
    ltr = "TG" + _rnd(296, 5) + "CA"
    S = ltr + _rnd(900, 6) + _evolve(ltr, 0.12, 7)
    _r, aln = _classify("x", S)
    cons = iupac_consensus(*aln)
    assert len(cons) == len(aln[0]) == len(aln[1])


# --------------------------------------------------------------------------- #
# Header coordinate correction
# --------------------------------------------------------------------------- #

def test_shift_locus_moves_both_ends_inwards():
    assert shift_locus("chr1:1000-2000", 40, 25) == "chr1:1040-1975"


def test_shift_locus_preserves_the_separator_dialect_and_the_class_tag():
    """LTR_retriever writes `start..end`; the arabidopsis set writes `start-end`.
    Rewriting one as the other would break whatever produced the header."""
    assert shift_locus("chr1:1000..2000#LTR/Gypsy", 5, 5) == "chr1:1005..1995#LTR/Gypsy"
    assert shift_locus("chr1:1000-2000#LTR/Copia", 0, 0) == "chr1:1000-2000#LTR/Copia"


def test_shift_locus_binds_the_last_colon_so_colons_in_chrom_names_survive():
    assert shift_locus("scaffold:a:100-200", 10, 10) == "scaffold:a:110-190"


@pytest.mark.parametrize("seq_id", [
    "elem00042",                 # no locus at all
    "chr1:1000",                 # no range
    "chr1:1000-2000-3000",       # ambiguous three-part range
    "chr1:abc-def",              # non-numeric
])
def test_shift_locus_declines_rather_than_guesses(seq_id):
    """This is the one place Kmer2LTR reads a header as anything but an opaque id.
    A header it cannot parse exactly is returned untouched."""
    assert shift_locus(seq_id, 10, 10) == seq_id


def test_shift_locus_declines_when_the_trim_would_invert_the_interval():
    assert shift_locus("chr1:100-200", 500, 0) == "chr1:100-200"


def test_fasta_record_wraps_and_survives_an_empty_sequence():
    assert fasta_record("h", "ACGTAC", wrap=4) == ">h\nACGT\nAC\n"
    assert fasta_record("h", "") == ">h\n\n"


# --------------------------------------------------------------------------- #
# Per-record extras
# --------------------------------------------------------------------------- #

_ALL = ExtraSpec(consensus=True, internal=True, trimmed=True,
                 perfect=("5p", "3p", "consensus"))


def test_build_declines_anything_but_a_pass():
    """weak_pair coordinates are real, but these outputs assert that the element
    IS one -- the claim the significance gate just declined to support."""
    S = _element(2)
    result, aln = _classify("x", S)
    assert result.status == "pass"
    for status in ("weak_pair", "k2p_undefined", "no_pair"):
        assert build(replace(result, status=status), aln, S, _ALL) is None


def test_trimmed_record_is_the_input_sliced_at_the_called_boundaries():
    S = _element(3, flank5=40, flank3=25)
    result, aln = _classify("chr5:1000-{}".format(999 + len(S)), S)
    e = build(result, aln, S, ExtraSpec(trimmed=True))
    body = "".join(e.trimmed.split("\n")[1:])
    assert body == S[result.ltr5_start - 1:result.ltr3_end]
    assert e.trimmed.startswith(
        f">chr5:{1000 + result.flank5_len}-{999 + len(S) - result.flank3_len}\n")


def test_slices_of_the_input_keep_the_input_s_own_characters():
    """Soft-masking and IUPAC codes must survive into every sliced output.

    The tool measures on a sanitised copy where lowercase is uppercased and
    every ambiguity code becomes N. Writing THAT back out as "your element,
    trimmed" would hand the user a worse file than the one they supplied.
    """
    S = _element(4, flank5=30, flank3=30)
    raw = S[:400].lower() + "R" + S[401:]          # soft-mask a stretch, add a code
    result, aln = _classify("x", S)
    e = build(result, aln, raw, ExtraSpec(trimmed=True, internal=True,
                                          perfect=("5p",)))
    trimmed = "".join(e.trimmed.split("\n")[1:])
    assert trimmed == raw[result.ltr5_start - 1:result.ltr3_end]
    assert trimmed != trimmed.upper(), "lowercase was lost"
    assert "R" in "".join(e.internal.split("\n")[1:]) or "R" in trimmed


def test_perfect_records_carry_one_ltr_copy_on_both_ends():
    S = _element(5)
    result, aln = _classify("x", S)
    e = build(result, aln, S, ExtraSpec(perfect=("5p", "3p", "consensus")))
    ltr5 = S[result.ltr5_start - 1:result.ltr5_end]
    ltr3 = S[result.ltr3_start - 1:result.ltr3_end]
    internal = S[result.ltr5_end:result.ltr3_start - 1]
    for mode, flank in (("5p", ltr5), ("3p", ltr3)):
        header, *body = e.perfect[mode].split("\n")
        assert "".join(body) == flank + internal + flank
        assert header == f">x~LTRlen:{len(flank)}"
    assert e.perfect["consensus"].startswith(">x~LTRlen:")


def test_consensus_is_built_even_when_only_a_perfect_record_needs_it():
    """`--perfect-ltr-rt consensus` without `--ltr-cluster` still needs the
    consensus sequence, but must not write the consensus FASTA."""
    S = _element(6)
    result, aln = _classify("x", S)
    e = build(result, aln, S, ExtraSpec(perfect=("consensus",)))
    assert e.consensus is None
    assert e.perfect["consensus"]


def test_an_element_with_no_internal_region_yields_no_perfect_record():
    S = _element(7)
    result, aln = _classify("x", S)
    abutting = replace(result, ltr3_start=result.ltr5_end + 1)
    assert build(abutting, aln, S, ExtraSpec(perfect=("5p",))).perfect == {}


def test_spec_is_falsy_when_nothing_was_asked_for():
    assert not ExtraSpec()
    assert ExtraSpec(trimmed=True)
    assert not ExtraSpec(trimmed=True).wants_consensus
    assert ExtraSpec(perfect=("consensus",)).wants_consensus


# --------------------------------------------------------------------------- #
# ExtraWriter
# --------------------------------------------------------------------------- #

def test_writer_opens_one_stream_per_requested_output(tmp_path):
    spec = ExtraSpec(consensus=True, internal=True, trimmed=True, perfect=("5p",))
    with ExtraWriter(str(tmp_path / "out"), spec) as w:
        pass
    assert w.paths == {}, "streams that received nothing should not survive"
    assert list(tmp_path.iterdir()) == [], f"left behind {list(tmp_path.iterdir())}"


def test_writer_keeps_only_the_streams_that_received_records(tmp_path):
    """An input with no passing element should leave no misleading zero-byte
    FASTA to be mistaken for a real, empty result."""
    S = _element(8)
    result, aln = _classify("x", S)
    spec = ExtraSpec(consensus=True, trimmed=True)
    w = ExtraWriter(str(tmp_path / "out"), spec)
    w.write(build(result, aln, S, ExtraSpec(trimmed=True)))   # consensus stays empty
    w.close()
    assert set(w.paths) == {"trimmed"}
    assert not (tmp_path / "out.consensus.fa").exists()
    assert (tmp_path / "out.trimmed.fa").read_text().startswith(">x\n")


def test_writer_warns_once_on_a_duplicate_id(tmp_path, capsys):
    """The TSV is positional so duplicate ids are harmless there. These files
    are keyed by id, and an mmseqs cluster naming a duplicated id cannot be
    joined back to one element -- so say so, exactly once."""
    S = _element(9)
    result, aln = _classify("dup", S)
    e = build(result, aln, S, ExtraSpec(trimmed=True))
    with ExtraWriter(str(tmp_path / "out"), ExtraSpec(trimmed=True)) as w:
        for _ in range(3):
            w.write(e)
    err = capsys.readouterr().err
    assert err.count("duplicate record id") == 1, err
    assert w.counts["trimmed"] == 3, "the records themselves must still be written"


def test_writer_ignores_records_from_a_non_passing_element(tmp_path):
    with ExtraWriter(str(tmp_path / "out"), ExtraSpec(trimmed=True)) as w:
        w.write(None)
    assert w.counts["trimmed"] == 0


# --------------------------------------------------------------------------- #
# Regressions
# --------------------------------------------------------------------------- #

def test_shift_locus_finds_a_locus_after_a_bedtools_name_prefix():
    """`bedtools getfasta -name` on a RepeatMasker-named BED writes
    `TE_1#LTR/Copia::chr1:1001-4860` -- the `#` is in the MIDDLE and the locus
    is past it. Reading only up to the first `#` finds nothing, and a header
    left unchanged over a sequence that WAS trimmed lies about its own span."""
    assert (shift_locus("TE_1#LTR/Copia::chr1:1001-4860", 30, 30)
            == "TE_1#LTR/Copia::chr1:1031-4830")
    assert shift_locus("name::chr1:100..200", 5, 5) == "name::chr1:105..195"


def test_shift_locus_does_not_swallow_a_trailing_newline():
    """`$` also matches just before a trailing newline; `\\Z` does not. Not
    reachable through `read_fasta_raw`, but a rewrite that silently eats a
    character is not a behaviour to leave available."""
    assert shift_locus("chr1:1000-2000\n", 10, 10) == "chr1:1000-2000\n"


def test_writer_never_deletes_a_file_it_did_not_create(tmp_path):
    """`close()` removes streams that received nothing. It must distinguish
    "empty file I just made" from "file that was already there and I truncated":
    deleting the user's file is a worse outcome than leaving an empty one."""
    pre = tmp_path / "out.trimmed.fa"
    pre.write_text(">existing\nACGT\n")
    with ExtraWriter(str(tmp_path / "out"), ExtraSpec(trimmed=True)) as w:
        pass
    assert pre.exists(), "an existing file was deleted"
    assert pre.read_text() == "", "the run wrote nothing, so it should read empty"
    assert w.counts["trimmed"] == 0


def test_writer_closes_what_it_opened_when_a_later_stream_fails(tmp_path):
    """A half-open writer is never handed back to the caller, so nothing else
    can close it. Reachable through EACCES, ENOSPC or EMFILE as well as this."""
    (tmp_path / "out.internal.fa").mkdir()
    spec = ExtraSpec(consensus=True, internal=True)
    with pytest.raises(OSError):
        ExtraWriter(str(tmp_path / "out"), spec)
    # The consensus stream opened first; it must not be left dangling.
    import gc
    gc.collect()
    assert (tmp_path / "out.consensus.fa").exists()


def test_stream_paths_lists_every_file_the_writer_would_open(tmp_path):
    """The CLI checks this set against the input before anything is opened,
    so it must not drift from what __init__ actually opens."""
    from kmer2ltr.extras import stream_paths
    spec = ExtraSpec(consensus=True, internal=True, trimmed=True, perfect=("5p", "3p"))
    declared = stream_paths(str(tmp_path / "o"), spec)
    with ExtraWriter(str(tmp_path / "o"), spec) as w:
        opened = dict(w.paths)
    assert declared == opened
