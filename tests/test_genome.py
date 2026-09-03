import random

import pytest

from dataclasses import fields

from kmer2ltr.align import Result
from kmer2ltr.genome import (Context, Options, Window, annotate, credit,
                             element_context, find_tsd, harvest, locus,
                             orient, revcomp)


def _rnd(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


def _genome(tmp_path, name="g.fa", wrap=10, **contigs):
    """A tiny reference FASTA. Wrapping is deliberately narrow so that every
    harvested window spans several lines and the streaming offsets are exercised."""
    p = tmp_path / name
    with open(p, "w") as fh:
        for cid, seq in contigs.items():
            fh.write(f">{cid} some description here\n")
            for i in range(0, len(seq), wrap):
                fh.write(seq[i:i + wrap] + "\n")
    return p


# --------------------------------------------------------------------------- #
# Locus parsing
# --------------------------------------------------------------------------- #

def test_locus_reads_both_placements_and_the_two_separators():
    assert locus("chr1:1000-2000#LTR/Gypsy") == ("chr1", 1000, 2000)
    assert locus("TE_1#LTR/Copia::chr1:1000-2000") == ("chr1", 1000, 2000)
    assert locus("chr1:1000..2000") == ("chr1", 1000, 2000)


def test_locus_keeps_a_colon_inside_a_contig_name():
    assert locus("scaf:2:1000-2000") == ("scaf:2", 1000, 2000)


@pytest.mark.parametrize("seq_id", ["just_an_id", "chr1:1000", "chr1:a-b",
                                    "chr1:1000-2000-3000", ""])
def test_locus_declines_rather_than_guesses(seq_id):
    assert locus(seq_id) is None


# --------------------------------------------------------------------------- #
# Harvest
# --------------------------------------------------------------------------- #

def test_harvest_cuts_the_four_windows_at_the_right_offsets(tmp_path):
    seq = _rnd(200, 1)
    g = _genome(tmp_path, c1=seq)
    w = harvest([g], [("c1", 51, 150)], pad=5, probe=8)[("c1", 51, 150)]
    assert w == Window(up=seq[45:50], head=seq[50:58],
                       tail=seq[142:150], down=seq[150:155])


def test_harvest_clips_at_the_contig_edges(tmp_path):
    g = _genome(tmp_path, c1="ACGTACGTAC")
    w = harvest([g], [("c1", 1, 10)], pad=5, probe=4)[("c1", 1, 10)]
    assert w == Window(up="", head="ACGT", tail="GTAC", down="")


def test_harvest_uppercases_a_soft_masked_reference(tmp_path):
    g = _genome(tmp_path, c1="acgtACGTacgtACGT")
    w = harvest([g], [("c1", 5, 12)], pad=4, probe=4)[("c1", 5, 12)]
    assert w == Window(up="ACGT", head="ACGT", tail="ACGT", down="ACGT")


def test_harvest_skips_a_locus_whose_contig_the_reference_lacks(tmp_path):
    g = _genome(tmp_path, c1="ACGTACGTACGT")
    assert harvest([g], [("nope", 1, 4)]) == {}


def test_harvest_reads_several_reference_files(tmp_path):
    a = _genome(tmp_path, "a.fa", c1="ACGTACGTACGTACGT")
    b = _genome(tmp_path, "b.fa", c2="TTTTGGGGCCCCAAAA")
    got = harvest([a, b], [("c1", 5, 8), ("c2", 5, 8)], pad=2, probe=2)
    assert got[("c1", 5, 8)].head == "AC"
    assert got[("c2", 5, 8)].head == "GG"


def test_harvest_handles_overlapping_and_nested_requests(tmp_path):
    """The cuts are opened and closed against one rolling offset, so a locus
    contained inside another must not close the outer one early."""
    seq = _rnd(140, 2)
    g = _genome(tmp_path, c1=seq)
    got = harvest([g], [("c1", 10, 100), ("c1", 20, 40), ("c1", 30, 110)],
                  pad=4, probe=6)
    assert got[("c1", 10, 100)] == Window(seq[5:9], seq[9:15], seq[94:100], seq[100:104])
    assert got[("c1", 20, 40)] == Window(seq[15:19], seq[19:25], seq[34:40], seq[40:44])
    assert got[("c1", 30, 110)] == Window(seq[25:29], seq[29:35], seq[104:110], seq[110:114])


def test_harvest_reads_a_gzipped_reference(tmp_path):
    import gzip
    seq = _rnd(100, 3)
    p = tmp_path / "g.fa.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(">c1\n" + seq + "\n")
    w = harvest([p], [("c1", 21, 60)], pad=4, probe=4)[("c1", 21, 60)]
    assert w == Window(seq[16:20], seq[20:24], seq[56:60], seq[60:64])


def test_harvest_head_and_tail_may_overlap_on_a_short_element(tmp_path):
    """`head` and `tail` are separate cuts, so an element shorter than 2*probe
    simply has them overlap rather than the second being truncated."""
    seq = _rnd(60, 4)
    g = _genome(tmp_path, c1=seq)
    w = harvest([g], [("c1", 21, 30)], pad=3, probe=8)[("c1", 21, 30)]
    assert w == Window(seq[17:20], seq[20:28], seq[22:30], seq[30:33])


def test_harvest_ignores_a_contig_nothing_asks_about(tmp_path):
    seq = _rnd(80, 5)
    g = _genome(tmp_path, c1="A" * 80, c2=seq)
    got = harvest([g], [("c2", 21, 40)], pad=2, probe=4)
    assert list(got) == [("c2", 21, 40)]
    assert got[("c2", 21, 40)].head == seq[20:24]


# --------------------------------------------------------------------------- #
# Orientation
# --------------------------------------------------------------------------- #

def _placed(tmp_path, n=400, seed=11, s1=51, e1=250, pad=6, probe=40):
    """A reference and the window around one locus inside it."""
    seq = _rnd(n, seed)
    g = _genome(tmp_path, c1=seq)
    return seq, harvest([g], [("c1", s1, e1)], pad=pad, probe=probe)[("c1", s1, e1)]


def test_revcomp_is_its_own_inverse():
    s = _rnd(50, 99)
    assert revcomp(revcomp(s)) == s
    assert revcomp("ACGTN") == "NACGT"


def test_a_forward_record_is_called_plus_and_keeps_its_own_pads(tmp_path):
    seq, w = _placed(tmp_path)
    c = orient(seq[50:250], w)
    assert c.orientation == "+"
    assert c.pad5 == seq[44:50] and c.pad3 == seq[250:256]
    assert c.anchored5 and c.anchored3


def test_a_reverse_complemented_record_is_called_minus_with_swapped_pads(tmp_path):
    """The record's 5' terminus is then at the header's `end`, so the bases
    abutting it are the reverse complement of what follows the header interval."""
    seq, w = _placed(tmp_path)
    c = orient(revcomp(seq[50:250]), w)
    assert c.orientation == "-"
    assert c.pad5 == revcomp(seq[250:256])
    assert c.pad3 == revcomp(seq[44:50])
    assert c.anchored5 and c.anchored3


def test_a_record_shorter_than_its_span_still_anchors_at_both_ends(tmp_path):
    """A nested inner element excised from the middle leaves both termini on the
    header's coordinates, which is why each end is anchored independently."""
    seq = _rnd(800, 12)
    g = _genome(tmp_path, c1=seq)
    w = harvest([g], [("c1", 51, 750)], pad=6, probe=40)[("c1", 51, 750)]
    excised = seq[50:250] + seq[600:750]
    c = orient(excised, w)
    assert c.orientation == "+" and c.anchored5 and c.anchored3


def test_the_same_holds_for_a_reverse_complemented_short_record(tmp_path):
    seq = _rnd(800, 13)
    g = _genome(tmp_path, c1=seq)
    w = harvest([g], [("c1", 51, 750)], pad=6, probe=40)[("c1", 51, 750)]
    c = orient(revcomp(seq[50:250] + seq[600:750]), w)
    assert c.orientation == "-" and c.anchored5 and c.anchored3


def test_an_end_that_lost_bases_is_reported_unanchored(tmp_path):
    seq, w = _placed(tmp_path, seed=14)
    c = orient(seq[50:245], w)          # five bases gone from the 3' end
    assert c.orientation == "+"
    assert c.anchored5 and not c.anchored3


def test_orient_declines_when_the_record_is_not_where_the_header_says(tmp_path):
    _seq, w = _placed(tmp_path, seed=15)
    assert orient(_rnd(200, 16), w) is None


def test_ambiguous_columns_are_not_evidence(tmp_path):
    g = _genome(tmp_path, c1="N" * 400)
    w = harvest([g], [("c1", 51, 250)], pad=6, probe=40)[("c1", 51, 250)]
    assert orient("N" * 200, w) is None


def test_orient_declines_an_empty_window():
    assert orient(_rnd(200, 17), Window()) is None


def test_a_few_substitutions_do_not_break_the_anchor(tmp_path):
    """A handful of real edits should not cost the orientation call; a wholesale
    disagreement should. `MIN_IDENTITY` over `PROBE` bases puts the boundary at
    four substitutions, so three must hold and five must not."""
    seq, w = _placed(tmp_path, seed=18)
    def edited(n):
        rec = list(seq[50:250])
        for i in range(n):
            rec[i * 5] = "ACGT"[("ACGT".index(rec[i * 5]) + 1) % 4]
        return "".join(rec)
    assert orient(edited(3), w).anchored5
    assert orient(edited(4), w).anchored5          # exactly 0.90
    assert not orient(edited(5), w).anchored5
    assert orient(edited(5), w).orientation == "+"  # the 3' end still anchors it


def test_ambiguous_columns_are_dropped_rather_than_scored(tmp_path):
    """`sanitize` turns every ambiguity code into N and a reference gap is N on
    both sides, so scoring N against N as agreement would let two gaps anchor a
    record anywhere. They are excluded from the denominator instead."""
    seq, w = _placed(tmp_path, seed=19)
    rec = list(seq[50:250])
    for i in (3, 11, 27):
        rec[i] = "N"
    c = orient("".join(rec), w)
    assert c.orientation == "+" and c.anchored5 and c.anchored3


# --------------------------------------------------------------------------- #
# Target-site duplications
# --------------------------------------------------------------------------- #

def test_a_tsd_is_the_repeat_immediately_flanking_the_element():
    ctx = "ACGTA" + "TGCCCCCA" + "ACGTA"
    assert find_tsd(ctx, 5, 13) == ("ACGTA", 0, 0)


def test_no_tsd_when_the_two_flanks_differ():
    assert find_tsd("ACGTA" + "TGCCCCCA" + "GGTCC", 5, 13) is None


def test_the_longest_duplication_wins():
    ctx = "TACGTA" + "TGCCCCCA" + "TACGTA"
    assert find_tsd(ctx, 6, 14) == ("TACGTA", 0, 0)


def test_a_boundary_one_base_too_far_out_is_found_at_shift_one():
    """The element's last base actually belongs to the downstream copy, so the
    match appears once the 3' boundary is pulled one base into the element."""
    ctx = "ACGTT" + "TGCCCCCAA" + "CGTTAAAA"
    assert find_tsd(ctx, 5, 14) == ("ACGTT", 0, 1)


def test_the_shift_window_can_be_switched_off():
    ctx = "ACGTT" + "TGCCCCCAA" + "CGTTAAAA"
    assert find_tsd(ctx, 5, 14, shifts=(0,)) is None


def test_homopolymer_flanks_are_refused():
    """A run of one base matches its own reflection at most positions in a
    genome, and would swamp the real signal."""
    assert find_tsd("AAAAAA" + "TGCCCCCA" + "AAAAAA", 6, 14, ks=(6,)) is None


def test_ambiguous_flanks_are_refused():
    assert find_tsd("ANTGA" + "TGCCCCCA" + "ANTGA", 5, 13, ks=(5,)) is None


def test_a_flank_shorter_than_the_kmer_yields_nothing():
    assert find_tsd("AC" + "TGCCCCCA" + "AC", 2, 10, ks=(5,), shifts=(0,)) is None


def test_element_context_places_the_record_between_its_pads():
    c = Context("+", "ACGTA", "TTTGG", anchored5=True, anchored3=True)
    assert element_context("GGGG", c) == ("ACGTA" + "GGGG" + "TTTGG", 5, 9)


def test_element_context_masks_a_pad_that_abuts_nothing():
    """An end that failed to anchor has a pad taken from the wrong place, so it
    is replaced by N and every k-mer touching it is refused on its own merits."""
    c = Context("+", "ACGTA", "TTTGG", anchored5=False, anchored3=True)
    assert element_context("GGGG", c) == ("NNNNN" + "GGGG" + "TTTGG", 5, 9)


# --------------------------------------------------------------------------- #
# Filling the columns
# --------------------------------------------------------------------------- #

def _result(**kw):
    base = {f.name: None for f in fields(Result)}
    base.update(seq_id="x", seq_len=20, status="pass", flank5_len=0, flank3_len=0)
    base.update(kw)
    return Result(**base)


# `shifts=(0,)` throughout: these check that `annotate` routes the right two
# boundaries to `find_tsd`, and the shift search is already covered above.
_EXACT = Options(shifts=(0,))


def test_annotate_fills_orientation_and_both_tsd_columns():
    c = Context("+", "ACGTA", "ACGTA", True, True)
    seq = "TG" + "C" * 16 + "CA"
    r = annotate(_result(seq_len=len(seq)), seq, c, _EXACT)
    assert r.orientation == "+"
    assert r.tsd == "ACGTA" and r.tsd_offset == "0,0" and r.tsd_input == "ACGTA"


def test_the_two_tsd_columns_separate_once_a_flank_is_called():
    """`tsd` describes the boundary Kmer2LTR settled on, `tsd_input` the record
    as it was handed in. They are the same measurement only when no flank was
    called, which is the majority of records but not the interesting ones."""
    c = Context("+", "TTGCA", "TTGCA", True, True)
    seq = "ACGTA" + "GG" + "C" * 10 + "AA" + "ACGTA"
    r = annotate(_result(seq_len=len(seq), flank5_len=5, flank3_len=5),
                 seq, c, _EXACT)
    assert r.tsd == "ACGTA" and r.tsd_offset == "0,0"
    assert r.tsd_input == "TTGCA"


def test_an_absent_duplication_is_a_dot_and_a_missing_one_is_None():
    c = Context("+", "ACGTA", "GGTCC", True, True)
    r = annotate(_result(), "TG" + "C" * 16 + "CA", c, _EXACT)
    assert r.tsd == "." and r.tsd_input == "."
    assert r.tsd_offset is None                 # rendered NA: nothing to offset


def test_a_row_with_no_pair_still_reports_what_the_record_itself_shows():
    """`orientation` and `tsd_input` are properties of the record, not of a pair
    that was never found -- and on such a row they are the only thing that can
    still say whether the annotator was pointing at a real insertion."""
    c = Context("-", "ACGTA", "ACGTA", True, True)
    r = annotate(_result(status="no_pair", flank5_len=None, flank3_len=None),
                 "TG" + "C" * 16 + "CA", c, _EXACT)
    assert r.orientation == "-" and r.tsd_input == "ACGTA"
    assert r.tsd is None and r.tsd_offset is None


def test_annotate_leaves_every_other_field_alone():
    c = Context("+", "ACGTA", "ACGTA", True, True)
    before = _result(seq_len=20, k2p=0.125, cigar="20=", motif="tg...ca")
    after = annotate(before, "TG" + "C" * 16 + "CA", c, _EXACT)
    assert (after.k2p, after.cigar, after.motif) == (0.125, "20=", "tg...ca")


def test_credit_is_paid_only_for_a_duplication_at_the_untrimmed_termini():
    opt = Options(anchor=12.0, shifts=(0,))
    yes = Context("+", "ACGTA", "ACGTA", True, True)
    no = Context("+", "ACGTA", "GGTCC", True, True)
    seq = "TG" + "C" * 16 + "CA"
    assert credit(seq, yes, opt) == 12.0
    assert credit(seq, no, opt) == 0.0
    assert credit(seq, yes, Options()) == 0.0       # anchor 0 is off
    assert credit(seq, None, opt) == 0.0


def test_credit_ignores_the_shift_window():
    """A duplication that only appears once a boundary is moved says the
    boundary is wrong, not that it is right, so it must not buy the record the
    benefit of the doubt."""
    opt = Options(anchor=12.0, shifts=(0, 1, -1))
    c = Context("+", "ACGTT", "CGTTAAAA", True, True)
    seq = "TGCCCCCAA"
    assert annotate(_result(seq_len=len(seq)), seq, c, opt).tsd_input == "ACGTT"
    assert credit(seq, c, opt) == 0.0


def test_an_unanchored_end_reports_NA_not_absence():
    """The pad on an end that did not anchor abuts nothing, so there is nothing
    to have found. `.` would claim a measurement nobody made."""
    c = Context("+", "ACGTA", "ACGTA", anchored5=False, anchored3=True)
    r = annotate(_result(), "TG" + "C" * 16 + "CA", c, _EXACT)
    assert r.orientation == "+"
    assert r.tsd_input is None and r.tsd is None


def test_an_element_at_a_contig_edge_reports_NA_not_absence():
    """`harvest` clips at the contig, so there is no pad to read a k-mer from."""
    c = Context("+", "", "ACGTA", anchored5=True, anchored3=True)
    r = annotate(_result(), "TG" + "C" * 16 + "CA", c, _EXACT)
    assert r.tsd_input is None and r.tsd is None


def test_a_reference_gap_beside_the_element_reports_NA_not_absence():
    c = Context("+", "ACNTA", "ACGTA", anchored5=True, anchored3=True)
    r = annotate(_result(), "TG" + "C" * 16 + "CA", c, _EXACT)
    assert r.tsd_input is None


def test_a_readable_pad_with_no_duplication_reports_absence():
    c = Context("+", "ACGTA", "GGTCC", anchored5=True, anchored3=True)
    r = annotate(_result(), "TG" + "C" * 16 + "CA", c, _EXACT)
    assert r.tsd_input == "." and r.tsd == "."


def test_a_called_boundary_far_from_the_pad_is_readable_regardless(tmp_path):
    """With a flank called, the k-mers come from the record's own bases, so an
    unanchored pad says nothing about whether the search could run there."""
    c = Context("+", "ACGTA", "ACGTA", anchored5=False, anchored3=False)
    seq = "AACCT" + "TG" + "C" * 10 + "CA" + "AACCT"
    r = annotate(_result(seq_len=len(seq), flank5_len=5, flank3_len=5),
                 seq, c, _EXACT)
    assert r.tsd == "AACCT"        # read entirely from the record
    assert r.tsd_input is None     # but the record's own termini are unreadable


def test_a_zero_shift_hit_beats_a_longer_one_at_a_shifted_boundary():
    """Both fit here: a 6 bp duplication if the 3' boundary is a base out, and a
    5 bp one if it is exactly right. Reporting the second makes the smaller
    claim, and the sweep found the choice never changes what is detected."""
    ctx = "TACGTA" + "TGCCCCCT" + "ACGTAGGG"
    assert find_tsd(ctx, 6, 14) == ("ACGTA", 0, 0)



def test_harvest_reads_an_unwrapped_reference(tmp_path):
    """One line per contig is what `seqkit seq -w 0` and several assemblers
    emit. Reading it a line at a time would cost twice the contig."""
    seq = _rnd(500, 21)
    p = tmp_path / "flat.fa"
    p.write_text(">c1 desc\n" + seq + "\n")
    w = harvest([p], [("c1", 101, 400)], pad=5, probe=8)[("c1", 101, 400)]
    assert w == Window(seq[95:100], seq[100:108], seq[392:400], seq[400:405])


def test_harvest_survives_a_reference_with_no_trailing_newline(tmp_path):
    seq = _rnd(200, 22)
    p = tmp_path / "x.fa"
    p.write_text(">c1\n" + seq)
    w = harvest([p], [("c1", 51, 150)], pad=4, probe=4)[("c1", 51, 150)]
    assert w == Window(seq[46:50], seq[50:54], seq[146:150], seq[150:154])


def test_whitespace_inside_a_sequence_line_does_not_shift_the_windows(tmp_path):
    """It would shift the base offset every cut is measured against, so the
    harvest would come back from the wrong place rather than fail."""
    seq = _rnd(120, 23)
    clean = tmp_path / "clean.fa"
    _genome(tmp_path, "clean.fa", wrap=10, c1=seq)
    spaced = tmp_path / "spaced.fa"
    with open(spaced, "w") as fh:
        fh.write(">c1\n")
        for i in range(0, len(seq), 20):
            fh.write(seq[i:i + 10] + " " + seq[i + 10:i + 20] + "\n")
    loci = [("c1", 31, 90)]
    assert harvest([spaced], loci, pad=6, probe=8) == harvest([clean], loci, pad=6, probe=8)


def test_a_repeated_contig_name_is_reported(tmp_path, capsys):
    """Reachable from `--genome primary.fa alts.fa`, or a badly concatenated
    multi-genome reference. The last copy wins, so say so."""
    a = _genome(tmp_path, "a.fa", c1="ACGT" * 20)
    b = _genome(tmp_path, "b.fa", c1="TTTT" * 20)
    w = harvest([a, b], [("c1", 21, 40)], pad=4, probe=4)[("c1", 21, 40)]
    assert w.head == "TTTT"
    assert "appears more than once" in capsys.readouterr().err


def test_shifts_are_tried_in_order_of_total_displacement():
    """`(1, 1)` moves two boundaries and `(-1, 0)` moves one, so a lexicographic
    pass over the shift set would report the larger claim first."""
    from kmer2ltr.genome import _shift_pairs
    got = _shift_pairs((0, 1, -1))
    assert got[0] == (0, 0)
    assert [abs(a) + abs(b) for a, b in got] == sorted(abs(a) + abs(b) for a, b in got)
    assert got.index((-1, 0)) < got.index((1, 1))


def test_a_zero_based_bedtools_header_still_anchors_both_ends(tmp_path):
    """`bedtools getfasta -name` writes the BED interval verbatim, so its
    `chr:1000-2000` means the 1-based `chr:1001-2000`. Read as 1-based, the
    `end`-anchored terminus still matches and the record is still oriented --
    the only symptom is an empty TSD column for the whole file."""
    seq = _rnd(400, 31)
    g = _genome(tmp_path, c1=seq)
    # the element is 1-based 51-250; a BED header names it 50-250
    w = harvest([g], [("c1", 50, 250)], pad=6, probe=40)[("c1", 50, 250)]
    c = orient(seq[50:250], w)
    assert c.orientation == "+"
    assert c.anchored5 and c.anchored3
    assert c.pad5 == seq[44:50] and c.pad3 == seq[250:256]


def test_a_zero_based_header_works_for_a_reversed_record_too(tmp_path):
    seq = _rnd(400, 32)
    g = _genome(tmp_path, c1=seq)
    w = harvest([g], [("c1", 50, 250)], pad=6, probe=40)[("c1", 50, 250)]
    c = orient(revcomp(seq[50:250]), w)
    assert c.orientation == "-" and c.anchored5 and c.anchored3
    assert c.pad5 == revcomp(seq[250:256]) and c.pad3 == revcomp(seq[44:50])


def test_a_one_based_header_is_not_dragged_off_by_the_fallback(tmp_path):
    seq = _rnd(400, 33)
    g = _genome(tmp_path, c1=seq)
    w = harvest([g], [("c1", 51, 250)], pad=6, probe=40)[("c1", 51, 250)]
    c = orient(seq[50:250], w)
    assert c.pad5 == seq[44:50]          # not seq[45:51]
