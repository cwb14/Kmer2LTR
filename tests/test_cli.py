import random, subprocess, sys
import pytest
from pathlib import Path
from kmer2ltr.cli import main
from kmer2ltr.runner import COLUMNS

PY = sys.executable

def _rnd(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))

def _fa(tmp_path, n=5):
    p = tmp_path / "in.fa"
    recs = []
    for i in range(n):
        ltr = _rnd(300, i)
        recs.append(f">r{i}\n{ltr + _rnd(600, i+50) + ltr}\n")
    p.write_text("".join(recs))
    return p

def test_writes_tsv_to_output_file(tmp_path):
    inp = _fa(tmp_path)
    out = tmp_path / "out.tsv"
    assert main([str(inp), "-o", str(out)]) == 0
    lines = out.read_text().rstrip("\n").split("\n")
    assert len(lines) == 6           # header + 5
    assert lines[0].startswith("seq_id\t")

def test_row_count_equals_record_count(tmp_path):
    inp = _fa(tmp_path, n=12)
    out = tmp_path / "out.tsv"
    main([str(inp), "-o", str(out)])
    n_in = sum(1 for l in inp.read_text().split("\n") if l.startswith(">"))
    n_out = len(out.read_text().rstrip("\n").split("\n")) - 1
    assert n_in == n_out

def test_cs_flag_changes_the_cigar_column(tmp_path):
    inp = _fa(tmp_path, n=2)
    a, b = tmp_path / "a.tsv", tmp_path / "b.tsv"
    main([str(inp), "-o", str(a)])
    main([str(inp), "--cs", "-o", str(b)])
    # Indexed by name, not by -1: `cigar` is no longer the last column.
    col = COLUMNS.index("cigar")
    cig_a = a.read_text().rstrip("\n").split("\n")[1].split("\t")[col]
    cig_b = b.read_text().rstrip("\n").split("\n")[1].split("\t")[col]
    assert cig_a != cig_b

def test_resume_appends_only_missing_rows(tmp_path):
    inp = _fa(tmp_path, n=10)
    out = tmp_path / "out.tsv"
    main([str(inp), "-o", str(out)])
    full = out.read_text()
    # truncate to header + 4 rows, then resume
    out.write_text("\n".join(full.rstrip("\n").split("\n")[:5]) + "\n")
    main([str(inp), "-o", str(out), "--resume"])
    assert out.read_text() == full

def test_flank_bits_flag_actually_reaches_classify(tmp_path):
    """Regression: tuning flags must be threaded, not assigned to module
    globals -- those are already bound as default args, so assignment is a
    silent no-op and the flag would appear to work while doing nothing."""
    ltr = _rnd(400, 40)
    # 40 bp of flank: detected at a normal threshold, absorbed at a huge one
    seq = _rnd(40, 41) + ltr + _rnd(900, 42) + ltr
    inp = tmp_path / "f.fa"
    inp.write_text(f">e\n{seq}\n")
    lo, hi = tmp_path / "lo.tsv", tmp_path / "hi.tsv"
    main([str(inp), "-o", str(lo), "--flank-bits", "5"])
    main([str(inp), "-o", str(hi), "--flank-bits", "10000"])
    f_lo = lo.read_text().rstrip("\n").split("\n")[1].split("\t")[9]
    f_hi = hi.read_text().rstrip("\n").split("\n")[1].split("\t")[9]
    assert f_lo != f_hi, f"--flank-bits had no effect (flank5_len {f_lo} both ways)"


def test_missing_input_fails_fast_with_clear_message(tmp_path, capsys):
    rc = main([str(tmp_path / "nope.fa"), "-o", str(tmp_path / "o.tsv")])
    assert rc != 0
    assert "not found" in capsys.readouterr().err.lower()

def test_console_entry_point_installed(tmp_path):
    inp = _fa(tmp_path, n=2)
    r = subprocess.run([PY, "-m", "kmer2ltr", str(inp)], capture_output=True, text=True)
    assert r.returncode == 0
    assert r.stdout.startswith("seq_id\t")

def test_extreme_flank_bits_does_not_crash(tmp_path):
    """--flank-bits >= ~2000 forced boundaries across kb of unrelated flank,
    producing a bitscore below -1024 and an OverflowError that aborted the run
    with a truncated output file. Needs a LONG record; a ~1kb fixture cannot
    reproduce it."""
    ltr = _rnd(400, 3)
    seq = _rnd(3000, 4) + ltr + _rnd(4000, 5) + ltr + _rnd(3000, 6)
    inp = tmp_path / "long.fa"
    inp.write_text(f">long\n{seq}\n")
    out = tmp_path / "o.tsv"
    assert main([str(inp), "-o", str(out), "--flank-bits", "5000"]) == 0
    rows = out.read_text().rstrip("\n").split("\n")
    assert len(rows) == 2          # header + exactly one record, not truncated


# --------------------------------------------------------------------------- #
# Mutation rate and the extra outputs
# --------------------------------------------------------------------------- #

def _fa_flanked(tmp_path, n=4, flank=45, name="f.fa"):
    """Elements with a genomic-coordinate header and known flank on both ends."""
    p = tmp_path / name
    recs = []
    for i in range(n):
        ltr = "TG" + _rnd(298, i) + "CA"
        seq = (_rnd(flank, 200 + i) + ltr + _rnd(700, 300 + i) + ltr
               + _rnd(flank, 400 + i))
        recs.append(f">chr1:{1000 + i * 9000}-{1000 + i * 9000 + len(seq) - 1}#LTR/Gypsy\n{seq}\n")
    p.write_text("".join(recs))
    return p


def _rows(path):
    lines = path.read_text().rstrip("\n").split("\n")
    return [dict(zip(COLUMNS, ln.split("\t"))) for ln in lines[1:]]


def test_mutation_rate_fills_the_k2p_time_column(tmp_path):
    inp = _fa_flanked(tmp_path, n=2)
    plain, dated = tmp_path / "a.tsv", tmp_path / "b.tsv"
    assert main([str(inp), "-o", str(plain)]) == 0
    assert main([str(inp), "-o", str(dated), "-u", "7e-9"]) == 0
    assert all(r["k2p_time"] == "NA" for r in _rows(plain))
    assert all(r["k2p_time"] != "NA" for r in _rows(dated) if r["status"] == "pass")


def test_a_nonpositive_mutation_rate_is_rejected(tmp_path, capsys):
    inp = _fa_flanked(tmp_path, n=1)
    assert main([str(inp), "-o", str(tmp_path / "o.tsv"), "-u", "0"]) == 2
    assert "mutation-rate" in capsys.readouterr().err


def test_extra_outputs_require_an_output_file(tmp_path, capsys):
    """They are named after -o, so there is nothing to name them after when the
    table goes to stdout."""
    inp = _fa_flanked(tmp_path, n=1)
    assert main([str(inp), "--trim-flanks"]) == 2
    assert "-o is required" in capsys.readouterr().err


def test_extra_outputs_refuse_to_resume(tmp_path, capsys):
    """The TSV resumes on its data-line count, but the auxiliary FASTAs hold
    only the passing subset -- that count cannot say where they stopped, so
    appending would duplicate or drop records in them."""
    inp = _fa_flanked(tmp_path, n=2)
    out = tmp_path / "o.tsv"
    assert main([str(inp), "-o", str(out)]) == 0
    assert main([str(inp), "-o", str(out), "--resume", "--trim-flanks"]) == 2
    assert "--resume cannot be combined" in capsys.readouterr().err


def test_min_seq_id_is_rejected_without_a_clustering_mode(tmp_path, capsys):
    inp = _fa_flanked(tmp_path, n=1)
    rc = main([str(inp), "-o", str(tmp_path / "o.tsv"), "--min-seq-id", "0.8"])
    assert rc == 2 and "no effect" in capsys.readouterr().err


def test_trim_flanks_writes_corrected_elements_and_corrected_headers(tmp_path):
    inp = _fa_flanked(tmp_path, n=3, flank=45)
    out = tmp_path / "o.tsv"
    assert main([str(inp), "-o", str(out), "--trim-flanks"]) == 0
    trimmed = (tmp_path / "o.trimmed.fa").read_text()
    src = {}
    for block in inp.read_text().split(">")[1:]:
        head, _, body = block.partition("\n")
        src[head] = body.replace("\n", "")
    written = [b for b in trimmed.split(">")[1:]]
    assert len(written) == len(_rows(out))
    for row, block in zip(_rows(out), written):
        head, _, body = block.partition("\n")
        chrom, span = row["seq_id"].split("#")[0].rsplit(":", 1)
        start, end = (int(v) for v in span.split("-"))
        expect = f"{chrom}:{start + int(row['flank5_len'])}-{end - int(row['flank3_len'])}#LTR/Gypsy"
        assert head == expect
        assert body.replace("\n", "") == src[row["seq_id"]][
            int(row["ltr5_start"]) - 1:int(row["ltr3_end"])]


def test_perfect_ltr_rt_writes_one_file_per_requested_mode(tmp_path):
    inp = _fa_flanked(tmp_path, n=2)
    out = tmp_path / "o.tsv"
    assert main([str(inp), "-o", str(out), "--perfect-ltr-rt", "5p", "consensus"]) == 0
    for mode in ("5p", "consensus"):
        text = (tmp_path / f"o.perfect_{mode}.fa").read_text()
        assert text.count(">") == 2 and "~LTRlen:" in text
    assert not (tmp_path / "o.perfect_3p.fa").exists()


def test_an_input_with_no_passing_element_leaves_no_empty_fasta(tmp_path):
    """A zero-byte FASTA reads as "an empty result" rather than "nothing
    qualified", and is the kind of leftover that outlives the run."""
    inp = tmp_path / "junk.fa"
    inp.write_text(f">a\n{_rnd(2000, 77)}\n")
    out = tmp_path / "o.tsv"
    assert main([str(inp), "-o", str(out), "--trim-flanks"]) == 0
    assert not (tmp_path / "o.trimmed.fa").exists()
    assert sorted(p.name for p in tmp_path.iterdir()) == ["junk.fa", "o.tsv"]


def test_plot_writes_a_density_pdf(tmp_path):
    pytest.importorskip("matplotlib")
    inp = _fa_flanked(tmp_path, n=6)
    out = tmp_path / "o.tsv"
    assert main([str(inp), "-o", str(out), "-u", "7e-9", "--plot"]) == 0
    assert (tmp_path / "o.density.pdf").read_bytes()[:5] == b"%PDF-"


def test_ltr_cluster_keeps_the_consensus_and_internal_cluster_drops_its_fasta(tmp_path):
    """The cluster table names sequences that live only in the consensus FASTA,
    so that file is a result. The internal FASTA is scratch for a second
    opinion on the same grouping, so it is not."""
    from kmer2ltr import cluster
    if not cluster.available():
        pytest.skip("mmseqs not on PATH")
    inp = _fa_flanked(tmp_path, n=5)
    out = tmp_path / "o.tsv"
    rc = main([str(inp), "-o", str(out), "--ltr-cluster", "--internal-cluster",
               "--min-seq-id", "0.8"])
    assert rc == 0
    assert (tmp_path / "o.consensus.fa").exists()
    assert not (tmp_path / "o.internal.fa").exists()
    assert (tmp_path / "o.consensus_id0.80_cluster.tsv").exists()
    assert (tmp_path / "o.internal_id0.80_cluster.tsv").exists()
    leftovers = [p.name for p in tmp_path.iterdir()
                 if p.is_dir() or p.name.endswith(("_all_seqs.fasta", "_rep_seq.fasta"))]
    assert leftovers == [], f"mmseqs scratch survived: {leftovers}"


def test_clustering_without_mmseqs_fails_before_doing_any_work(tmp_path, capsys, monkeypatch):
    from kmer2ltr import cluster
    monkeypatch.setattr(cluster, "available", lambda: False)
    inp = _fa_flanked(tmp_path, n=1)
    out = tmp_path / "o.tsv"
    assert main([str(inp), "-o", str(out), "--ltr-cluster"]) == 2
    assert "mmseqs" in capsys.readouterr().err
    assert not out.exists(), "failed fast, so nothing should have been written"


# --------------------------------------------------------------------------- #
# Regressions
# --------------------------------------------------------------------------- #

def test_an_extra_output_is_refused_when_it_would_land_on_the_input(tmp_path, capsys):
    """Every output is opened for WRITING before the input is read, so a derived
    path that IS the input truncates the file about to be parsed -- and the
    emptied stream is then deleted as "received no records". Chaining a
    --trim-flanks output back in is the obvious way to land here."""
    inp = _fa_flanked(tmp_path, n=2, name="elements.trimmed.fa")
    before = inp.read_bytes()
    rc = main([str(inp), "-o", str(tmp_path / "elements.tsv"), "--trim-flanks"])
    assert rc == 2
    assert "would be destroyed" in capsys.readouterr().err
    assert inp.read_bytes() == before, "the input was modified"


def test_output_is_refused_when_it_is_the_input(tmp_path, capsys):
    inp = _fa_flanked(tmp_path, n=2)
    before = inp.read_bytes()
    assert main([str(inp), "-o", str(inp)]) == 2
    assert "would be destroyed" in capsys.readouterr().err
    assert inp.read_bytes() == before


def test_a_directory_input_fails_cleanly(tmp_path, capsys):
    assert main([str(tmp_path), "-o", str(tmp_path / "o.tsv")]) == 2
    assert "is a directory" in capsys.readouterr().err


@pytest.mark.parametrize("out", ["", "."])
def test_an_unusable_output_path_fails_cleanly(tmp_path, capsys, out):
    """`Path('.').with_suffix('')` raises; it must not reach the user as a
    traceback from the middle of a run."""
    inp = _fa_flanked(tmp_path, n=1)
    assert main([str(inp), "-o", out]) == 2
    assert "not a usable path" in capsys.readouterr().err


def test_min_seq_id_finer_than_the_naming_is_refused(tmp_path, capsys):
    """Identities are used to two decimals in both the mmseqs argument and the
    file name, so 0.999 would run at 1.00 -- one cluster per sequence -- in a
    file named for the value it was not run at."""
    inp = _fa_flanked(tmp_path, n=1)
    rc = main([str(inp), "-o", str(tmp_path / "o.tsv"), "--ltr-cluster",
               "--min-seq-id", "0.999"])
    assert rc == 2
    assert "two decimal places" in capsys.readouterr().err


def test_resume_onto_a_missing_output_still_writes_the_header(tmp_path):
    """`--resume` in a restartable job wrapper hits this on the FIRST run. A
    headerless TSV costs every downstream csv.DictReader its first data row."""
    inp = _fa_flanked(tmp_path, n=4)
    out = tmp_path / "fresh.tsv"
    assert main([str(inp), "-o", str(out), "--resume"]) == 0
    lines = out.read_text().rstrip("\n").split("\n")
    assert lines[0].split("\t") == COLUMNS
    assert len(lines) == 5, "header + 4 records"
    # A second resume must now find 4 done and add nothing.
    assert main([str(inp), "-o", str(out), "--resume"]) == 0
    lines = out.read_text().rstrip("\n").split("\n")
    assert len(lines) == 5
    ids = [ln.split("\t")[0] for ln in lines[1:]]
    assert len(set(ids)) == 4, f"records duplicated: {ids}"


def test_resume_discards_a_partial_final_line(tmp_path):
    """A run killed mid-write leaves one. Counting it as complete skips the
    record it belongs to and then appends onto the fragment, producing one lost
    record and one row holding two."""
    inp = _fa_flanked(tmp_path, n=5)
    full = tmp_path / "full.tsv"
    assert main([str(inp), "-o", str(full)]) == 0
    complete = full.read_text()
    n_full = len(complete.rstrip("\n").split("\n"))

    trunc = tmp_path / "trunc.tsv"
    cut = complete.index("\n", complete.index("\n") + 1) + 40   # mid row 2
    trunc.write_text(complete[:cut])
    assert main([str(inp), "-o", str(trunc), "--resume"]) == 0

    lines = trunc.read_text().rstrip("\n").split("\n")
    assert len(lines) == n_full, "a record was lost or duplicated"
    assert {len(ln.split("\t")) for ln in lines} == {len(COLUMNS)}
    assert len({ln.split("\t")[0] for ln in lines[1:]}) == 5


def test_an_incomplete_clustering_is_visible_in_the_exit_code(tmp_path, monkeypatch):
    """Otherwise a shell pipeline carries on into tables that are not there."""
    from kmer2ltr import cluster
    inp = _fa_flanked(tmp_path, n=3)
    monkeypatch.setattr(cluster, "available", lambda: True)
    monkeypatch.setattr(cluster, "cluster", lambda *a, **k: [])
    assert main([str(inp), "-o", str(tmp_path / "o.tsv"), "--ltr-cluster"]) == 1


def test_the_internal_fasta_survives_an_incomplete_clustering(tmp_path, monkeypatch):
    """It is the expensive half of the run. Dropping it when some identities
    failed makes those tables unrecoverable without re-aligning everything."""
    from kmer2ltr import cluster
    inp = _fa_flanked(tmp_path, n=3)
    monkeypatch.setattr(cluster, "available", lambda: True)

    def partial(fasta, *a, **k):
        out = Path(f"{Path(fasta).with_suffix('')}_id0.70_cluster.tsv")
        out.write_text("a\ta\n")
        return [out]
    monkeypatch.setattr(cluster, "cluster", partial)
    rc = main([str(inp), "-o", str(tmp_path / "o.tsv"), "--internal-cluster"])
    assert rc == 1
    assert (tmp_path / "o.internal.fa").exists(), "scratch dropped mid-sweep"


# --------------------------------------------------------------------------- #
# --genome
# --------------------------------------------------------------------------- #

_COMP = str.maketrans("ACGT", "TGCA")


def _rc(s):
    return s.translate(_COMP)[::-1]


def _placed(tmp_path, tsd="ACGTA", pad=(0, 0), reverse=False, n=3):
    """A reference, and an element FASTA cut out of it at known loci.

    Each planted record is `pad[0]` bases of genome, an element, and `pad[1]`
    more, with a real target-site duplication wrapped around the whole record --
    so the duplication is readable only from the reference, which is the point
    of the flag. Pass a `(left, right)` pair for `tsd` to plant flanks that do
    not match. `reverse` stores the record on the other strand while leaving
    forward coordinates in its header, as annotation pipelines do.
    """
    left, right = (tsd, tsd) if isinstance(tsd, str) else tsd
    parts, loci = [], []
    pos = 0
    for i in range(n):
        lead = _rnd(500, 900 + i)
        ltr = _rnd(300, 100 + i)
        rec = (_rnd(pad[0], 300 + i) + ltr + _rnd(600, 200 + i) + ltr
               + _rnd(pad[1], 400 + i))
        parts.append(lead + left + rec + right)
        start = pos + len(lead) + len(left) + 1       # 1-based, record only
        loci.append((start, start + len(rec) - 1))
        pos += len(lead) + len(left) + len(right) + len(rec)
    ref_seq = "".join(parts)

    g = tmp_path / "ref.fa"
    with open(g, "w") as fh:
        fh.write(">c1 test contig\n")
        for i in range(0, len(ref_seq), 60):
            fh.write(ref_seq[i:i + 60] + "\n")
    fa = tmp_path / "elements.fa"
    with open(fa, "w") as fh:
        for s1, e1 in loci:
            piece = ref_seq[s1 - 1:e1]
            fh.write(f">c1:{s1}-{e1}#LTR/Copia\n"
                     f"{_rc(piece) if reverse else piece}\n")
    return fa, g


def _col(path, name):
    lines = path.read_text().rstrip("\n").split("\n")
    i = lines[0].split("\t").index(name)
    return [l.split("\t")[i] for l in lines[1:]]


def test_genome_adds_the_columns_without_moving_any_other_field(tmp_path):
    """The reference adds information; at the default --tsd-anchor it must not
    move a single boundary."""
    fa, g = _placed(tmp_path)
    a, b = tmp_path / "a.tsv", tmp_path / "b.tsv"
    assert main([str(fa), "-o", str(a)]) == 0
    assert main([str(fa), "-o", str(b), "--genome", str(g)]) == 0
    la = a.read_text().rstrip("\n").split("\n")
    lb = b.read_text().rstrip("\n").split("\n")
    assert len(la) == len(lb)
    for x, y in zip(la, lb):
        assert x.split("\t")[:25] == y.split("\t")[:25]
    assert _col(a, "tsd") == ["NA"] * 3
    assert all("ACGTA" in v for v in _col(b, "tsd"))
    assert _col(b, "tsd") == _col(b, "tsd_input")     # no flank was called
    assert _col(b, "orientation") == ["+"] * 3


def test_a_reverse_complemented_record_is_reported_as_such(tmp_path):
    fa, g = _placed(tmp_path, reverse=True)
    out = tmp_path / "o.tsv"
    assert main([str(fa), "-o", str(out), "--genome", str(g)]) == 0
    assert _col(out, "orientation") == ["-"] * 3
    assert all(_rc("ACGTA") in v for v in _col(out, "tsd"))


def test_the_two_tsd_columns_separate_when_a_flank_is_called(tmp_path):
    fa, g = _placed(tmp_path, pad=(40, 40))
    out = tmp_path / "o.tsv"
    assert main([str(fa), "-o", str(out), "--genome", str(g)]) == 0
    assert all(int(v) > 0 for v in _col(out, "flank5_len"))
    assert all("ACGTA" in v for v in _col(out, "tsd_input"))   # as supplied
    assert _col(out, "tsd") == ["."] * 3                       # not where it cut


def test_trim_flanks_takes_the_five_prime_trim_off_end_when_reversed(tmp_path):
    """The header carries forward coordinates; a reversed record's 5' terminus
    is at `end`, so an asymmetric trim must move the far coordinate."""
    from kmer2ltr.extras import shift_locus
    fa, g = _placed(tmp_path, pad=(60, 15), reverse=True, n=1)
    out = tmp_path / "o.tsv"
    assert main([str(fa), "-o", str(out), "--genome", str(g), "--trim-flanks"]) == 0
    header = [l[1:] for l in (tmp_path / "o.trimmed.fa").read_text().split("\n")
              if l.startswith(">")][0]
    orig = [l[1:] for l in fa.read_text().split("\n") if l.startswith(">")][0]
    f5, f3 = int(_col(out, "flank5_len")[0]), int(_col(out, "flank3_len")[0])
    assert f5 != f3 and f5 and f3
    assert header == shift_locus(orig, f5, f3, "-")
    assert header != shift_locus(orig, f5, f3, "+")


def test_records_with_no_locus_in_the_header_report_NA(tmp_path):
    fa, g = _placed(tmp_path, n=1)
    plain = tmp_path / "plain.fa"
    plain.write_text(fa.read_text().replace(">c1:", ">elem_"))
    out = tmp_path / "o.tsv"
    assert main([str(plain), "-o", str(out), "--genome", str(g)]) == 0
    assert _col(out, "orientation") == ["NA"]
    assert _col(out, "tsd") == ["NA"] and _col(out, "tsd_input") == ["NA"]


def test_a_locus_on_a_contig_the_reference_lacks_reports_NA(tmp_path):
    fa, g = _placed(tmp_path, n=1)
    other = tmp_path / "other.fa"
    other.write_text(fa.read_text().replace(">c1:", ">nope:"))
    out = tmp_path / "o.tsv"
    assert main([str(other), "-o", str(out), "--genome", str(g)]) == 0
    assert _col(out, "orientation") == ["NA"]


def test_tsd_anchor_needs_a_genome(tmp_path, capsys):
    fa, _g = _placed(tmp_path, n=1)
    assert main([str(fa), "-o", str(tmp_path / "o.tsv"), "--tsd-anchor", "8"]) == 2
    assert "--tsd-anchor needs --genome" in capsys.readouterr().err


def test_a_negative_tsd_anchor_is_refused(tmp_path, capsys):
    fa, g = _placed(tmp_path, n=1)
    assert main([str(fa), "-o", str(tmp_path / "o.tsv"),
                 "--genome", str(g), "--tsd-anchor", "-1"]) == 2
    assert "--tsd-anchor must be >= 0" in capsys.readouterr().err


def test_a_missing_genome_is_refused(tmp_path, capsys):
    fa, _g = _placed(tmp_path, n=1)
    assert main([str(fa), "-o", str(tmp_path / "o.tsv"),
                 "--genome", str(tmp_path / "nope.fa")]) == 2
    assert "--genome file not found" in capsys.readouterr().err


def test_an_output_that_would_land_on_the_reference_is_refused(tmp_path, capsys):
    """Outputs are opened for writing before the reference is read, so naming
    one as the other would destroy it."""
    fa, g = _placed(tmp_path, n=1)
    assert main([str(fa), "-o", str(g), "--genome", str(g)]) == 2
    assert "would be written to an input file" in capsys.readouterr().err


def test_the_genome_columns_are_thread_count_independent(tmp_path):
    fa, g = _placed(tmp_path, n=6)
    a, b = tmp_path / "a.tsv", tmp_path / "b.tsv"
    assert main([str(fa), "-o", str(a), "--genome", str(g), "-t", "1"]) == 0
    assert main([str(fa), "-o", str(b), "--genome", str(g), "-t", "3"]) == 0
    assert a.read_text() == b.read_text()


def test_tsd_anchor_suppresses_a_flank_the_duplication_argues_against(tmp_path):
    fa, g = _placed(tmp_path, pad=(40, 40), n=2)
    off, on = tmp_path / "off.tsv", tmp_path / "on.tsv"
    assert main([str(fa), "-o", str(off), "--genome", str(g)]) == 0
    assert main([str(fa), "-o", str(on), "--genome", str(g),
                 "--tsd-anchor", "1e6"]) == 0
    assert all(int(v) > 0 for v in _col(off, "flank5_len"))
    assert all(int(v) == 0 for v in _col(on, "flank5_len"))
    assert all(int(v) == 0 for v in _col(on, "flank3_len"))


def test_tsd_anchor_does_nothing_where_no_duplication_backs_the_record(tmp_path):
    """The credit is paid for evidence, not for the flag being on."""
    fa, g = _placed(tmp_path, tsd=("ACGTA", "GGTCC"), pad=(40, 40), n=2)
    off, on = tmp_path / "off.tsv", tmp_path / "on.tsv"
    assert main([str(fa), "-o", str(off), "--genome", str(g)]) == 0
    assert main([str(fa), "-o", str(on), "--genome", str(g),
                 "--tsd-anchor", "1e6"]) == 0
    assert _col(off, "tsd_input") == ["."] * 2
    assert _col(off, "flank5_len") == _col(on, "flank5_len")
    assert _col(off, "flank3_len") == _col(on, "flank3_len")


def test_resume_reproduces_an_uninterrupted_run_with_a_genome(tmp_path):
    """The reference is harvested for every locus before any record is skipped,
    so a resumed run sees the same context an uninterrupted one did."""
    fa, g = _placed(tmp_path, pad=(20, 20), n=5)
    whole, part = tmp_path / "whole.tsv", tmp_path / "part.tsv"
    assert main([str(fa), "-o", str(whole), "--genome", str(g)]) == 0
    lines = whole.read_text().split("\n")
    part.write_text("\n".join(lines[:3]) + "\n")          # header + 2 records
    assert main([str(fa), "-o", str(part), "--genome", str(g), "--resume"]) == 0
    assert part.read_text() == whole.read_text()
