import random, subprocess, sys
import pytest
from pathlib import Path
from ltrk2p.cli import main
from ltrk2p.runner import COLUMNS

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
    r = subprocess.run([PY, "-m", "ltrk2p", str(inp)], capture_output=True, text=True)
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
    from ltrk2p import cluster
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
    from ltrk2p import cluster
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
    from ltrk2p import cluster
    inp = _fa_flanked(tmp_path, n=3)
    monkeypatch.setattr(cluster, "available", lambda: True)
    monkeypatch.setattr(cluster, "cluster", lambda *a, **k: [])
    assert main([str(inp), "-o", str(tmp_path / "o.tsv"), "--ltr-cluster"]) == 1


def test_the_internal_fasta_survives_an_incomplete_clustering(tmp_path, monkeypatch):
    """It is the expensive half of the run. Dropping it when some identities
    failed makes those tables unrecoverable without re-aligning everything."""
    from ltrk2p import cluster
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
