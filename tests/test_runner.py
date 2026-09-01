import io, random, gzip
import pytest
from ltrk2p.runner import COLUMNS, format_row, run, count_data_lines
from ltrk2p.align import classify

def _rnd(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))

def _fasta(tmp_path, records, name="in.fa"):
    p = tmp_path / name
    p.write_text("".join(f">{i}\n{s}\n" for i, s in records))
    return p

def test_columns_match_spec_order_and_count():
    assert len(COLUMNS) == 25
    assert COLUMNS[0] == "seq_id" and COLUMNS[2] == "status"
    assert COLUMNS[3:7] == ["ltr5_start", "ltr5_end", "ltr3_start", "ltr3_end"]
    assert COLUMNS[18] == "k2p" and COLUMNS[22] == "cigar"


def test_motif_and_k2p_time_were_appended_after_the_columns_that_existed_before():
    """`motif` and `k2p_time` sit AFTER `cigar` so that every column predating
    them keeps its index -- the README's own `cut -f1,4-7,19,23` still selects
    the same six fields, and no downstream script shifts by two."""
    assert COLUMNS[23:] == ["motif", "k2p_time"]

def test_format_row_renders_none_as_NA():
    r = classify("x", "ACGT" * 10)          # too_short -> all None
    fields = format_row(r).split("\t")
    assert len(fields) == 25
    assert fields[3] == "NA"

def test_one_row_per_record_always(tmp_path):
    ltr = _rnd(300, 1)
    recs = [("good", ltr + _rnd(800, 2) + ltr),
            ("random", _rnd(2000, 3)),
            ("short", "ACGT" * 5),
            ("allN", "N" * 500)]
    p = _fasta(tmp_path, recs)
    out = io.StringIO()
    n = run(p, out, threads=1)
    lines = out.getvalue().rstrip("\n").split("\n")
    assert n == 4
    assert len(lines) == 5                   # header + 4 records
    assert lines[0].split("\t") == COLUMNS

def test_output_order_matches_input_order_with_threads(tmp_path):
    recs = []
    for i in range(40):
        ltr = _rnd(200, i)
        recs.append((f"r{i}", ltr + _rnd(400, i + 100) + ltr))
    p = _fasta(tmp_path, recs)
    out = io.StringIO()
    run(p, out, threads=4)
    ids = [l.split("\t")[0] for l in out.getvalue().rstrip("\n").split("\n")[1:]]
    assert ids == [f"r{i}" for i in range(40)]

def test_duplicate_ids_both_appear(tmp_path):
    ltr = _rnd(300, 5)
    p = _fasta(tmp_path, [("dup", ltr + _rnd(600, 6) + ltr),
                          ("dup", ltr + _rnd(600, 7) + ltr)])
    out = io.StringIO()
    run(p, out, threads=1)
    ids = [l.split("\t")[0] for l in out.getvalue().rstrip("\n").split("\n")[1:]]
    assert ids == ["dup", "dup"]

def test_resume_skips_completed_records(tmp_path):
    recs = [(f"r{i}", _rnd(200, i) + _rnd(400, i+50) + _rnd(200, i)) for i in range(10)]
    p = _fasta(tmp_path, recs)
    out = io.StringIO()
    run(p, out, threads=1, resume_skip=6)
    lines = out.getvalue().rstrip("\n").split("\n")
    assert lines[0].split("\t")[0] == "r6"   # no header when resuming
    assert len(lines) == 4

def test_count_data_lines_excludes_header(tmp_path):
    p = tmp_path / "o.tsv"
    p.write_text("\t".join(COLUMNS) + "\nrow1\nrow2\n")
    assert count_data_lines(p) == 2
    assert count_data_lines(tmp_path / "missing.tsv") == 0

def test_gzip_input(tmp_path):
    ltr = _rnd(300, 8)
    p = tmp_path / "in.fa.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(f">g\n{ltr + _rnd(600, 9) + ltr}\n")
    out = io.StringIO()
    assert run(p, out, threads=1) == 1


def test_resume_with_zero_completed_rows_does_not_duplicate_header(tmp_path):
    """A job killed after the header flushed but before the first record has
    resume_skip == 0 while genuinely resuming. Inferring 'fresh run' from the
    count writes a second header into the middle of the data."""
    ltr = _rnd(300, 1)
    p = _fasta(tmp_path, [(f"r{i}", ltr + _rnd(600, i) + ltr) for i in range(5)])
    full = io.StringIO()
    run(p, full, threads=1)
    expected = full.getvalue()

    out = io.StringIO()
    out.write("\t".join(COLUMNS) + "\n")          # header already flushed, 0 data rows
    run(p, out, threads=1, resume_skip=0, resuming=True)
    assert out.getvalue() == expected
    assert sum(1 for l in out.getvalue().splitlines() if l.startswith("seq_id\t")) == 1


def test_parallel_path_streams_and_does_not_materialise_input(tmp_path):
    """ProcessPoolExecutor.map drains its input generator entirely before
    returning, which read the whole FASTA into memory. Bounded submission must
    consume the input incrementally."""
    import ltrk2p.runner as R
    ltr = _rnd(200, 2)
    p = _fasta(tmp_path, [(f"r{i}", ltr + _rnd(300, i) + ltr) for i in range(60)])
    orig = R.read_fasta
    pulled = [0]
    seen_at_first_row = []

    def counting(path):
        for rec in orig(path):
            pulled[0] += 1
            yield rec

    class H(io.StringIO):
        def write(self, s):
            if not s.startswith("seq_id"):
                seen_at_first_row.append(pulled[0])
            return super().write(s)

    R.read_fasta = counting
    try:
        run(p, H(), threads=4)
    finally:
        R.read_fasta = orig
    assert seen_at_first_row[0] < 60, (
        f"{seen_at_first_row[0]}/60 records consumed before the first row was written "
        "-- the input is being materialised, not streamed")


# --------------------------------------------------------------------------- #
# Auxiliary output streams
# --------------------------------------------------------------------------- #

def _elements(n, flank=40):
    """n flanked elements, plus the raw text so slices can be checked."""
    recs = []
    for i in range(n):
        ltr = "TG" + _rnd(298, i) + "CA"
        recs.append((f"chr1:{1000 + i * 9000}-x",
                     _rnd(flank, 500 + i) + ltr + _rnd(800, 600 + i) + ltr
                     + _rnd(flank, 700 + i)))
    return recs


def test_the_tsv_is_unchanged_by_asking_for_auxiliary_output(tmp_path):
    """The extras path reads the raw sequence and runs `_classify` instead of
    `classify`. Neither may perturb a single byte of the table -- otherwise a
    run with `--ltr-cluster` would not be comparable to one without."""
    from ltrk2p.extras import ExtraSpec, ExtraWriter
    inp = _fasta(tmp_path, _elements(6))
    plain = io.StringIO()
    run(inp, plain, threads=1)
    with_extras = io.StringIO()
    spec = ExtraSpec(consensus=True, trimmed=True, perfect=("5p",))
    with ExtraWriter(str(tmp_path / "aux"), spec) as w:
        run(inp, with_extras, threads=1, spec=spec, writer=w)
    assert with_extras.getvalue() == plain.getvalue()


def test_auxiliary_records_stay_in_input_order_across_thread_counts(tmp_path):
    from ltrk2p.extras import ExtraSpec, ExtraWriter
    inp = _fasta(tmp_path, _elements(12))
    seen = []
    for threads in (1, 4):
        out = io.StringIO()
        spec = ExtraSpec(trimmed=True)
        base = str(tmp_path / f"aux{threads}")
        with ExtraWriter(base, spec) as w:
            run(inp, out, threads=threads, spec=spec, writer=w)
        seen.append((out.getvalue(), (tmp_path / f"aux{threads}.trimmed.fa").read_text()))
    assert seen[0] == seen[1]
    ids = [ln[1:] for ln in seen[0][1].split("\n") if ln.startswith(">")]
    assert ids == sorted(ids, key=lambda s: int(s.split(":")[1].split("-")[0]))


def test_auxiliary_streams_only_receive_passing_records(tmp_path):
    from ltrk2p.extras import ExtraSpec, ExtraWriter
    good = _elements(3)
    records = good + [("junk", _rnd(2000, 90)), ("tiny", "ACGT" * 10)]
    inp = _fasta(tmp_path, records)
    out = io.StringIO()
    spec = ExtraSpec(trimmed=True)
    with ExtraWriter(str(tmp_path / "aux"), spec) as w:
        n = run(inp, out, threads=1, spec=spec, writer=w)
    assert n == len(records), "the TSV still reports one row per input record"
    statuses = [r.split("\t")[2] for r in out.getvalue().rstrip("\n").split("\n")[1:]]
    assert w.counts["trimmed"] == statuses.count("pass") == 3


def test_scan_output_counts_only_newline_terminated_lines(tmp_path):
    """A run killed mid-write leaves a partial final line. Counting it as
    complete makes --resume skip a record that was never written, and the
    returned offset is what lets the caller truncate the fragment away."""
    from ltrk2p.runner import scan_output, count_data_lines
    cases = {
        "missing.tsv": None,
        "empty.tsv": "",
        "header_only.tsv": "h\n",
        "two_rows.tsv": "h\na\nb\n",
        "truncated.tsv": "h\na\nbb",
    }
    expect = {
        "missing.tsv": (0, 0, False),
        "empty.tsv": (0, 0, False),
        "header_only.tsv": (0, 2, True),
        "two_rows.tsv": (2, 6, True),
        "truncated.tsv": (1, 4, True),
    }
    for name, text in cases.items():
        p = tmp_path / name
        if text is not None:
            p.write_text(text)
        assert scan_output(p) == expect[name], name
        assert count_data_lines(p) == expect[name][0], name
