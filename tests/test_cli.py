import random, subprocess, sys
import pytest
from ltrk2p.cli import main

PY = "/anvil/projects/x-bio250178/conda/envs/ltrrt4/bin/python"

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

def test_cs_flag_changes_last_column(tmp_path):
    inp = _fa(tmp_path, n=2)
    a, b = tmp_path / "a.tsv", tmp_path / "b.tsv"
    main([str(inp), "-o", str(a)])
    main([str(inp), "--cs", "-o", str(b)])
    last_a = a.read_text().rstrip("\n").split("\n")[1].split("\t")[-1]
    last_b = b.read_text().rstrip("\n").split("\n")[1].split("\t")[-1]
    assert last_a != last_b

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
