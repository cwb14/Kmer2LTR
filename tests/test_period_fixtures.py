"""bench/period_fixtures.py: the driver that scores --period-rule on real data.

Exercised on a synthetic fixture directory built to the same layout, so the
driver is tested without shipping real genomic sequence into the repo.
"""
import random

import pytest

from bench.period_fixtures import load, score

HEADER = ("set\tid\tlen\tltrlen\texpect_ltr5\texpect_ltr3\texpect_internal\t"
          "kmer2ltr_current_call\ttrf_period\ttrf_copies\tsrc_library_entry\tclass\n")


def _rnd(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


def _hard(seed):
    """An element whose longest self-alignment overlaps, not contains, its LTRs."""
    C, M, E = _rnd(100, seed), _rnd(500, seed + 1), _rnd(200, seed + 2)
    D = C + M + C
    return M[400:500] + D + E + D + M[0:100]          # 1800 bp, LTRs 1-300/1501-1800


def _easy(seed):
    ltr = _rnd(300, seed)
    return ltr + _rnd(1200, seed + 1) + ltr           # 1800 bp, LTRs 1-300/1501-1800


def _row(sset, sid):
    return (f"{sset}\t{sid}\t1800\t300\t1-300\t1501-1800\t1200\t"
            f"(not read)\t0\t0.00\tsynthetic\tLTR/unknown/unknown\n")


@pytest.fixture
def fixture_dir(tmp_path):
    (tmp_path / "rebound_failures.fa").write_text(
        "".join(f">hard{i}\n{_hard(300 + 10 * i)}\n" for i in range(2)))
    (tmp_path / "controls_must_not_regress.fa").write_text(
        "".join(f">easy{i}\n{_easy(400 + 10 * i)}\n" for i in range(2)))
    (tmp_path / "expected.tsv").write_text(
        HEADER
        + "".join(_row("rebound_failures", f"hard{i}") for i in range(2))
        + "".join(_row("controls_must_not_regress", f"easy{i}") for i in range(2)))
    return tmp_path


def test_load_reads_both_fastas_and_the_table(fixture_dir):
    seqs, rows = load(fixture_dir)
    assert set(seqs) == {"hard0", "hard1", "easy0", "easy1"}
    assert len(rows) == 4


def test_load_refuses_a_table_naming_an_absent_record(fixture_dir):
    p = fixture_dir / "expected.tsv"
    p.write_text(p.read_text() + _row("rebound_failures", "ghost"))
    with pytest.raises(KeyError, match="ghost"):
        load(fixture_dir)


def test_load_refuses_a_missing_fasta(tmp_path):
    with pytest.raises(FileNotFoundError, match="rebound_failures.fa"):
        load(tmp_path)


def test_outermost_fixes_the_failures_and_leaves_the_controls_alone(fixture_dir):
    s = score(fixture_dir)["summary"]
    assert s["rebound_failures"] == {"n": 2, "best-score": 0, "outermost": 2}
    assert s["controls_must_not_regress"] == {"n": 2, "best-score": 2, "outermost": 2}


def test_the_recorded_baseline_column_is_not_trusted(fixture_dir):
    """`kmer2ltr_current_call` is recomputed, never read: a hand-recorded
    baseline drifts from the code and then becomes the thing being measured."""
    p = fixture_dir / "expected.tsv"
    p.write_text(p.read_text().replace("(not read)", "1-300/1501-1800"))
    s = score(fixture_dir)["summary"]
    assert s["rebound_failures"]["best-score"] == 0
