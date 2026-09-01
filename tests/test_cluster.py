import subprocess
from pathlib import Path

import pytest

from ltrk2p import cluster

_FASTA = ">a\nACGTACGTACGTACGTACGTACGTACGTACGT\n>b\nACGTACGTACGTACGTACGTACGTACGTACGA\n"


def _fa(tmp_path, name="out.consensus.fa", text=_FASTA):
    p = tmp_path / name
    p.write_text(text)
    return p


def _stub(calls, fail_at=(), write_tsv=True):
    """Stand in for mmseqs: record the command and lay down its output files."""
    def run(cmd, **kw):
        calls.append(cmd)
        prefix = Path(cmd[3])
        tag = cmd[cmd.index("--min-seq-id") + 1]
        if tag in fail_at:
            raise subprocess.CalledProcessError(1, cmd)
        for aux in ("_all_seqs.fasta", "_rep_seq.fasta"):
            Path(f"{prefix}{aux}").write_text(">a\nACGT\n")
        if write_tsv:
            Path(f"{prefix}_cluster.tsv").write_text("a\ta\na\tb\n")
        return subprocess.CompletedProcess(cmd, 0)
    return run


def test_sweep_runs_every_identity_and_names_each_table(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(cluster.subprocess, "run", _stub(calls))
    tsvs = cluster.cluster(_fa(tmp_path), threads=4)
    assert len(tsvs) == len(cluster.SWEEP)
    assert [t.name for t in tsvs] == [
        f"out.consensus_id{v:.2f}_cluster.tsv" for v in cluster.SWEEP]


def test_a_pinned_identity_runs_once(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(cluster.subprocess, "run", _stub(calls))
    tsvs = cluster.cluster(_fa(tmp_path), threads=1, min_seq_id=0.85)
    assert len(calls) == 1 and len(tsvs) == 1
    assert tsvs[0].name == "out.consensus_id0.85_cluster.tsv"


def test_the_grid_searched_parameters_are_passed_verbatim(tmp_path, monkeypatch):
    """These were chosen by a grid search on Arabidopsis annotations against
    singleton count and family mixing. Drifting any of them silently changes
    what a "family" means in the output."""
    calls = []
    monkeypatch.setattr(cluster.subprocess, "run", _stub(calls))
    cluster.cluster(_fa(tmp_path), threads=7, min_seq_id=0.9)
    cmd = calls[0]
    assert cmd[:2] == ["mmseqs", "easy-cluster"]
    for flag, value in (("--min-seq-id", "0.90"), ("-c", "0.5"), ("--cov-mode", "0"),
                        ("--cluster-mode", "1"), ("--mask", "0"), ("-s", "7.5"),
                        ("--threads", "7")):
        assert cmd[cmd.index(flag) + 1] == value, f"{flag} changed"


def test_mmseqs_scratch_never_survives_a_run(tmp_path, monkeypatch):
    """The tmp directory and the two auxiliary FASTAs are mmseqs' working
    files, not results. Leaving them turns one run into nine files per
    identity."""
    monkeypatch.setattr(cluster.subprocess, "run", _stub([]))
    cluster.cluster(_fa(tmp_path), min_seq_id=0.8)
    left = sorted(p.name for p in tmp_path.iterdir())
    assert left == ["out.consensus.fa", "out.consensus_id0.80_cluster.tsv"], left


def test_scratch_is_removed_even_when_mmseqs_fails(tmp_path, monkeypatch):
    monkeypatch.setattr(cluster.subprocess, "run", _stub([], fail_at=("0.80",)))
    assert cluster.cluster(_fa(tmp_path), min_seq_id=0.8) == []
    assert [p.name for p in tmp_path.iterdir()] == ["out.consensus.fa"]


def test_one_failing_identity_does_not_abort_the_sweep(tmp_path, monkeypatch):
    """The identities are independent, and the run has already paid for the
    alignment -- losing six good tables to one bad one would be the expensive
    mistake."""
    calls = []
    monkeypatch.setattr(cluster.subprocess, "run", _stub(calls, fail_at=("0.80", "0.95")))
    tsvs = cluster.cluster(_fa(tmp_path))
    assert len(calls) == len(cluster.SWEEP)
    assert len(tsvs) == len(cluster.SWEEP) - 2


def test_a_missing_or_empty_fasta_is_a_warning_not_a_crash(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(cluster.subprocess, "run", _stub([]))
    assert cluster.cluster(tmp_path / "absent.fa") == []
    assert cluster.cluster(_fa(tmp_path, "empty.fa", "")) == []
    assert capsys.readouterr().err.count("nothing to cluster") == 2


def test_a_run_that_writes_no_table_is_reported(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(cluster.subprocess, "run", _stub([], write_tsv=False))
    assert cluster.cluster(_fa(tmp_path), min_seq_id=0.7) == []
    assert "wrote no" in capsys.readouterr().err


@pytest.mark.skipif(not cluster.available(), reason="mmseqs not on PATH")
def test_real_mmseqs_clusters_identical_sequences_together(tmp_path):
    """The stubbed tests pin the command; this one pins that the command
    actually does what the rest of the tool assumes it does."""
    ltr = "ACGTTGCATTACGGATCCATTGACCAGTTACGATCGGATTACAGTCCATGGATCAGTTAC" * 4
    fa = _fa(tmp_path, text="".join(
        f">e{i}\n{ltr}\n" for i in range(3)) + f">far\n{'AT' * 120}\n")
    tsvs = cluster.cluster(fa, threads=1, min_seq_id=0.9)
    assert len(tsvs) == 1
    members = {}
    for line in tsvs[0].read_text().splitlines():
        rep, member = line.split("\t")
        members.setdefault(rep, set()).add(member)
    triple = [v for v in members.values() if {"e0", "e1", "e2"} <= v]
    assert triple, f"identical sequences did not co-cluster: {members}"
    assert [p.name for p in tmp_path.iterdir() if p.suffix == ".fasta"] == []


def test_a_failed_identity_does_not_leave_the_previous_run_s_table(tmp_path, monkeypatch):
    """A stale table beside six fresh ones makes the directory look like a
    complete sweep of the current input while one identity silently describes
    a different one."""
    fa = _fa(tmp_path)
    monkeypatch.setattr(cluster.subprocess, "run", _stub([]))
    first = cluster.cluster(fa)
    assert len(first) == len(cluster.SWEEP)
    for t in first:
        t.write_text("STALE\tSTALE\n")

    monkeypatch.setattr(cluster.subprocess, "run", _stub([], fail_at=("0.80", "0.90")))
    second = cluster.cluster(fa)
    assert len(second) == len(cluster.SWEEP) - 2
    stale = [p.name for p in tmp_path.glob("*_cluster.tsv")
             if "STALE" in p.read_text()]
    assert stale == [], f"tables left over from the previous run: {stale}"


def test_identities_reports_what_the_sweep_will_attempt():
    """The CLI compares this against what came back to decide whether the sweep
    finished -- and, from that, whether the internal FASTA is safe to drop."""
    assert cluster.identities(None) == cluster.SWEEP
    assert cluster.identities(0.85) == (0.85,)
