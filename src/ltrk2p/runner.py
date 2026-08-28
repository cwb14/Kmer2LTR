"""Parallel driver: streams records, preserves input order, one row per record."""
from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor
from dataclasses import astuple, fields
from itertools import islice
from pathlib import Path

from .align import Result, classify
from .fasta import read_fasta

COLUMNS = [f.name for f in fields(Result)]
CHUNK = 64


def _fmt(v) -> str:
    if v is None:
        return "NA"
    if isinstance(v, float):
        return f"{v:.6g}"
    return str(v)


def format_row(r: Result) -> str:
    return "\t".join(_fmt(v) for v in astuple(r))


def count_data_lines(path) -> int:
    p = Path(path)
    if not p.exists():
        return 0
    with open(p) as fh:
        return max(0, sum(1 for _ in fh) - 1)


def _work(args):
    seq_id, seq, kw = args
    return classify(seq_id, seq, **kw)


def run(input_path, out_handle, threads: int = 1, cs: bool = False,
        resume_skip: int = 0, verbose: bool = False, **classify_kw) -> int:
    records = read_fasta(input_path)
    if resume_skip:
        records = islice(records, resume_skip, None)
    else:
        out_handle.write("\t".join(COLUMNS) + "\n")

    kw = {"cs": cs, **classify_kw}
    tasks = ((sid, seq, kw) for sid, seq in records)
    n = 0
    if threads <= 1:
        for t in tasks:
            out_handle.write(format_row(_work(t)) + "\n")
            n += 1
            if verbose and n % 1000 == 0:
                print(f"  {n} records", file=sys.stderr)
    else:
        with ProcessPoolExecutor(max_workers=threads) as pool:
            # imap-style ordered results; chunking amortises IPC per record
            for res in pool.map(_work, tasks, chunksize=CHUNK):
                out_handle.write(format_row(res) + "\n")
                n += 1
                if verbose and n % 1000 == 0:
                    print(f"  {n} records", file=sys.stderr)
    return n
