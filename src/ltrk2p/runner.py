"""Parallel driver: streams records, preserves input order, one row per record."""
from __future__ import annotations

import sys
from collections import deque
from concurrent.futures import ProcessPoolExecutor
from dataclasses import astuple, fields
from itertools import islice
from pathlib import Path

from .align import Result, classify
from .fasta import read_fasta

COLUMNS = [f.name for f in fields(Result)]


def _fmt(v) -> str:
    if v is None:
        return "NA"
    if isinstance(v, float):
        # `+ 0.0` normalises negative zero. K2P of an identical pair evaluates to
        # -0.5*log(1) - 0.25*log(1) == -0.0, which formats as the string "-0" --
        # numerically equal to zero but a needless surprise in a data column
        # (1,767 of 10,307 arabidopsis rows in the Task 16 run).
        return f"{v + 0.0:.6g}"
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
        resume_skip: int = 0, verbose: bool = False, resuming: bool = False,
        **classify_kw) -> int:
    records = read_fasta(input_path)
    if resume_skip:
        records = islice(records, resume_skip, None)
    # resume_skip > 0 can only happen when resuming, so it implies it. The explicit
    # flag additionally covers resume_skip == 0, where the header was already flushed
    # but no record finished -- the case that previously wrote a duplicate header.
    if not (resuming or resume_skip > 0):
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
            elif not verbose and n % 25000 == 0:
                print(f"  {n} records", file=sys.stderr)
    else:
        max_inflight = max(1, threads) * 4
        with ProcessPoolExecutor(max_workers=threads) as pool:
            it = iter(tasks)
            pending = deque(pool.submit(_work, t) for t in islice(it, max_inflight))
            while pending:
                out_handle.write(format_row(pending.popleft().result()) + "\n")
                n += 1
                if verbose and n % 1000 == 0:
                    print(f"  {n} records", file=sys.stderr)
                elif not verbose and n % 25000 == 0:
                    print(f"  {n} records", file=sys.stderr)
                nxt = next(it, None)
                if nxt is not None:
                    pending.append(pool.submit(_work, nxt))
    return n
