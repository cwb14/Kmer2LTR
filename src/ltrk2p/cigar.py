"""Alignment string emission: extended CIGAR and minimap2 cs.

Orientation convention, fixed here and relied on everywhere downstream:
query = 5' LTR, ref = 3' LTR. I = base in the 5' LTR absent from the 3' LTR.
"""
from __future__ import annotations

# pywfa op codes. Note pywfa emits 0 ("M") for TRUE matches and 8 ("X") for
# mismatches -- a mixed convention. Because X is emitted separately, 0
# unambiguously means "equal" here, contrary to the SAM spec's meaning of M.
# NOTE: When calling pywfa as WavefrontAligner(query)(ref), binding pattern=query
# and text=ref, pywfa's op codes are: 1=deletion (in query), 2=insertion (in query).
# This is backwards from the SAM spec (where 1=I, 2=D). The mapping holds for this
# calling convention; a caller passing sequences in the opposite order would flip it.
_OP_EQ, _OP_INS, _OP_DEL, _OP_X = 0, 2, 1, 8


def aligned_pair_from_cigartuples(tuples, query: str, ref: str) -> tuple[str, str]:
    qi = ri = 0
    qa: list[str] = []
    ra: list[str] = []
    for op, n in tuples:
        if op in (_OP_EQ, _OP_X):
            qa.append(query[qi:qi + n]); ra.append(ref[ri:ri + n]); qi += n; ri += n
        elif op == _OP_INS:
            qa.append(query[qi:qi + n]); ra.append("-" * n); qi += n
        elif op == _OP_DEL:
            qa.append("-" * n); ra.append(ref[ri:ri + n]); ri += n
        else:
            raise ValueError(f"unexpected cigar op code {op}")
    return "".join(qa), "".join(ra)


def _ops(a: str, b: str):
    if len(a) != len(b):
        raise ValueError(f"aligned strings differ in length: {len(a)} vs {len(b)}")
    for x, y in zip(a, b):
        if x == "-":
            yield "D", x, y
        elif y == "-":
            yield "I", x, y
        elif x == y:
            yield "=", x, y
        else:
            yield "X", x, y


def extended_cigar(a: str, b: str) -> str:
    out: list[str] = []
    run_op: str | None = None
    run_n = 0
    for op, _, _ in _ops(a, b):
        if op == run_op:
            run_n += 1
        else:
            if run_op is not None:
                out.append(f"{run_n}{run_op}")
            run_op, run_n = op, 1
    if run_op is not None:
        out.append(f"{run_n}{run_op}")
    return "".join(out)


def cs_string(a: str, b: str) -> str:
    """minimap2 short cs. Substitutions are *<ref><query>, both lowercase."""
    out: list[str] = []
    match_run = 0
    pending_op: str | None = None
    pending: list[str] = []

    def flush_pending():
        nonlocal pending_op, pending
        if pending_op is not None:
            out.append(pending_op + "".join(pending).lower())
            pending_op, pending = None, []

    for op, x, y in _ops(a, b):
        if op == "=":
            flush_pending()
            match_run += 1
            continue
        if match_run:
            out.append(f":{match_run}")
            match_run = 0
        if op == "X":
            flush_pending()
            out.append(f"*{y.lower()}{x.lower()}")
        elif op == "I":
            if pending_op != "+":
                flush_pending(); pending_op = "+"
            pending.append(x)
        else:  # D
            if pending_op != "-":
                flush_pending(); pending_op = "-"
            pending.append(y)
    flush_pending()
    if match_run:
        out.append(f":{match_run}")
    return "".join(out)
