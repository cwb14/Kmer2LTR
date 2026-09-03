"""Parallel driver: streams records, preserves input order, one row per record."""
from __future__ import annotations

import sys
from collections import deque
from concurrent.futures import ProcessPoolExecutor
from dataclasses import astuple, fields
from itertools import islice
from pathlib import Path

from .align import Result, classify, _classify
from .extras import ExtraSpec, ExtraWriter, build
from .fasta import read_fasta, read_fasta_raw, read_headers
from .genome import Options, annotate, credit, harvest, locus, orient

COLUMNS = [f.name for f in fields(Result)]


def _fmt(v) -> str:
    if v is None:
        return "NA"
    if isinstance(v, float):
        # `+ 0.0` normalises negative zero. K2P of an identical pair evaluates to
        # -0.5*log(1) - 0.25*log(1) == -0.0, which formats as the string "-0" --
        # numerically equal to zero but a needless surprise in a data column
        # (1,767 of 10,307 rows on the arabidopsis dataset).
        return f"{v + 0.0:.6g}"
    return str(v)


def format_row(r: Result) -> str:
    return "\t".join(_fmt(v) for v in astuple(r))


def scan_output(path) -> tuple[int, int, bool]:
    """`(complete_data_lines, byte_offset_after_them, has_header)` for a TSV.

    Only NEWLINE-TERMINATED lines are counted. A run killed mid-write leaves a
    partial final line; counting it as complete makes `--resume` skip a record
    that was never written and then append onto the fragment, producing one
    lost record and one row of two records concatenated. The offset lets the
    caller truncate that fragment away before appending.

    `has_header` is false for a missing or empty file, which is what stops a
    `--resume` against a file that does not exist yet from writing a headerless
    TSV -- silently costing the first data row to every `csv.DictReader`
    downstream.
    """
    p = Path(path)
    if not p.exists():
        return 0, 0, False
    lines = offset = 0
    with open(p, "rb") as fh:
        for raw in fh:
            if not raw.endswith(b"\n"):
                break
            lines += 1
            offset += len(raw)
    return max(0, lines - 1), offset, lines > 0


def count_data_lines(path) -> int:
    """Complete data lines already written to an output TSV."""
    return scan_output(path)[0]


def _work(args):
    """One record: a TSV row, plus the auxiliary records if any were asked for.

    The extras are built HERE rather than in the parent so only the records
    actually wanted travel back; sending the aligned pair and the raw sequence
    home to build them there would be the larger payload of the two.

    Order matters at the end: `annotate` settles the record's orientation, and
    `build` needs it to shift a `--trim-flanks` header the right way round.
    """
    seq_id, seq, raw, spec, window, gopt, kw = args
    ctx = orient(seq, window) if window is not None else None
    if ctx is not None:
        kw = {**kw, "tsd_credit": credit(seq, ctx, gopt)}
    if spec is None:
        result, aln = classify(seq_id, seq, **kw), None
    else:
        result, aln = _classify(seq_id, seq, **kw)
    if ctx is not None:
        result = annotate(result, seq, ctx, gopt)
    return result, (build(result, aln, raw, spec) if spec is not None else None)


def genome_windows(input_path, genome) -> dict:
    """One streaming pass over the reference, for every locus the input names.

    A header-only pre-pass first, because the reference has to be read in one
    sweep -- a gzipped FASTA cannot be seeked -- and so every locus must be known
    before it starts. The pre-pass builds no sequences, and the reference costs
    one sequential read at a bounded memory whatever its size.

    What is held is the harvest, and it scales with the number of loci rather
    than with the reference: about 700 bytes per record, so a 10,000-record set
    costs single-digit megabytes and a million-record library about 700 MB.
    """
    loci = {locus(sid) for sid in read_headers(input_path)}
    loci.discard(None)
    if not loci:
        print("Kmer2LTR: warning: --genome was given but no header carries a "
              "chrom:start-end locus; the genome columns will all be NA",
              file=sys.stderr)
        return {}
    windows = harvest(genome, loci)
    # Unconditional: reading a 3 Gbp reference is a milestone, and the ratio is
    # how a user finds out that their headers and their reference disagree about
    # sequence names -- which otherwise shows up only as a column of NA.
    print(f"Kmer2LTR: located {len(windows)}/{len(loci)} loci in the reference",
          file=sys.stderr)
    if not windows:
        print("Kmer2LTR: warning: no locus matched a contig in the reference; "
              "check that the headers and the reference share sequence names",
              file=sys.stderr)
    return windows


def run(input_path, out_handle, threads: int = 1, cs: bool = False,
        resume_skip: int = 0, verbose: bool = False, resuming: bool = False,
        spec: ExtraSpec | None = None, writer: ExtraWriter | None = None,
        genome=None, gopt: Options | None = None, **classify_kw) -> int:
    # The raw sequence is read, and sent to the workers, only when something
    # needs it -- a run with no auxiliary output moves exactly the bytes it
    # always did.
    spec = spec or None
    windows = genome_windows(input_path, genome) if genome else {}
    gopt = gopt or Options()
    records = (read_fasta_raw(input_path) if spec
               else ((sid, seq, "") for sid, seq in read_fasta(input_path)))
    if resume_skip:
        records = islice(records, resume_skip, None)
    # resume_skip > 0 can only happen when resuming, so it implies it. The explicit
    # flag additionally covers resume_skip == 0, where the header was already flushed
    # but no record finished -- the case that previously wrote a duplicate header.
    if not (resuming or resume_skip > 0):
        out_handle.write("\t".join(COLUMNS) + "\n")

    kw = {"cs": cs, **classify_kw}
    tasks = ((sid, seq, raw, spec, windows.get(locus(sid)) if windows else None,
              gopt, kw) for sid, seq, raw in records)
    n = 0

    def emit(payload) -> None:
        nonlocal n
        result, extras = payload
        out_handle.write(format_row(result) + "\n")
        if writer is not None:
            writer.write(extras)
        n += 1
        every = 1000 if verbose else 25000
        if n % every == 0:
            print(f"  {n} records", file=sys.stderr)

    if threads <= 1:
        for t in tasks:
            emit(_work(t))
    else:
        max_inflight = max(1, threads) * 4
        with ProcessPoolExecutor(max_workers=threads) as pool:
            it = iter(tasks)
            pending = deque(pool.submit(_work, t) for t in islice(it, max_inflight))
            while pending:
                emit(pending.popleft().result())
                nxt = next(it, None)
                if nxt is not None:
                    pending.append(pool.submit(_work, nxt))
    return n
