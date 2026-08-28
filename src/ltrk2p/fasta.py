"""Streaming FASTA input and sequence sanitisation."""
from __future__ import annotations

import gzip
from pathlib import Path
from typing import Iterator

class _Table(dict):
    """Translation table covering every codepoint.

    ACGT (either case) -> uppercase; any whitespace -> dropped; everything
    else (IUPAC codes, *, -, digits, non-ASCII) -> N. Using __missing__
    rather than a prebuilt 256-entry table means codepoints >= 256 and the
    whitespace-classified control codes are covered too.
    """

    def __missing__(self, key):
        return None if chr(key).isspace() else "N"


_TABLE = _Table({ord(c): c.upper() for c in "ACGTacgt"})


def sanitize(seq: str) -> str:
    """Uppercase, drop whitespace, map every non-ACGT character to N."""
    return seq.translate(_TABLE)


def _open(path) -> Iterator[str]:
    path = Path(path)
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as fh:
            yield from fh
    else:
        with open(path, "rt") as fh:
            yield from fh


def read_fasta(path) -> Iterator[tuple[str, str]]:
    """Yield (seq_id, sanitized_seq). seq_id is the header up to first whitespace.

    Streaming: only one record is held in memory at a time. Duplicate IDs are
    yielded unchanged -- the one-row-per-record output contract makes them safe.
    """
    seq_id: str | None = None
    chunks: list[str] = []
    for line in _open(path):
        if line.startswith(">"):
            if seq_id is not None:
                yield seq_id, sanitize("".join(chunks))
            seq_id = line[1:].strip().split(None, 1)[0] if line[1:].strip() else ""
            chunks = []
        elif seq_id is not None:
            chunks.append(line)
    if seq_id is not None:
        yield seq_id, sanitize("".join(chunks))
