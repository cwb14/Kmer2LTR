"""Streaming FASTA input and sequence sanitisation."""
from __future__ import annotations

import gzip
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


def is_gzip(path) -> bool:
    """True when `path` holds gzip (or bgzip) data, whatever it is named.

    Sniffed rather than read off the suffix because callers rename: LTRquest
    stages every genome as `<prefix>.input_genome.fa` whether the file behind it
    is a plain FASTA or the caller's original `.fa.gz`. Handing a compressed
    stream to a text reader does not fail where the mistake was made -- it
    sanitises binary to N, or raises UnicodeDecodeError several steps later on
    whichever DEFLATE byte happened to look like a header.

    A file that cannot be opened is not gzip; the caller's own open raises the
    error, which says more than one from a two-byte probe would.
    """
    try:
        with open(path, "rb") as fh:
            return fh.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def _open(path) -> Iterator[str]:
    op = gzip.open if is_gzip(path) else open
    with op(path, "rt") as fh:
        yield from fh


def _despace(seq: str) -> str:
    """Drop every whitespace character, keeping all others untouched.

    Uses the same notion of whitespace as `sanitize` (`str.isspace`, which is
    what `str.split()` splits on), so `sanitize(_despace(s))` is the same
    length as `_despace(s)` and the two index identically.
    """
    return "".join(seq.split())


def read_fasta_raw(path) -> Iterator[tuple[str, str, str]]:
    """Yield (seq_id, sanitized_seq, original_seq).

    `original_seq` has had whitespace removed and nothing else: case, IUPAC
    codes and every other character survive. It indexes identically to the
    sanitized sequence, so a coordinate measured on one slices the other.

    That is what the sequence-slicing outputs are written from, so
    `--trim-flanks` and `--perfect-ltr-rt` never silently uppercase a
    soft-masked genome or collapse an ambiguity code to N.

    Streaming: at most one record is held in memory, in both forms.
    """
    seq_id: str | None = None
    chunks: list[str] = []
    for line in _open(path):
        if line.startswith(">"):
            if seq_id is not None:
                raw = _despace("".join(chunks))
                yield seq_id, sanitize(raw), raw
            seq_id = line[1:].strip().split(None, 1)[0] if line[1:].strip() else ""
            chunks = []
        elif seq_id is not None:
            chunks.append(line)
    if seq_id is not None:
        raw = _despace("".join(chunks))
        yield seq_id, sanitize(raw), raw


def read_fasta(path) -> Iterator[tuple[str, str]]:
    """Yield (seq_id, sanitized_seq). seq_id is the header up to first whitespace.

    Streaming: only one record is held in memory at a time. Duplicate IDs are
    yielded unchanged -- the one-row-per-record output contract makes them safe.
    """
    for seq_id, seq, _raw in read_fasta_raw(path):
        yield seq_id, seq


def read_headers(path) -> Iterator[str]:
    """Yield just the record ids, assembling no sequence at all.

    `--genome` has to know every locus before the first record is classified,
    and the sequences are the expensive part of that pass: on a 68 Mbp element
    set they would be built, joined and immediately discarded. Ids are split
    exactly as `read_fasta_raw` splits them, so the two passes agree record for
    record -- a disagreement would harvest one element's locus for another's
    sequence.
    """
    for line in _open(path):
        if line.startswith(">"):
            yield line[1:].strip().split(None, 1)[0] if line[1:].strip() else ""
