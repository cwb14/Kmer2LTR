"""Per-element derived sequences: consensus LTR, internal region, perfect
elements, boundary-corrected element.

Everything here is built from measurements `align` has already made -- the
settled boundaries and the final aligned LTR pair -- so nothing is re-aligned
and no derived record can disagree with the TSV row beside it.

**Slices of the input reproduce the input verbatim**, case and IUPAC codes
intact, because they are cut from `fasta.read_fasta_raw`'s original string
rather than from the sanitised one the alignment ran on. A `--trim-flanks`
output that silently uppercased a soft-masked genome would be a worse file than
the one it replaced. The IUPAC consensus is the single synthesised sequence, and
it is uppercase.
"""
from __future__ import annotations

import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

WRAP = 60

# Two-base IUPAC ambiguity codes. There is no three- or four-base case: a
# pairwise alignment column holds at most two observations.
_IUPAC_PAIR = {
    frozenset("AG"): "R",   # purine
    frozenset("CT"): "Y",   # pyrimidine
    frozenset("GC"): "S",   # strong
    frozenset("AT"): "W",   # weak
    frozenset("GT"): "K",   # keto
    frozenset("AC"): "M",   # amino
}
_BASES = frozenset("ACGT")

PERFECT_MODES = ("5p", "3p", "consensus")


def iupac_consensus(a: str, b: str) -> str:
    """IUPAC consensus of two equal-length gapped aligned strings.

    Per column: equal bases keep the base; two different bases become the
    two-way ambiguity code; a gap on one side yields the other side's base;
    anything else yields `N`. A column gapped on BOTH sides cannot occur in a
    pairwise alignment and so is not a case here.

    The result is as long as the alignment, which is at least as long as either
    LTR -- an indel in one copy contributes the other copy's base rather than
    being dropped, so no observed base is discarded.
    """
    if len(a) != len(b):
        raise ValueError(f"aligned strings differ in length: {len(a)} vs {len(b)}")
    out: list[str] = []
    for x, y in zip(a, b):
        if x == y:
            out.append(x if x in _BASES else "N")
        elif x == "-":
            out.append(y if y in _BASES else "N")
        elif y == "-":
            out.append(x if x in _BASES else "N")
        elif x in _BASES and y in _BASES:
            out.append(_IUPAC_PAIR[frozenset((x, y))])
        else:
            out.append("N")
    return "".join(out)


# `chrom:start-end`, also accepting LTR_retriever's `chrom:start..end`. `(.+)` is
# greedy so it binds the LAST colon, which is what keeps a chromosome name that
# itself contains a colon intact. `\Z` rather than `$`, which would also match
# before a trailing newline and silently eat it.
_LOCUS = re.compile(r"(.+):(\d+)(-|\.\.)(\d+)\Z")


def _locate_locus(seq_id: str):
    """`(match, start, end)` for the locus field inside a header, or None.

    Two placements are recognised, which between them cover what the common
    producers emit:

        chr1:1000-2000#LTR/Gypsy          locus first, class tag after
        TE_1#LTR/Copia::chr1:1000-2000    `bedtools getfasta -name`, locus last

    The second matters because the `#` sits in the MIDDLE there: reading only up
    to the first `#` finds no locus, and a header left unchanged over a sequence
    that WAS trimmed is a header that lies about its own span.
    """
    head_end = seq_id.find("#")
    if head_end == -1:
        head_end = len(seq_id)
    m = _LOCUS.match(seq_id, 0, head_end)
    if m:
        return m, 0, head_end
    tail = seq_id.rfind("::")
    if tail != -1:
        m = _LOCUS.match(seq_id, tail + 2)
        if m:
            return m, tail + 2, len(seq_id)
    return None


def shift_locus(seq_id: str, trim5: int, trim3: int) -> str:
    """Move a `chrom:start-end` header inwards by the trimmed flank lengths.

    Everything around the locus -- a RepeatMasker-style class tag, a
    `bedtools`-style name prefix -- is carried through untouched. Headers with
    no parseable locus are returned unchanged: this is the one place in Kmer2LTR
    that reads a header as anything but an opaque id, and it declines rather
    than guesses.
    """
    found = _locate_locus(seq_id)
    if found is None:
        return seq_id
    m, lo, hi = found
    chrom, start, dash, end = m.group(1), int(m.group(2)), m.group(3), int(m.group(4))
    start += trim5
    end -= trim3
    if start > end:
        return seq_id
    return f"{seq_id[:lo]}{chrom}:{start}{dash}{end}{seq_id[hi:]}"


def fasta_record(header: str, seq: str, wrap: int = WRAP) -> str:
    """One FASTA record, sequence wrapped at `wrap` columns."""
    lines = [seq[i:i + wrap] for i in range(0, len(seq), wrap)] or [""]
    return ">" + header + "\n" + "\n".join(lines) + "\n"


@dataclass(frozen=True)
class ExtraSpec:
    """Which auxiliary records the run wants."""
    consensus: bool = False
    internal: bool = False
    trimmed: bool = False
    perfect: tuple[str, ...] = ()

    def __bool__(self) -> bool:
        return bool(self.consensus or self.internal or self.trimmed or self.perfect)

    @property
    def wants_consensus(self) -> bool:
        """The consensus sequence is needed, whether or not it is written out."""
        return self.consensus or "consensus" in self.perfect


@dataclass(frozen=True)
class Extras:
    """The auxiliary FASTA records for one element, ready to write."""
    seq_id: str
    consensus: str | None = None
    internal: str | None = None
    trimmed: str | None = None
    perfect: dict[str, str] = field(default_factory=dict)


def build(result, aln, raw: str, spec: ExtraSpec) -> Extras | None:
    """Auxiliary records for one classified element, or None if it has none.

    Only `status == "pass"` contributes. A weak or saturated pair still carries
    real coordinates, but the claim these outputs make -- this is an element,
    here is its family, here is when it inserted -- is exactly the claim the
    significance gate declined to support.

    `result` coordinates are 1-based inclusive into the sequence; `raw` indexes
    identically to the sanitised sequence they were measured on.
    """
    if result.status != "pass":
        return None
    l5s, l5e = result.ltr5_start, result.ltr5_end
    l3s, l3e = result.ltr3_start, result.ltr3_end

    # Both are sliced/derived only when something asks for them: on a long
    # element the internal region is most of the record.
    internal = raw[l5e:l3s - 1] if (spec.internal or spec.perfect) else ""
    cons = iupac_consensus(*aln) if (aln and spec.wants_consensus) else None

    perfect: dict[str, str] = {}
    if spec.perfect and internal:
        flanks = {"5p": raw[l5s - 1:l5e], "3p": raw[l3s - 1:l3e], "consensus": cons}
        for mode in spec.perfect:
            flank = flanks.get(mode)
            if flank:
                # `~LTRlen:<n>` is Kmer2LTR's convention, kept so its
                # `lib_mutator.py` still reads these files.
                perfect[mode] = fasta_record(
                    f"{result.seq_id}~LTRlen:{len(flank)}", flank + internal + flank)

    return Extras(
        seq_id=result.seq_id,
        consensus=(fasta_record(result.seq_id, cons)
                   if (spec.consensus and cons) else None),
        internal=(fasta_record(result.seq_id, internal)
                  if (spec.internal and internal) else None),
        trimmed=(fasta_record(
            shift_locus(result.seq_id, result.flank5_len, result.flank3_len),
            raw[l5s - 1:l3e]) if spec.trimmed else None),
        perfect=perfect,
    )


def stream_paths(base: str, spec: ExtraSpec) -> dict[str, Path]:
    """Every file `ExtraWriter(base, spec)` would open, without opening any.

    Separate from the writer so the CLI can check the set for collisions with
    the input before a single handle is created -- opening truncates, and
    truncating the file you are about to read destroys it.
    """
    paths: dict[str, Path] = {}
    if spec.consensus:
        paths["consensus"] = Path(f"{base}.consensus.fa")
    if spec.internal:
        paths["internal"] = Path(f"{base}.internal.fa")
    if spec.trimmed:
        paths["trimmed"] = Path(f"{base}.trimmed.fa")
    for mode in spec.perfect:
        paths[f"perfect_{mode}"] = Path(f"{base}.perfect_{mode}.fa")
    return paths


class ExtraWriter:
    """Open handles for the auxiliary FASTA streams, written in input order.

    A stream this writer CREATED and then never wrote to is deleted on close,
    so an input with no passing element leaves no misleading zero-byte FASTA.
    A file that already existed is only truncated, never removed: deleting a
    file the run did not create is a worse outcome than leaving an empty one,
    and the empty file at least describes the run that just happened.
    """

    def __init__(self, base: str, spec: ExtraSpec):
        self.paths = stream_paths(base, spec)
        self.counts = dict.fromkeys(self.paths, 0)
        self._ours = {k: not p.exists() for k, p in self.paths.items()}
        self._fh: dict[str, object] = {}
        try:
            for key, path in self.paths.items():
                self._fh[key] = open(path, "w")
        except OSError:
            # Half-open is not a state anything downstream can use, and the
            # caller never gets the object back to close it.
            for fh in self._fh.values():
                fh.close()
            self._fh.clear()
            raise
        self._seen: set[str] = set()
        self._warned = False

    def _put(self, key: str, record: str | None) -> None:
        if record and key in self._fh:
            self._fh[key].write(record)
            self.counts[key] += 1

    def write(self, extras: Extras | None) -> None:
        if extras is None:
            return
        # The TSV is positional, so duplicate ids are harmless there. These
        # files are keyed by id -- an mmseqs cluster naming a duplicated id
        # cannot be joined back to one element -- so say so once.
        if not self._warned:
            if extras.seq_id in self._seen:
                print(f"Kmer2LTR: warning: duplicate record id {extras.seq_id!r} in the "
                      f"auxiliary FASTA output; ids there are not unique",
                      file=sys.stderr)
                self._warned = True
            else:
                self._seen.add(extras.seq_id)
        self._put("consensus", extras.consensus)
        self._put("internal", extras.internal)
        self._put("trimmed", extras.trimmed)
        for mode, record in extras.perfect.items():
            self._put(f"perfect_{mode}", record)

    def close(self) -> None:
        for fh in self._fh.values():
            fh.close()
        self._fh.clear()
        self._seen.clear()
        for key, path in list(self.paths.items()):
            if not self.counts[key] and self._ours[key]:
                path.unlink(missing_ok=True)
                del self.paths[key]

    def __enter__(self) -> "ExtraWriter":
        return self

    def __exit__(self, *exc) -> None:
        self.close()
