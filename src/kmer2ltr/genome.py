"""Genomic context for records whose headers carry a locus.

Everything here is optional and off unless `--genome` is given. Without it no
header is parsed for coordinates, no reference is opened, and the four columns
this module fills stay `NA`.

Three things the element FASTA cannot say on its own, and one reference can:

* **Which strand the record is stored on.** Annotation pipelines routinely
  reverse-complement an element while leaving forward genomic coordinates in its
  header, and nothing in the sequence reveals that. It matters because
  `--trim-flanks` shifts those coordinates inwards, and on a reversed record the
  5' trim belongs to the header's `end`.
* **Target-site duplications.** A TSD lies *outside* the element, so without the
  reference it can only be looked for when a flank was called -- exactly the
  records where the boundary is least certain.
* **Whether the supplied record was already perfectly bounded.**

The reference is read as a plain FASTA: one pass per file, in bounded blocks
against a rolling offset, cutting only the few dozen bases wanted around each
locus. A contig is never held whole even when the reference is written unwrapped,
so a 3.2 Gbp reference costs what a 100 Mbp one does -- 27 s at a 38 MB peak.
Whitespace inside a sequence line is dropped rather than kept, because it shifts
the offset every cut is measured against.

The header parser is imported from `extras` rather than rewritten. Those are the
only two places in Kmer2LTR that read a header as anything but an opaque id, and
two regexes for one format is a bug waiting for one of them to be improved.
"""
from __future__ import annotations

import gzip
import sys
from dataclasses import dataclass, replace

from .extras import _locate_locus
from .fasta import is_gzip

PAD = 32          # reference bases kept outside each end of the record
PROBE = 40        # bases compared against the record to settle its orientation
MIN_PROBE = 16    # below this there is not enough sequence to settle anything
MIN_IDENTITY = 0.9


def locus(seq_id: str) -> tuple[str, int, int] | None:
    """`(chrom, start, end)` from a header, 1-based inclusive, or None.

    Recognises exactly what `extras.shift_locus` recognises, because it asks the
    same parser: `chr1:1000-2000#LTR/Gypsy`, `bedtools getfasta -name`'s
    `TE_1#LTR/Copia::chr1:1000-2000`, and `..` for `-`.
    """
    found = _locate_locus(seq_id)
    if found is None:
        return None
    m = found[0]
    return m.group(1), int(m.group(2)), int(m.group(4))


@dataclass(frozen=True)
class Window:
    """The four reference cuts around one locus, uppercased.

    `up` and `down` are the flanking bases the element itself does not contain;
    `head` and `tail` are its own first and last bases as the reference holds
    them, which is what the record is matched against to settle its orientation.
    All four are clipped at the contig edges, so any of them may be short or
    empty.
    """
    up: str = ""
    head: str = ""
    tail: str = ""
    down: str = ""


def _cuts(start: int, end: int, pad: int, probe: int):
    """The four half-open 0-based intervals `Window` is built from."""
    return ((start - 1 - pad, start - 1), (start - 1, start - 1 + probe),
            (end - probe, end), (end, end + pad))


BLOCK = 1 << 20


def _chunks(path, size: int = BLOCK):
    """Yield `(is_header, bytes)` for a FASTA, in bounded pieces.

    Reading by line would be simpler, but a reference written unwrapped -- one
    line per contig, which `seqkit seq -w 0` and several assemblers emit -- makes
    one line a whole chromosome, and reading it costs twice the contig. Blocks
    are split at newlines where there is one and handed over as fragments where
    there is not, so the peak here is a megabyte whatever the line width.

    Bytes rather than text throughout: decoding 3 Gbp to slice a few kilobases
    out of it is the most expensive thing this module could do, so the extracted
    windows are decoded and nothing else is.
    """
    # By content, not by name: the reference this is handed is routinely a
    # gzipped original renamed `.fa` by the pipeline that staged it.
    op = gzip.open if is_gzip(path) else open
    with op(path, "rb") as fh:
        buf = b""
        while True:
            block = fh.read(size)
            if not block:
                break
            # One split per block, not one slice per line: consuming lines by
            # re-slicing the buffer copies its tail every time, which is
            # quadratic in the block size and dominates everything else.
            parts = (buf + block).split(b"\n")
            buf = parts.pop()
            for line in parts:
                yield line[:1] == b">", line
            # What is left begins a line. A partial header has to wait for its
            # newline; a partial sequence line can go out as a fragment, which
            # is what keeps an unwrapped contig off the heap.
            if buf and buf[:1] != b">":
                yield False, buf
                buf = b""
        if buf:
            yield buf[:1] == b">", buf


def _despace(raw: bytes) -> bytes:
    """Drop every whitespace byte. Returns `raw` itself when there is none.

    Not cosmetic: whitespace inside a sequence line shifts the base offset every
    cut is measured against, so keeping it would silently harvest the wrong
    window rather than fail. The single-part fast path is what stops this from
    copying the whole reference.
    """
    parts = raw.split()
    if len(parts) == 1:
        return parts[0]
    return b"".join(parts)


def _flush(live, out) -> None:
    """Close every cut still open at a contig boundary, keeping what it got."""
    for lo, hi, key, slot, buf in live:
        out[key][slot] = b"".join(buf).upper().decode()


def _scan(path, wants, out, seen) -> None:
    """One streaming pass over one reference file."""
    cuts, live, i, pos = None, [], 0, 0
    for is_header, raw in _chunks(path):
        if is_header:
            _flush(live, out)
            live, i, pos = [], 0, 0
            head = raw[1:].split()
            name = head[0].decode() if head else ""
            if name in seen:
                print(f"Kmer2LTR: warning: reference sequence {name!r} appears "
                      f"more than once; the last copy wins", file=sys.stderr)
            cuts = wants.get(name)
            if cuts is not None:
                seen.add(name)
            continue
        if cuts is None:
            continue
        s = _despace(raw)
        if not s:
            continue
        end = pos + len(s)
        # Cuts are sorted by start, so everything reachable on this line is a
        # prefix of what is left. A cut is opened once and stays open across as
        # many lines as it spans.
        while i < len(cuts) and cuts[i][0] < end:
            lo, hi, key, slot = cuts[i]
            live.append((lo, hi, key, slot, []))
            i += 1
        if live:
            still = []
            for rec in live:
                lo, hi, key, slot, buf = rec
                a, b = max(lo, pos), min(hi, end)
                if b > a:
                    buf.append(s[a - pos:b - pos])
                if hi <= end:
                    out[key][slot] = b"".join(buf).upper().decode()
                else:
                    still.append(rec)
            live = still
        pos = end
    _flush(live, out)


def harvest(paths, loci, pad: int = PAD, probe: int = PROBE) -> dict:
    """`{(chrom, start, end): Window}` for every locus the reference holds.

    `loci` is any iterable of `(chrom, start, end)`; `None` entries and repeats
    are ignored. Loci on contigs no reference file contains are simply absent
    from the result -- a missing contig is a fact about the input, not an error,
    and the caller reports `NA` for those records.
    """
    wants: dict[str, list] = {}
    out: dict[tuple, list[str]] = {}
    for key in loci:
        # An inverted interval is a malformed header, not a locus on the other
        # strand -- strand lives in the sequence here, never in the coordinates.
        if key is None or key in out or key[2] < key[1]:
            continue
        out[key] = ["", "", "", ""]
        for slot, (lo, hi) in enumerate(_cuts(key[1], key[2], pad, probe)):
            lo = max(0, lo)
            if hi > lo:
                wants.setdefault(key[0], []).append((lo, hi, key, slot))
    if not wants:
        return {}
    for v in wants.values():
        v.sort()
    seen: set[str] = set()
    for path in paths:
        _scan(path, wants, out, seen)
    # A locus whose four cuts all fell past the end of its contig is not
    # located, whatever the contig's name said: `orient` would decline it, and
    # counting it would overstate what the caller's summary line reports.
    return {k: Window(*v) for k, v in out.items() if k[0] in seen and any(v)}


# Full IUPAC, both cases: the reference is uppercased on harvest but this is
# also applied to record sequences, which reach it however the caller had them.
_COMPLEMENT = str.maketrans("ACGTURYSWKMBDHVNacgturyswkmbdhvn",
                            "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn")


def revcomp(seq: str) -> str:
    """Reverse complement, leaving unknown characters as themselves."""
    return seq.translate(_COMPLEMENT)[::-1]


def _identity(a: str, b: str) -> float:
    """Fraction of agreeing columns, or -1.0 when too few can be compared.

    Columns where either side is `N` are dropped from both numerator and
    denominator rather than scored: `sanitize` turns every ambiguity code into
    `N`, and a reference assembly gap is `N` on both sides, so counting `N`
    against `N` as agreement would let two gaps anchor a record anywhere.
    """
    n = ok = 0
    for x, y in zip(a, b):
        if x == "N" or y == "N":
            continue
        n += 1
        ok += x == y
    return ok / n if n >= MIN_PROBE else -1.0


@dataclass(frozen=True)
class Context:
    """One record's genomic surroundings, turned to face the record.

    `pad5` abuts the record's 5' terminus and `pad3` its 3' terminus whichever
    way round the record is stored, so everything downstream works in record
    coordinates and never has to think about strand again.

    The two `anchored` flags are per end on purpose. A record whose middle was
    excised still sits exactly on its header's coordinates at both ends, while a
    record that lost bases off one terminus sits on only the other -- and there
    the pad on the bad end abuts nothing, so nothing may be read from it.
    """
    orientation: str
    pad5: str
    pad3: str
    anchored5: bool
    anchored3: bool


def _strand(seq: str, up: str, head: str, tail: str, down: str,
            min_identity: float):
    """`(Context, identity at the `start`-anchored end)` for the better strand.

    Both hypotheses are scored at both ends and the better one wins, so a record
    that lost bases off one terminus is still oriented by the other. Returning
    None rather than guessing is the point: this feeds `--trim-flanks`, and a
    header confidently shifted the wrong way is worse than one left alone.
    """
    p = min(len(head), len(tail), len(seq))
    if p < MIN_PROBE:
        return None
    h, t = head[:p], tail[-p:]
    f5, f3 = _identity(seq[:p], h), _identity(seq[-p:], t)
    r5, r3 = _identity(seq[:p], revcomp(t)), _identity(seq[-p:], revcomp(h))
    fwd, rev = max(f5, f3), max(r5, r3)
    if fwd >= min_identity and fwd > rev:
        return Context("+", up, down, f5 >= min_identity, f3 >= min_identity), f5
    if rev >= min_identity and rev > fwd:
        return (Context("-", revcomp(down), revcomp(up),
                        r5 >= min_identity, r3 >= min_identity), r3)
    return None


def orient(seq: str, window: Window,
           min_identity: float = MIN_IDENTITY) -> Context | None:
    """Which way round `seq` is stored, and its genomic pads. None if undecidable.

    The 5' end is anchored from the header's `start` and the 3' end from its
    `end`, independently and with no reference to the middle: in real annotation
    sets 15-30% of records are shorter than their header span, because a nested
    inner element was excised while the header kept the outer interval.

    **Both coordinate conventions are tried.** `bedtools getfasta -name` writes
    the BED interval verbatim, so its `chr:1000-2000` is 0-based half-open and
    means the 1-based `chr:1001-2000` -- the `end` is the same number and only
    the `start` is one lower. Guessing wrong there does not fail loudly: the
    `end`-anchored terminus still matches, so the record is still oriented, and
    the only symptom is that the `start`-side pad abuts the wrong base and every
    TSD column comes back empty for the whole file. Since `up` ends exactly where
    `head` begins, the two can be re-split to move the start by one base without
    re-reading the reference, and the `start`-anchored identity -- the only
    comparison that can tell the conventions apart -- picks the winner.
    """
    joined = window.up + window.head
    base = len(window.up)
    best = None
    for off in (0, 1):
        # `joined[off : base + off]` is always `base` bases wide: sliding the
        # start by one moves which bases the pad is, not how many.
        got = _strand(seq, joined[off:base + off], joined[base + off:],
                      window.tail, window.down, min_identity)
        if got is not None and (best is None or got[1] > best[1]):
            best = got
    return best[0] if best else None


# --------------------------------------------------------------------------- #
# Target-site duplications
# --------------------------------------------------------------------------- #

# Provisional until the sweep in docs/benchmarks.md sets them.
TSD_K = (6, 5)            # lengths tried, longest first
TSD_SHIFTS = (0, 1, -1)   # boundary offsets tried, smallest first


def _pair(ctx: str, b5: int, b3: int, k: int):
    """The two k-mers flanking `ctx[b5:b3]`, or None if either runs off the end."""
    if b5 - k < 0 or b3 + k > len(ctx) or b5 >= b3:
        return None
    return ctx[b5 - k:b5], ctx[b3:b3 + k]


def _shift_pairs(shifts):
    """Boundary shifts in order of increasing total displacement.

    Sorted rather than nested: `(1, 1)` moves two boundaries and `(-1, 0)` moves
    one, so a lexicographic pass would try the larger claim first. `sorted` is
    stable, so within a displacement the caller's own order decides.
    """
    return sorted(((d5, d3) for d5 in shifts for d3 in shifts),
                  key=lambda p: (abs(p[0]) + abs(p[1]), abs(p[0]), abs(p[1])))


def readable(ctx: str, b5: int, b3: int, ks=None, shifts=None) -> bool:
    """Whether any candidate pair here is real sequence that could be compared.

    This is what separates the output's `.` from its `NA`: `.` says the search
    ran and there is no duplication, `NA` says it could not run. It cannot run
    where the context runs out -- an element at a contig edge -- or where the
    bases are unknown, which covers both a reference assembly gap and the `N`
    that `element_context` writes over the pad of a terminus that failed to
    anchor. Reporting `.` for those would assert a measurement nobody made.
    """
    ks = TSD_K if ks is None else ks
    shifts = TSD_SHIFTS if shifts is None else shifts
    for d5, d3 in _shift_pairs(shifts):
        for k in ks:
            pair = _pair(ctx, b5 + d5, b3 - d3, k)
            if pair and "N" not in pair[0] and "N" not in pair[1]:
                return True
    return False


def find_tsd(ctx: str, b5: int, b3: int, ks=TSD_K, shifts=TSD_SHIFTS):
    """`(motif, d5, d3)` for the target-site duplication flanking `ctx[b5:b3]`.

    A retrotransposon inserts into a staggered cut, so the same few genomic
    bases end up on both sides of it. This looks for that: the k bases
    immediately before the element and the k immediately after, identical.

    `d5` and `d3` shift the respective boundary, positive being *into* the
    element, which is how a boundary called one base wide is recovered.

    The smallest total boundary displacement is tried first and the longest k
    within it, so a reported hit makes the smallest claim that fits. The order matters only for what gets
    reported -- swept over 48 parameter cells it changed the detection rate in
    none of them -- but when a longer duplication at a shifted boundary and a
    shorter one at the boundary as called both fit, saying the boundary is right
    is the smaller claim than saying it is off by a base in both directions.

    Two refusals. A k-mer containing `N` is not evidence of anything. Nor is one
    with fewer than two distinct bases: a homopolymer run matches its own
    reflection at a large fraction of genomic positions, and would swamp the
    signal it is supposed to carry.
    """
    for d5, d3 in _shift_pairs(shifts):
        for k in ks:
            pair = _pair(ctx, b5 + d5, b3 - d3, k)
            if pair is None or pair[0] != pair[1]:
                continue
            left = pair[0]
            if "N" in left or len(set(left)) < 2:
                continue
            return left, d5, d3
    return None


def element_context(seq: str, context: Context) -> tuple[str, int, int]:
    """`(padded, b5, b3)`: the record with its genomic pads, in record coordinates.

    A pad on an end that failed to anchor is replaced by `N`, which is not a
    special case downstream: `find_tsd` already refuses any k-mer containing one,
    so a boundary near that end simply reports no TSD rather than reporting one
    read from the wrong place.
    """
    pad5 = context.pad5 if context.anchored5 else "N" * len(context.pad5)
    pad3 = context.pad3 if context.anchored3 else "N" * len(context.pad3)
    return pad5 + seq + pad3, len(pad5), len(pad5) + len(seq)


# --------------------------------------------------------------------------- #
# Filling the columns
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class Options:
    """What `--genome` was asked to do.

    `anchor` is `--tsd-anchor`, in bits: how much a target-site duplication at
    the record's own termini is worth as evidence against calling a flank. Zero
    -- the default -- means the reference is read but no boundary can move, so
    the TSD stays an external check on the tool's output rather than an input to
    it. See docs/benchmarks.md for why that is the shipped setting.
    """
    anchor: float = 0.0
    ks: tuple[int, ...] = TSD_K
    shifts: tuple[int, ...] = TSD_SHIFTS


def credit(seq: str, context: Context | None, options: Options) -> float:
    """Bits of external evidence that this record's termini are its boundaries.

    Deliberately blind to `options.shifts`. A duplication that appears only once
    a boundary is moved says that boundary is *wrong*; letting it argue against
    trimming would be reading the evidence backwards.
    """
    if context is None or not options.anchor:
        return 0.0
    ctx, b5, b3 = element_context(seq, context)
    return options.anchor if find_tsd(ctx, b5, b3, options.ks, (0,)) else 0.0


def annotate(result, seq: str, context: Context, options: Options):
    """Return `result` with its four genome columns filled in.

    `tsd` is read off the boundary Kmer2LTR settled on, exactly as `motif` is,
    so it describes the row it sits in. `tsd_input` is read off the record's own
    termini, so it describes what was handed in. They coincide whenever no flank
    was called, and separate exactly where the annotator and the tool disagree.

    `orientation` and `tsd_input` are properties of the record rather than of a
    pair, so they are reported even on a row where no pair was found -- there
    they are the only remaining evidence about whether the annotator was
    pointing at a real insertion.

    A duplication that is absent reads `.`; one that could not be looked for
    reads `NA`. `readable` decides which, and the cases it catches are real: an
    element at a contig edge, a reference gap beside it, and a terminus that
    failed to anchor -- where the pad abuts nothing, so reporting "no
    duplication here" would assert a measurement nobody made.

    Motifs are uppercase. `motif` is lowercased because it is a tag compared
    against a constant; a TSD is a sequence, and the sequence it came from is
    uppercase.
    """
    ctx, b5, b3 = element_context(seq, context)

    def look(i: int, j: int):
        """`(motif or '.' or None, offset or None)` for one pair of boundaries."""
        hit = find_tsd(ctx, i, j, options.ks, options.shifts)
        if hit:
            return hit[0], f"{hit[1]},{hit[2]}"
        return ("." if readable(ctx, i, j, options.ks, options.shifts) else None), None

    at_input, _unused = look(b5, b3)
    tsd = offset = None
    if isinstance(result.flank5_len, int) and isinstance(result.flank3_len, int):
        tsd, offset = look(b5 + result.flank5_len, b3 - result.flank3_len)
    return replace(result, orientation=context.orientation, tsd=tsd,
                   tsd_offset=offset, tsd_input=at_input)
