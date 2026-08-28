"""Construct ground-truth full-length LTR-RTs from split library entries.

Libraries store many families as separate X-LTR and X-I records. Concatenating
LTR + I + LTR yields an element whose boundaries are known EXACTLY -- the truth
is a construction, so there is no annotation error to confound the benchmark.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ltrk2p.fasta import read_fasta   # noqa: E402

_LTR_SUFFIX = re.compile(r"^(?P<fam>.+?)[-_]LTR(?P<tail>_.*)?$", re.IGNORECASE)
_INT_SUFFIX = re.compile(r"^(?P<fam>.+?)[-_](?:I|INT)(?P<tail>_.*)?$", re.IGNORECASE)
_NON_LTR = re.compile(r"#(DNA|LINE|SINE|RC|Helitron|Satellite|rRNA|tRNA|Simple)", re.I)

# Unambiguous non-LTR classes only. A negative control has to be a class that
# certainly lacks terminal DIRECT repeats, so anything ambiguous is excluded from
# BOTH truth and negatives rather than guessed at -- notably DIRS, Ngaro, Viper
# (tyrosine-recombinase elements with split or inverted terminal repeats),
# Penelope, and Troyka. Excluding them costs a little data and avoids
# contaminating the false-positive estimate with elements that may legitimately
# carry a terminal repeat.
_NEG_CLASSES = re.compile(
    r"^(hAT|Mariner|Tc1|MuDR|Harbinger|EnSpm|CACTA|PiggyBac|P$|Transib|Kolobok|"
    r"Merlin|Sola|Zator|Academ|Ginger|Novosib|Helitron|Crypton|"
    r"L1|L2|CR1|RTE|RTEX|Jockey|R1|R2|I$|LOA|Rex|Proto|Tad1|"
    r"SINE|Alu|B4|MIR|tRNA|5S|7SL|"
    r"DNA transposon|Satellite|Simple|rRNA|snRNA)", re.I)

_AMBIGUOUS_CLASSES = re.compile(r"^(DIRS|Ngaro|Viper|Penelope|Troyka)", re.I)


def family_key(seq_id: str) -> tuple[str, str] | None:
    """Return (family, 'LTR'|'I') for a split-element part, else None."""
    name = seq_id.split("#", 1)[0]
    if _NON_LTR.search(seq_id) or _AMBIGUOUS_CLASSES.search(seq_id):
        return None
    for rx, part in ((_INT_SUFFIX, "I"), (_LTR_SUFFIX, "LTR")):
        m = rx.match(name)
        if m:
            fam = m.group("fam") + (m.group("tail") or "")
            return fam, part
    return None


def pair_library(records):
    """Yield (family, ltr_seq, internal_seq) for families having both parts."""
    parts: dict[str, dict[str, str]] = {}
    order: list[str] = []
    for seq_id, seq in records:
        key = family_key(seq_id)
        if key is None:
            continue
        fam, part = key
        slot = parts.setdefault(fam, {})
        if fam not in order:
            order.append(fam)
        slot.setdefault(part, seq)      # first occurrence wins
    for fam in order:
        slot = parts[fam]
        if "LTR" in slot and "I" in slot:
            yield fam, slot["LTR"], slot["I"]


def build_truth(inputs, out_fasta, out_tsv, min_ltr: int = 80, min_int: int = 200) -> int:
    n = 0
    with open(out_fasta, "w") as fa, open(out_tsv, "w") as tsv:
        tsv.write("elem_id\ttotal_len\tltr_len\tltr5_start\tltr5_end\tltr3_start\tltr3_end\n")
        for path in inputs:
            for fam, ltr, internal in pair_library(read_fasta(path)):
                if len(ltr) < min_ltr or len(internal) < min_int:
                    continue
                if ltr.count("N") > 0.05 * len(ltr) or internal.count("N") > 0.05 * len(internal):
                    continue
                elem = ltr + internal + ltr
                eid = f"{Path(path).stem}|{fam}"
                fa.write(f">{eid}\n{elem}\n")
                tsv.write(f"{eid}\t{len(elem)}\t{len(ltr)}\t1\t{len(ltr)}"
                          f"\t{len(ltr) + len(internal) + 1}\t{len(elem)}\n")
                n += 1
    return n


def _records_with_header(path):
    """Yield (full_header_line, sanitized_seq).

    repbase encodes its class in tab-delimited field 2 of the header
    (`>MARINERN10_AG\\tMariner/Tc1\\tAnopheles gambiae`), which read_fasta
    discards when it truncates the id at the first whitespace. Negative-control
    selection needs that field, so this reader keeps the whole header. Streaming,
    one record at a time -- repbase.fa is 369 MB.
    """
    import gzip
    from ltrk2p.fasta import sanitize
    opener = gzip.open if str(path).endswith(".gz") else open
    header, chunks = None, []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, sanitize("".join(chunks))
                header, chunks = line[1:].rstrip("\n"), []
            elif header is not None:
                chunks.append(line)
    if header is not None:
        yield header, sanitize("".join(chunks))


def build_negatives(inputs, out_fasta, min_len: int = 500) -> int:
    """Write non-LTR entries >= min_len as negative controls.

    A record counts as a negative if its class is recognisable AND
    unambiguously non-LTR, from either header convention:
      - Dfam/MTEC/riceTElib ``name#class`` form -- checked with the existing
        _NON_LTR regex against the truncated id (same id read_fasta yields).
      - repbase's tab-delimited ``id<TAB>class<TAB>species`` form -- checked
        with _NEG_CLASSES against field 2, which read_fasta's id-only
        contract discards (id is truncated at the first whitespace, which
        includes tabs).
    Records whose class matches _AMBIGUOUS_CLASSES (DIRS/Ngaro/Viper/
    Penelope/Troyka -- tyrosine-recombinase elements that may carry a
    terminal repeat) are skipped entirely, in either form, rather than
    guessed at -- the same guard family_key applies to keep them out of
    the truth set.
    """
    n = 0
    with open(out_fasta, "w") as fa:
        for path in inputs:
            for header, seq in _records_with_header(path):
                if len(seq) < min_len:
                    continue
                stripped = header.strip()
                id_part = stripped.split(None, 1)[0] if stripped else ""
                fields = header.split("\t")
                field2 = fields[1].strip() if len(fields) > 1 else ""

                if _AMBIGUOUS_CLASSES.search(id_part) or (
                    field2 and _AMBIGUOUS_CLASSES.search(field2)
                ):
                    continue

                is_neg = bool(_NON_LTR.search(id_part)) or bool(
                    field2 and _NEG_CLASSES.search(field2)
                )
                if not is_neg:
                    continue

                fa.write(f">{Path(path).stem}|{id_part}\n{seq}\n")
                n += 1
    return n


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Build ground-truth and negative sets.")
    ap.add_argument("inputs", nargs="+")
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-negatives")
    a = ap.parse_args()
    print(f"truth elements: {build_truth(a.inputs, a.out_fasta, a.out_tsv)}")
    if a.out_negatives:
        print(f"negatives: {build_negatives(a.inputs, a.out_negatives)}")
