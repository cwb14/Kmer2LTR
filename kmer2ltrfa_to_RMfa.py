#!/usr/bin/env python3
"""
kmer2ltrfa_to_RMfa.py

Example:
python kmer2ltrfa_to_RMfa.py \
  --fasta Athal.fa_Kmer2LTR_TSD_class_purge.fa \
  --tsv   Athal.fa_Kmer2LTR_TSD_class_purge \
  --clean-non-bases \
  --unknown-class LTR \
  --filterin-class LTR \
  --make-perfect-repeat \
  > Athal.fa_Kmer2LTR_TSD_class_purge.clean.fa

Convert a FASTA + TSV annotation into a RepeatMasker-style FASTA.

- FASTA headers: >NAME
- TSV (tab-delimited, NO header):
    col1  = NAME (matches FASTA header, without '>')
    col2  = LEN1
    col3  = LEN2
    col15 = Class
    col16 = Superfamily

Output headers:
    >NAME#Class/Superfamily

Optional:
- Clean non-ACGT bases
- Handle NA values for class/superfamily
- Filter sequences by class (filter stage before/after NA conversion)
- --make-perfect-repeat:
    LONGER_LEN = max(LEN1, LEN2)
    Replace the rightmost (3') LONGER_LEN bases with the leftmost (5') LONGER_LEN bases.
"""

import sys
import argparse
import csv
import re
from typing import Dict, Tuple, Optional, Set


# name -> (raw_class, raw_superfamily, final_class, final_superfamily, len1, len2)
AnnTuple = Tuple[str, str, str, str, Optional[int], Optional[int]]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="kmer2ltrfa_to_RMfa.py",
        description=(
            "Convert a FASTA and a TSV annotation table to a RepeatMasker-style FASTA.\n\n"
            "TSV is tab-delimited with NO header:\n"
            "  - Column 1  (index 1): Name (matches FASTA header without '>')\n"
            "  - Column 2  (index 2): LEN1\n"
            "  - Column 3  (index 3): LEN2\n"
            "  - Column 15 (index 15): Class\n"
            "  - Column 16 (index 16): Superfamily\n\n"
            "Output headers:\n"
            "  >NAME#Class/Superfamily\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument("--fasta", required=True, help="Input FASTA file.")
    parser.add_argument("--tsv", required=True, help="Input TSV file (no header).")

    parser.add_argument(
        "--clean-non-bases",
        action="store_true",
        help="Convert any character not in A,C,G,T (case-insensitive) to 'T' and uppercase the sequence.",
    )

    parser.add_argument(
        "--unknown-class",
        default="Unknown",
        help="Label to use when Class (col15) is 'NA' or empty. Default: 'Unknown'.",
    )
    parser.add_argument(
        "--unknown-superfamily",
        default="Unknown",
        help="Label to use when Superfamily (col16) is 'NA' or empty. Default: 'Unknown'.",
    )

    parser.add_argument(
        "--filterin-class",
        default=None,
        help="Comma-separated list of class names to RETAIN (case-insensitive).",
    )
    parser.add_argument(
        "--filterout-class",
        default=None,
        help="Comma-separated list of class names to OMIT (case-insensitive).",
    )
    parser.add_argument(
        "--filter-stage",
        choices=["before", "after"],
        default="after",
        help="Filter on raw class (before) or NA-converted class (after). Default: after.",
    )

    parser.add_argument(
        "--make-perfect-repeat",
        action="store_true",
        help=(
            "Use Kmer2LTR repeat classification to build a perfect (unmutated) TE. "
            "Replace the rightmost repeat (3') with the leftmost repeat (5')."
        ),
    )

    args = parser.parse_args()

    if args.filterin_class and args.filterout_class:
        parser.error("You cannot use --filterin-class and --filterout-class at the same time.")

    return args


def read_tsv_annotations(
    tsv_path: str,
    unknown_class: str,
    unknown_superfamily: str,
) -> Dict[str, AnnTuple]:
    """
    Read TSV and return a mapping:
      name -> (raw_class, raw_superfamily, final_class, final_superfamily, len1, len2)

    TSV is NO-header, tab-delimited.
    """
    annotations: Dict[str, AnnTuple] = {}

    def is_na(s: str) -> bool:
        return s == "" or s.upper() == "NA"

    def parse_int_maybe(s: str) -> Optional[int]:
        s = s.strip()
        if s == "" or s.upper() == "NA":
            return None
        # allow "123.0" or "123" (but not arbitrary text)
        try:
            return int(float(s))
        except ValueError:
            return None

    with open(tsv_path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for lineno, row in enumerate(reader, start=1):
            if not row:
                continue

            # Need at least 16 columns for class/superfamily, plus cols 2/3 for LEN1/LEN2.
            if len(row) < 16:
                sys.stderr.write(
                    f"WARNING: Skipping TSV line {lineno} (expected >=16 columns, got {len(row)}):\n"
                    f"  {row}\n"
                )
                continue

            name = row[0].strip()
            len1 = parse_int_maybe(row[1])  # col2 (1-based)
            len2 = parse_int_maybe(row[2])  # col3 (1-based)

            raw_class = row[14].strip()  # col15 (1-based)
            raw_super = row[15].strip()  # col16 (1-based)

            final_class = raw_class if not is_na(raw_class) else unknown_class
            final_super = raw_super if not is_na(raw_super) else unknown_superfamily

            if name in annotations:
                sys.stderr.write(
                    f"WARNING: Duplicate annotation for '{name}' in TSV (line {lineno}); "
                    f"keeping the first occurrence.\n"
                )
                continue

            annotations[name] = (raw_class, raw_super, final_class, final_super, len1, len2)

    return annotations


def fasta_reader(fasta_path: str):
    """Yield (header_without_gt, sequence_string)."""
    header: Optional[str] = None
    seq_chunks = []

    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())

    if header is not None:
        yield header, "".join(seq_chunks)


def clean_sequence(seq: str) -> str:
    """Uppercase and replace any non-ACGT with 'T'."""
    seq = seq.upper()
    return re.sub(r"[^ACGT]", "T", seq)


def parse_class_list(arg: Optional[str]) -> Optional[Set[str]]:
    if arg is None:
        return None
    items = [x.strip() for x in arg.split(",") if x.strip() != ""]
    return {x.upper() for x in items}


def normalize_class_for_compare(cls: Optional[str]) -> str:
    if cls is None:
        return ""
    return cls.strip().upper()


def should_keep_record(
    raw_class: str,
    final_class: str,
    filterin: Optional[Set[str]],
    filterout: Optional[Set[str]],
    filter_stage: str,
) -> bool:
    if filterin is None and filterout is None:
        return True

    cls_for_filter = normalize_class_for_compare(raw_class if filter_stage == "before" else final_class)

    if filterin is not None:
        return cls_for_filter in filterin
    if filterout is not None:
        return cls_for_filter not in filterout
    return True


def make_perfect_repeat(seq: str, longer_len: int) -> str:
    """
    Replace the last longer_len bases with the first longer_len bases.
    If longer_len <= 0: return seq unchanged.
    If longer_len > len(seq): clamp to len(seq).
    """
    if longer_len <= 0:
        return seq
    L = len(seq)
    if L == 0:
        return seq
    if longer_len > L:
        longer_len = L
    left = seq[:longer_len]
    return seq[:-longer_len] + left


def main():
    args = parse_args()

    ann = read_tsv_annotations(
        args.tsv,
        unknown_class=args.unknown_class,
        unknown_superfamily=args.unknown_superfamily,
    )

    filterin_set = parse_class_list(args.filterin_class)
    filterout_set = parse_class_list(args.filterout_class)

    n_total = 0
    n_no_annotation = 0
    n_filtered_out = 0
    n_written = 0
    n_perfect_applied = 0
    n_perfect_skipped_no_len = 0

    for name, seq in fasta_reader(args.fasta):
        n_total += 1

        if name not in ann:
            sys.stderr.write(
                f"WARNING: No TSV annotation found for '{name}'; "
                f"using unknown-class '{args.unknown_class}' and "
                f"unknown-superfamily '{args.unknown_superfamily}'.\n"
            )
            raw_class = ""
            raw_super = ""
            final_class = args.unknown_class
            final_super = args.unknown_superfamily
            len1 = None
            len2 = None
            n_no_annotation += 1
        else:
            raw_class, raw_super, final_class, final_super, len1, len2 = ann[name]

        if not should_keep_record(
            raw_class=raw_class,
            final_class=final_class,
            filterin=filterin_set,
            filterout=filterout_set,
            filter_stage=args.filter_stage,
        ):
            n_filtered_out += 1
            continue

        if args.clean_non_bases:
            seq = clean_sequence(seq)

        if args.make_perfect_repeat:
            if len1 is None and len2 is None:
                n_perfect_skipped_no_len += 1
                sys.stderr.write(
                    f"WARNING: --make-perfect-repeat set but LEN1/LEN2 missing/unparseable for '{name}'; leaving sequence unchanged.\n"
                )
            else:
                l1 = len1 if len1 is not None else -1
                l2 = len2 if len2 is not None else -1
                longer_len = max(l1, l2)
                if longer_len < 0:
                    n_perfect_skipped_no_len += 1
                    sys.stderr.write(
                        f"WARNING: --make-perfect-repeat set but LEN1/LEN2 invalid for '{name}'; leaving sequence unchanged.\n"
                    )
                else:
                    which = "LEN1" if l1 >= l2 else "LEN2"
#                    sys.stderr.write(
#                        f"INFO: {name}: LEN1={len1} LEN2={len2} -> LONGER_LEN={longer_len} ({which})\n"
#                    )
                    seq = make_perfect_repeat(seq, longer_len)
                    n_perfect_applied += 1

        cls_header = final_class.strip().replace(" ", "_")
        super_header = final_super.strip().replace(" ", "_")

        header = f"{name}#{cls_header}/{super_header}"

        sys.stdout.write(f">{header}\n")
        sys.stdout.write(f"{seq}\n")
        n_written += 1

    sys.stderr.write(
        "Done.\n"
        f"  Total FASTA records:        {n_total}\n"
        f"  Records with no annotation: {n_no_annotation}\n"
        f"  Records filtered out:       {n_filtered_out}\n"
        f"  Records written:            {n_written}\n"
        f"  Perfect-repeat applied:     {n_perfect_applied}\n"
        f"  Perfect-repeat skipped:     {n_perfect_skipped_no_len}\n"
    )


if __name__ == "__main__":
    main()
