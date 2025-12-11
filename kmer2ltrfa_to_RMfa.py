#!/usr/bin/env python3
"""
kmer2ltrfa_to_RMfa.py

python kmer2ltrfa_to_RMfa.py --fasta Athal.fa_Kmer2LTR_TSD_class_purge.fa --tsv Athal.fa_Kmer2LTR_TSD_class_purge --clean-non-bases --unknown-class LTR  --filterin-class LTR  > Athal.fa_Kmer2LTR_TSD_class_purge.clean.fa

Convert a FASTA + TSV annotation into a RepeatMasker-style FASTA.

- FASTA headers: >NAME
- TSV columns:
    col1  = NAME (matches FASTA header, without '>')
    col15 = Class
    col16 = Superfamily

Output headers:
    >NAME#Class/Superfamily

With options for:
- Cleaning non-ACGT bases
- Handling NA values for class/superfamily
- Filtering sequences by class, with control over whether filtering
  uses the raw class values or the post-NA-conversion class values.
"""

import sys
import argparse
import csv
import re
from typing import Dict, Tuple, Optional, Set


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="kmer2ltrfa_to_RMfa.py",
        description=(
            "Convert a FASTA and a TSV annotation table to a RepeatMasker-style FASTA.\n\n"
            "Each FASTA header should be a name like:\n"
            "  >Athal_chr2:3945512-3959397\n\n"
            "The TSV should be tab-delimited, with:\n"
            "  - Column 1  (index 1): Name (matches FASTA header without '>')\n"
            "  - Column 15 (index 15): Class\n"
            "  - Column 16 (index 16): Superfamily\n\n"
            "The output FASTA headers will look like:\n"
            "  >NAME#Class/Superfamily\n\n"
            "Class and Superfamily are taken from the TSV annotation, with\n"
            "optional handling for NA values, and optional filtering by Class."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "EXAMPLES\n"
            "========\n\n"
            "1) Basic usage (default unknown labels):\n"
            "   python kmer2ltrfa_to_RMfa.py \\\n"
            "       --fasta multiseq.fa \\\n"
            "       --tsv file.tsv \\\n"
            "       > output.fa\n\n"
            "   If a row has NA in Class and Superfamily (col15/col16), the\n"
            "   header will be, e.g.:\n"
            "   >Athal_chr5:8032695-8034748#Unknown/Unknown\n\n"
            "2) Custom labels for NA class/superfamily:\n"
            "   python kmer2ltrfa_to_RMfa.py \\\n"
            "       --fasta multiseq.fa \\\n"
            "       --tsv file.tsv \\\n"
            "       --unknown-class fun \\\n"
            "       --unknown-superfamily boring \\\n"
            "       > output.fa\n\n"
            "   NA class  -> 'fun'\n"
            "   NA superfamily -> 'boring'\n"
            "   Header example:\n"
            "   >Athal_chr5:8032695-8034748#fun/boring\n\n"
            "3) Clean non-ACGT bases (e.g. N -> T, X -> T):\n"
            "   python kmer2ltrfa_to_RMfa.py \\\n"
            "       --fasta multiseq.fa \\\n"
            "       --tsv file.tsv \\\n"
            "       --clean-non-bases \\\n"
            "       > output.fa\n\n"
            "4) Keep only specific classes (filter-in):\n"
            "   python kmer2ltrfa_to_RMfa.py \\\n"
            "       --fasta multiseq.fa \\\n"
            "       --tsv file.tsv \\\n"
            "       --filterin-class LTR,TIR \\\n"
            "       > output.fa\n\n"
            "   This retains only sequences where Class is LTR or TIR.\n\n"
            "5) Drop specific classes (filter-out):\n"
            "   python kmer2ltrfa_to_RMfa.py \\\n"
            "       --fasta multiseq.fa \\\n"
            "       --tsv file.tsv \\\n"
            "       --filterout-class LTR,TIR \\\n"
            "       > output.fa\n\n"
            "   This retains only sequences whose Class is NOT LTR or TIR.\n\n"
            "6) Control when filtering happens relative to NA handling:\n"
            "   By default, filtering uses the Class AFTER NA conversion\n"
            "   ('Unknown' etc.). To filter on the raw TSV values instead,\n"
            "   use:\n\n"
            "   python kmer2ltrfa_to_RMfa.py \\\n"
            "       --fasta multiseq.fa \\\n"
            "       --tsv file.tsv \\\n"
            "       --filterin-class NA \\\n"
            "       --filter-stage before \\\n"
            "       > only_NA.fa\n"
        )
    )

    parser.add_argument(
        "--fasta",
        required=True,
        help="Input FASTA file with headers like '>Athal_chr2:3945512-3959397'."
    )
    parser.add_argument(
        "--tsv",
        required=True,
        help="Input TSV file with name+annotation (name in column 1, class in column 15, superfamily in column 16)."
    )
    parser.add_argument(
        "--clean-non-bases",
        action="store_true",
        help=(
            "Convert any character not in A,C,G,T (case-insensitive) to 'T'. "
            "The entire sequence is also converted to uppercase. "
            "For example: ATCGGXTCGC -> ATCGGTTCGC."
        )
    )
    parser.add_argument(
        "--unknown-class",
        default="Unknown",
        help=(
            "Label to use when Class (col 15) is 'NA' or empty. "
            "Default: 'Unknown'."
        )
    )
    parser.add_argument(
        "--unknown-superfamily",
        default="Unknown",
        help=(
            "Label to use when Superfamily (col 16) is 'NA' or empty. "
            "Default: 'Unknown'."
        )
    )
    parser.add_argument(
        "--filterin-class",
        default=None,
        help=(
            "Comma-separated list of class names to RETAIN. "
            "Filtering is case-insensitive. "
            "Example: --filterin-class LTR,TIR\n"
            "NOTE: You cannot use --filterin-class and --filterout-class together."
        )
    )
    parser.add_argument(
        "--filterout-class",
        default=None,
        help=(
            "Comma-separated list of class names to OMIT. "
            "Filtering is case-insensitive. "
            "Example: --filterout-class LTR,TIR\n"
            "NOTE: You cannot use --filterin-class and --filterout-class together."
        )
    )
    parser.add_argument(
        "--filter-stage",
        choices=["before", "after"],
        default="after",
        help=(
            "When to apply class-based filtering relative to NA conversion.\n"
            "  before: filter on the RAW class from the TSV (e.g. 'NA', 'LTR').\n"
            "  after : filter on the class AFTER NA conversion (e.g. 'Unknown', 'LTR').\n"
            "Default: after."
        )
    )

    args = parser.parse_args()

    # Basic sanity check: cannot combine filterin and filterout
    if args.filterin_class and args.filterout_class:
        parser.error(
            "You cannot use --filterin-class and --filterout-class at the same time.\n"
            "Please choose one type of filtering."
        )

    return args


def read_tsv_annotations(
    tsv_path: str,
    unknown_class: str,
    unknown_superfamily: str
) -> Dict[str, Tuple[str, str, str, str]]:
    """
    Read TSV and return a mapping:
        name -> (raw_class, raw_superfamily, final_class, final_superfamily)

    raw_class / raw_superfamily come directly from TSV (stripped).
    final_class / final_superfamily have NA/empty converted to unknown_*.
    """
    annotations: Dict[str, Tuple[str, str, str, str]] = {}

    with open(tsv_path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for lineno, row in enumerate(reader, start=1):
            if not row:
                continue

            # Require at least 16 columns (1-based index 16 for superfamily)
            if len(row) < 16:
                sys.stderr.write(
                    f"WARNING: Skipping TSV line {lineno} (expected >=16 columns, got {len(row)}):\n"
                    f"  {row}\n"
                )
                continue

            name = row[0].strip()
            raw_class = row[14].strip()  # column 15 (1-based)
            raw_super = row[15].strip()  # column 16 (1-based)

            # NA/empty handling
            def is_na(s: str) -> bool:
                return s == "" or s.upper() == "NA"

            final_class = raw_class if not is_na(raw_class) else unknown_class
            final_super = raw_super if not is_na(raw_super) else unknown_superfamily

            if name in annotations:
                sys.stderr.write(
                    f"WARNING: Duplicate annotation for '{name}' in TSV (line {lineno}); "
                    f"keeping the first occurrence.\n"
                )
                continue

            annotations[name] = (raw_class, raw_super, final_class, final_super)

    return annotations


def fasta_reader(fasta_path: str):
    """
    Simple FASTA reader yielding (header_without_gt, sequence_string).

    Header is returned WITHOUT leading '>'.
    Sequences are concatenated across multiple lines.
    """
    header: Optional[str] = None
    seq_chunks = []

    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                # yield previous
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())

    # last record
    if header is not None:
        yield header, "".join(seq_chunks)


def clean_sequence(seq: str) -> str:
    """
    Convert sequence to uppercase and replace any non-ACGT with 'T'.
    """
    seq = seq.upper()
    # Anything not A/C/G/T becomes T
    return re.sub(r"[^ACGT]", "T", seq)


def parse_class_list(arg: Optional[str]) -> Optional[Set[str]]:
    """
    Parse a comma-separated class list into a set of UPPERCASE strings.
    """
    if arg is None:
        return None
    items = [x.strip() for x in arg.split(",") if x.strip() != ""]
    return {x.upper() for x in items}


def normalize_class_for_compare(cls: Optional[str]) -> str:
    """
    Normalize a class string for comparison: strip, uppercase.
    """
    if cls is None:
        return ""
    return cls.strip().upper()


def should_keep_record(
    raw_class: str,
    final_class: str,
    filterin: Optional[Set[str]],
    filterout: Optional[Set[str]],
    filter_stage: str
) -> bool:
    """
    Decide whether to keep a record based on filterin/filterout and filter_stage.

    filter_stage:
        'before' -> use raw_class
        'after'  -> use final_class
    """
    if filterin is None and filterout is None:
        return True  # no filtering

    if filter_stage == "before":
        cls_for_filter = normalize_class_for_compare(raw_class)
    else:
        cls_for_filter = normalize_class_for_compare(final_class)

    if filterin is not None:
        # Keep only if class is in filterin set
        return cls_for_filter in filterin

    if filterout is not None:
        # Drop if class is in filterout set
        return cls_for_filter not in filterout

    # Should not reach here
    return True


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

    for name, seq in fasta_reader(args.fasta):
        n_total += 1

        if name not in ann:
            # No annotation: use unknown labels directly
            sys.stderr.write(
                f"WARNING: No TSV annotation found for '{name}'; "
                f"using unknown-class '{args.unknown_class}' and "
                f"unknown-superfamily '{args.unknown_superfamily}'.\n"
            )
            raw_class = ""
            raw_super = ""
            final_class = args.unknown_class
            final_super = args.unknown_superfamily
            n_no_annotation += 1
        else:
            raw_class, raw_super, final_class, final_super = ann[name]

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

        # Sanitize class and superfamily for header (avoid spaces)
        cls_header = final_class.strip().replace(" ", "_")
        super_header = final_super.strip().replace(" ", "_")

        # Build RepeatMasker header
        header = f"{name}#{cls_header}/{super_header}"

        # Write FASTA entry
        sys.stdout.write(f">{header}\n")
        # print sequence as-is on a single line; user can wrap later if desired
        sys.stdout.write(f"{seq}\n")
        n_written += 1

    # Stats to stderr
    sys.stderr.write(
        f"Done.\n"
        f"  Total FASTA records:        {n_total}\n"
        f"  Records with no annotation: {n_no_annotation}\n"
        f"  Records filtered out:       {n_filtered_out}\n"
        f"  Records written:            {n_written}\n"
    )


if __name__ == "__main__":
    main()
