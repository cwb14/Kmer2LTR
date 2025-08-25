#!/usr/bin/env python3
import sys
import re

if len(sys.argv) != 2:
    print(f"Usage: {sys.argv[0]} <input.tsv>", file=sys.stderr)
    sys.exit(1)

path = sys.argv[1]

def parse_start_from_ltr_loc(ltr_loc: str) -> int:
    # Example: "CMHB_chr1:6139054..6147318" or "CMHB_chr1:7618381..7613627"
    m = re.search(r":(\d+)\.\.(\d+)$", ltr_loc)
    if not m:
        raise ValueError(f"Could not parse LTR_loc from: {ltr_loc}")
    coord1, coord2 = int(m.group(1)), int(m.group(2))
    return min(coord1, coord2)

def parse_start_from_internal(internal: str) -> int:
    # Example: "IN:6141817..6144540" -> 6141817
    m = re.search(r"IN:(\d+)\.\.\d+$", internal)
    if not m:
        raise ValueError(f"Could not parse Internal from: {internal}")
    return int(m.group(1))

with open(path, "r", encoding="utf-8") as f:
    header_idx = None
    internal_idx = None

    for line in f:
        line = line.rstrip("\n")
        if not line:
            continue

        # Header line starts with '#'
        if line.startswith("#"):
            hdr = line.lstrip("#")
            cols = hdr.split("\t")
            try:
                header_idx = cols.index("LTR_loc")
                internal_idx = cols.index("Internal")
            except ValueError:
                header_idx = 0 if header_idx is None else header_idx
                internal_idx = 6 if internal_idx is None else internal_idx
            continue

        parts = line.split("\t")
        if header_idx is None or internal_idx is None:
            header_idx, internal_idx = 0, 6

        try:
            ltr_loc = parts[header_idx]
            internal = parts[internal_idx]
            start_ltr = parse_start_from_ltr_loc(ltr_loc)
            start_int = parse_start_from_internal(internal)
            dist = start_int - start_ltr
            print(f"{ltr_loc}\t{dist}")
        except Exception as e:
            print(f"Warning: {e}", file=sys.stderr)
            continue
