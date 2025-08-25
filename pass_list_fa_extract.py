#!/usr/bin/env python3
import sys
import argparse
import gzip
import re

def open_maybe_gzip(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode=mode)
    return open(path, mode=mode)

def read_fasta(path):
    """Load a FASTA into a dict {header: sequence} (header = first token after '>')."""
    seqs = {}
    with open_maybe_gzip(path, "rt") as fh:
        name = None
        chunks = []
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).replace(" ", "").replace("\n", "")
                name = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            seqs[name] = "".join(chunks).replace(" ", "").replace("\n", "")
    return seqs

coord_re = re.compile(r"^(?P<chrom>[^:]+):(?P<start>\d+)\.\.(?P<end>\d+)$")

def parse_coords_from_tsv(tsv_path):
    """
    Yield (chrom, start, end, original_label) for each non-header line.
    Assumes column 1 looks like 'chrom:start..end'. Coordinates are 1-based inclusive.
    """
    with open_maybe_gzip(tsv_path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if not fields:
                continue
            label = fields[0].strip()
            m = coord_re.match(label)
            if not m:
                print(f"WARNING: couldn't parse coord in first column: {label}", file=sys.stderr)
                continue
            chrom = m.group("chrom")
            start = int(m.group("start"))
            end = int(m.group("end"))
            # normalize if start > end (rare; just swap)
            if start > end:
                start, end = end, start
            yield chrom, start, end, label

def main():
    ap = argparse.ArgumentParser(description="Extract sequences from FASTA using coords in TSV col1 'chrom:start..end'.")
    ap.add_argument("-fa", "--fasta", required=True, help="Genome FASTA (optionally .gz)")
    ap.add_argument("-tsv", "--tsv", required=True, help="TSV with coords in column 1 (header lines start with '#'; optionally .gz)")
    args = ap.parse_args()

    genome = read_fasta(args.fasta)
    if not genome:
        print("ERROR: No sequences loaded from FASTA.", file=sys.stderr)
        sys.exit(1)

    out = sys.stdout
    n_ok = 0
    n_skip = 0

    for chrom, start, end, label in parse_coords_from_tsv(args.tsv):
        if chrom not in genome:
            print(f"WARNING: chromosome '{chrom}' not found in FASTA for {label}; skipping.", file=sys.stderr)
            n_skip += 1
            continue
        seq = genome[chrom]
        # bounds check (1-based inclusive -> Python slice [start-1:end])
        if start < 1 or end > len(seq):
            print(f"WARNING: out-of-bounds for {label} (len={len(seq)}); skipping.", file=sys.stderr)
            n_skip += 1
            continue
        subseq = seq[start-1:end]
        # Output exactly in requested format:
        out.write(f">{label}\n")
        out.write(f"{subseq}\n")
        n_ok += 1

    print(f"# Done. Extracted {n_ok} regions; skipped {n_skip}.", file=sys.stderr)

if __name__ == "__main__":
    main()
