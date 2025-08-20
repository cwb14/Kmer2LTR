#!/usr/bin/env python3
"""
extract_ltrs_adjust.py

Usage:
    python extract_ltrs_adjust.py kmer.filtered input.fa output.fa [-e EXTEND]

Arguments:
    kmer.filtered    : TSV file with columns kmer, count, start1, end1, start2, end2
    input.fa         : FASTA file containing a single sequence
    output.fa        : Output FASTA with 5' and 3' LTRs
    -e, --extend N   : Optional extension (default 0) to extend both LTRs by N bp

Example:
    python extract_ltrs_adjust.py kmer.filtered LTRRT.fa LTRs.fa -e 10
"""

import sys
import argparse

def parse_kmer(kmer_tsv):
    max_end1 = 0
    min_start2 = float('inf')
    with open(kmer_tsv) as f:
        next(f)  # skip header
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 6:
                continue
            end1 = int(cols[3])
            start2 = int(cols[4])
            max_end1 = max(max_end1, end1)
            min_start2 = min(min_start2, start2)
    return max_end1, min_start2

def parse_fasta(fasta_file):
    with open(fasta_file) as f:
        header_line = f.readline().strip()
        if not header_line.startswith('>'):
            sys.exit("ERROR: FASTA does not begin with '>'")
        header = header_line[1:]
        seq = ''.join(line.strip() for line in f if not line.startswith('>'))
    return header, seq

def main():
    parser = argparse.ArgumentParser(description="Extract adjusted and equal-length 5' and 3' LTRs.")
    parser.add_argument("kmer_tsv", help="Input kmer.filtered TSV file")
    parser.add_argument("fasta", help="Input FASTA file (single sequence)")
    parser.add_argument("output", help="Output FASTA file for LTRs")
    parser.add_argument("-e", "--extend", type=int, default=0, help="Number of bases to extend both LTRs (default: 0)")
    args = parser.parse_args()

    max_end1, min_start2 = parse_kmer(args.kmer_tsv)
    header, seq = parse_fasta(args.fasta)
    seq_len = len(seq)

    # Unadjusted lengths (1-based inclusive coordinates from TSV)
    len5p = max_end1
    len3p = seq_len - (min_start2 - 1)

    # Equalize LTR lengths by expanding the shorter inward (toward the center)
    new_end1 = max_end1
    new_start2 = min_start2

    if len5p < len3p:
        diff = len3p - len5p
        new_end1 += diff
        if new_end1 > seq_len:
            new_end1 = seq_len
    elif len3p < len5p:
        diff = len5p - len3p
        new_start2 -= diff
        if new_start2 < 1:
            new_start2 = 1

    # Apply symmetric extension
    new_end1 += args.extend
    if new_end1 > seq_len:
        new_end1 = seq_len

    new_start2 -= args.extend
    if new_start2 < 1:
        new_start2 = 1

    # ---- SAFETY CHECK: prevent overlap by enforcing a maximum equal length ----
    # Proposed lengths after extension/clipping
    prop_len_5p = new_end1  # 1..new_end1  (inclusive)
    prop_len_3p = seq_len - new_start2 + 1  # new_start2..seq_len (inclusive)

    # Absolute maximum to avoid overlap is half the sequence length (floor)
    max_len_no_overlap = seq_len // 2

    # Final equal LTR length we can support without overlap and without exceeding proposed lengths
    final_len = min(max_len_no_overlap, prop_len_5p, prop_len_3p)

    # If we had to snip anything down, warn to stderr
    if final_len < prop_len_5p or final_len < prop_len_3p:
        print(
            f"WARNING: LTRs would overlap after extension. "
            f"Snipping to {final_len} bp each (max allowed: {max_len_no_overlap}).",
            file=sys.stderr
        )

    # Recompute coordinates symmetrically to guarantee equal length and no overlap
    # 5' LTR: 1..final_len
    # 3' LTR: (seq_len - final_len + 1)..seq_len
    new_end1 = final_len
    new_start2 = seq_len - final_len + 1

    # Extract LTR sequences
    ltr_5p = seq[0:new_end1]  # Python slice end is exclusive; new_end1 is 1-based inclusive
    ltr_3p = seq[new_start2 - 1:]  # start index is 0-based

    # Write output
    with open(args.output, 'w') as out:
        out.write(f">{header}\t5p:1-{new_end1}\n")
        out.write(f"{ltr_5p}\n")
        out.write(f">{header}\t3p:{new_start2}-{seq_len}\n")
        out.write(f"{ltr_3p}\n")

if __name__ == "__main__":
    main()
