#!/usr/bin/env python3

import sys
import argparse
from Bio import SeqIO

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Map kmers to positions in FASTA and apply distance and positional filters.")
    parser.add_argument("kmer_file", help="Kmer count file from Jellyfish (e.g., kmer_twice.txt)")
    parser.add_argument("fasta_file", help="FASTA file to search (e.g., LTRRT.fa)")
    parser.add_argument("-d", "--min_distance", type=int, default=0,
                        help="Minimum distance (bp) required between two kmer matches (default: 0 = no filter)")
    return parser.parse_args()

def main():
    args = parse_arguments()
    kmer_dict = {}

    # Read kmers with exactly 2 counts
    with open(args.kmer_file) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                kmer, count = parts
                if count == "2":
                    kmer_dict[kmer.upper()] = []

    # Read FASTA (assuming only 1 sequence)
    fasta_records = list(SeqIO.parse(args.fasta_file, "fasta"))
    if len(fasta_records) != 1:
        print("Error: Expected exactly one sequence in FASTA file.", file=sys.stderr)
        sys.exit(1)

    record = fasta_records[0]
    seq = str(record.seq).upper()
    seq_len = len(seq)
    half_len = seq_len // 2

    # Find kmers
    for kmer in kmer_dict:
        k = len(kmer)
        start = 0
        while True:
            pos = seq.find(kmer, start)
            if pos == -1:
                break
            kmer_dict[kmer].append(pos)
            start = pos + 1

    # Output
    print("kmer\tcount\tstart1\tend1\tstart2\tend2")
    for kmer, positions in kmer_dict.items():
        if len(positions) == 2:
            pos1, pos2 = sorted(positions)
            distance = abs(pos2 - pos1)
            start1 = pos1 + 1
            start2 = pos2 + 1
            end1 = pos1 + len(kmer)
            end2 = pos2 + len(kmer)

            # Filter by min_distance
            if distance < args.min_distance:
                continue

            # Filter by half-sequence mapping
            if start1 > half_len or start2 <= half_len:
                continue

            print(f"{kmer}\t2\t{start1}\t{end1}\t{start2}\t{end2}")

if __name__ == "__main__":
    main()
