#!/usr/bin/env python3
"""
LTR_mutator.py

Mutates only the LTR regions of intact LTR-RT sequences.
Usage:
    python LTR_mutator.py \
        -i LTR_RT_lib.fa \
        -o LTR_RT_lib_mutated.fa \
        -mp 15 \
        -TiTv 2.0 \
        [--seed 42]
"""

import argparse
import re
import random
import sys
from Bio import SeqIO
from Bio.Seq import Seq

def parse_args():
    parser = argparse.ArgumentParser(
        description="Mutate LTR regions of intact LTR-RT sequences"
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Input FASTA file of intact LTR-RTs"
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output FASTA file for mutated sequences"
    )
    parser.add_argument(
        "-mp", "--mutation_percent", type=float, required=True,
        help="Percent of LTR bases to mutate (e.g. 15 for 15%%)"
    )
    parser.add_argument(
        "-TiTv", "--titv_ratio", type=float, default=1.0,
        help="Transition/transversion ratio (e.g. 2.0 for twice as many transitions)"
    )
    parser.add_argument(
        "--seed", type=int, default=None,
        help="Random seed for reproducibility (optional)"
    )
    return parser.parse_args()

def mutate_ltr(seq_str, ltr_len, mp, titv_ratio):
    """
    Mutate the two LTRs (first ltr_len bases and last ltr_len bases)
    by sampling mutation events WITH replacement.
    Returns (mutated_seq_str, measured_percent_unique_sites).
    """
    seq = list(seq_str)
    seq_len = len(seq)
    total_ltr_bases = 2 * ltr_len
    # Number of mutation *events*
    n_events = int(round(total_ltr_bases * mp / 100.0))
    
    # Define transition / transversion mappings
    transitions = {'A':'G','G':'A','C':'T','T':'C'}
    transversions = {
        'A':['C','T'],
        'C':['A','G'],
        'G':['C','T'],
        'T':['A','G']
    }
    # Probability of choosing a transition vs. a transversion
    p_tr = titv_ratio / (titv_ratio + 2.0)
    
    # Collect LTR positions (0-based)
    ltr_positions = list(range(0, ltr_len)) + list(range(seq_len - ltr_len, seq_len))
    # Filter to only canonical bases
    eligible = [i for i in ltr_positions if seq[i].upper() in transitions]
    if not eligible or n_events == 0:
        return seq_str, 0.0
    
    mutated_sites = set()
    for _ in range(n_events):
        pos = random.choice(eligible)
        base = seq[pos].upper()
        if random.random() < p_tr:
            new = transitions[base]
        else:
            new = random.choice(transversions[base])
        # preserve case
        seq[pos] = new.lower() if seq[pos].islower() else new
        mutated_sites.add(pos)
    
    measured_pct = len(mutated_sites) / total_ltr_bases * 100.0
    return "".join(seq), measured_pct

def main():
    args = parse_args()
    if args.seed is not None:
        random.seed(args.seed)
    
    #pattern = re.compile(r'LTRlen(\d+)')
    pattern = re.compile(r'LTRlen:(\d+)')
    out_records = []
    stats = []

    for rec in SeqIO.parse(args.input, "fasta"):
        m = pattern.search(rec.description)
        if not m:
            print(f"Warning: no LTRlen found in header of {rec.id}, skipping.", file=sys.stderr)
            continue
        ltr_len = int(m.group(1))
        mutated_seq, meas_pct = mutate_ltr(str(rec.seq), ltr_len,
                                           args.mutation_percent,
                                           args.titv_ratio)
        rec.seq = Seq(mutated_seq)
        out_records.append(rec)
        stats.append((rec.id, meas_pct))

    # write output FASTA
    SeqIO.write(out_records, args.output, "fasta")

    # print measured mutation table
    print("SequenceID\tMeasured_LTR_Mutation(%)")
    for seqid, pct in stats:
        print(f"{seqid}\t{pct:.2f}")

if __name__ == "__main__":
    main()
