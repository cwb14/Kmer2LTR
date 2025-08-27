#!/usr/bin/env python3
"""
LTR_mutator.py

Mutates only the LTR regions of intact LTR-RT sequences.
Now supports SNPs (with Ti/Tv control) and small indels limited to LTRs.

Usage:
    python LTR_mutator.py \
        -i LTR_RT_lib.fa \
        -o LTR_RT_lib_mutated.fa \
        -mp 15 \
        -TiTv 2.0 \
        -idel_size 1,10 \
        -idel_freq 0.1 \
        [--seed 42]

Notes:
- -mp is the overall mutation percentage over LTR bases (events with replacement).
- -idel_freq is the fraction of those mutation *events* that are indels.
  e.g., -mp 10 -idel_freq 0.2 → 8% SNP events & 2% indel events.
- Indels occur only inside LTRs; internal regions are untouched.
- Insertions and deletions are equally likely (50/50).
"""

import argparse
import re
import random
import sys
from Bio import SeqIO
from Bio.Seq import Seq

def parse_idel_size(s):
    try:
        parts = s.split(",")
        if len(parts) != 2:
            raise ValueError
        a, b = int(parts[0]), int(parts[1])
        if a <= 0 or b <= 0 or a > b:
            raise ValueError
        return (a, b)
    except Exception:
        raise argparse.ArgumentTypeError(
            "Invalid -idel_size. Use 'min,max' with positive ints and min<=max, e.g. '1,10'."
        )

def parse_args():
    parser = argparse.ArgumentParser(
        description="Mutate LTR regions (SNPs + optional small indels) of intact LTR-RT sequences"
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
        "-idel_size", "--idel_size", type=parse_idel_size, default=(1,1),
        help="Comma-separated min,max size for indels (e.g. '1,10')."
    )
    parser.add_argument(
        "-idel_freq", "--idel_freq", type=float, default=0.0,
        help="Fraction of mutation events that are indels (e.g. 0.1 for 10%%)."
    )
    parser.add_argument(
        "--seed", type=int, default=None,
        help="Random seed for reproducibility (optional)"
    )
    return parser.parse_args()

def mutate_ltr(seq_str, ltr_len, mp, titv_ratio, idel_size_range, idel_freq):
    """
    Apply SNP and indel mutation events to the two LTRs (first ltr_len bases and last ltr_len bases).
    Events are sampled WITH replacement (i.e., positions can be hit multiple times).
    Indels are restricted to occur within LTR boundaries; internal region is not modified.

    Returns (mutated_seq_str, measured_percent_unique_SNP_sites).
    - measured_percent_unique_SNP_sites is computed as unique SNP-mutated sites / (2*ltr_len) * 100
      (for backward compatibility with previous script behavior).
    """
    # Short-circuit if nothing to do
    if mp <= 0:
        return seq_str, 0.0

    seq = list(seq_str)
    original_len = len(seq)
    total_ltr_bases_initial = 2 * ltr_len

    # Number of total mutation events (SNP + indel)
    n_events = int(round(total_ltr_bases_initial * mp / 100.0))
    if n_events == 0:
        return seq_str, 0.0

    # Partition events into indels vs SNPs
    indel_events = int(round(n_events * max(0.0, min(1.0, idel_freq))))
    snp_events = n_events - indel_events

    # Define transition / transversion mappings
    transitions = {'A':'G','G':'A','C':'T','T':'C'}
    transversions = {'A':['C','T'], 'C':['A','G'], 'G':['C','T'], 'T':['A','G']}
    p_tr = titv_ratio / (titv_ratio + 2.0) if titv_ratio >= 0 else 0.0

    # Track dynamic LTR spans as sequence length changes with indels
    # Left LTR: [0, left_len), Right LTR: [len(seq)-right_len, len(seq))
    left_len = ltr_len
    right_len = ltr_len

    # For measuring SNP coverage like before
    snp_mutated_sites = set()

    def pick_ltr_side():
        """Randomly choose Left or Right LTR with probability proportional to current LTR lengths."""
        nonlocal left_len, right_len
        total = left_len + right_len
        if total <= 0:
            return None
        return 'L' if random.randrange(total) < left_len else 'R'

    def pick_pos_in_ltr(side):
        """Pick a random 0-based position within the chosen LTR in CURRENT sequence coordinates."""
        nonlocal left_len, right_len, seq
        if side == 'L':
            if left_len <= 0:
                return None
            return random.randrange(0, left_len)
        else:
            if right_len <= 0:
                return None
            # right LTR occupies the tail region
            start = len(seq) - right_len
            return random.randrange(start, len(seq))

    # SNP events
    for _ in range(snp_events):
        side = pick_ltr_side()
        if side is None:
            break
        attempts = 0
        while attempts < 20:
            pos = pick_pos_in_ltr(side)
            if pos is None:
                break
            base = seq[pos]
            up = base.upper()
            if up in transitions:
                if random.random() < p_tr:
                    new = transitions[up]
                else:
                    new = random.choice(transversions[up])
                # preserve case
                seq[pos] = new.lower() if base.islower() else new
                snp_mutated_sites.add(pos)
                break
            attempts += 1
        # if we fail to find a canonical base after several attempts, skip this event silently

    # Indel events
    min_size, max_size = idel_size_range
    for _ in range(indel_events):
        side = pick_ltr_side()
        if side is None:
            break

        # Choose insertion or deletion (50/50)
        is_insertion = (random.random() < 0.5)

        if is_insertion:
            # Choose insert position within chosen LTR
            pos = pick_pos_in_ltr(side)
            if pos is None:
                continue
            size = random.randint(min_size, max_size)

            # Insert a random DNA string of given size; uppercase by default
            ins = [random.choice("ACGT") for _ in range(size)]
            seq[pos:pos] = ins  # insertion

            # Adjust LTR lengths
            if side == 'L':
                left_len += size
            else:
                right_len += size

        else:
            # deletion
            # Pick a start position where at least 1 bp can be deleted within the LTR
            # Then cap deletion length so it doesn't cross out of the LTR.
            if side == 'L':
                if left_len <= 0:
                    continue
                start = random.randrange(0, left_len)
                # max deletable inside left LTR:
                max_del = min(max_size, left_len - start)
                if max_del <= 0:
                    continue
                size = random.randint(min_size, max_del)
                del seq[start:start+size]
                left_len -= size
            else:
                if right_len <= 0:
                    continue
                # compute current right LTR start after any prior edits
                rstart = len(seq) - right_len
                start = random.randrange(rstart, len(seq))
                # max deletable inside right LTR:
                max_del = min(max_size, len(seq) - start)
                # but also ensure we don't cross outside the right LTR end:
                max_del = min(max_del, (len(seq) - start))  # already ensured
                # and ensure we don't cross the left boundary of right LTR:
                max_del = min(max_del, len(seq) - start)  # safe guard
                # also clamp to stay within right_len
                max_del = min(max_del, (len(seq) - start))  # redundant but explicit
                # Explicit cap by right LTR boundary:
                max_del = min(max_del, (len(seq) - start))  # (kept simple; rstart recomputed post edits)
                # Real cap: number of bases to the end of sequence (which equals end of right LTR)
                max_del = min(max_size, len(seq) - start)
                # Also ensure we don't delete beyond the right LTR:
                max_del = min(max_del, (len(seq) - start))  # right LTR ends at seq end
                # Finally ensure we don't delete more than right_len total:
                max_del = min(max_del, right_len - (start - (len(seq) - right_len)))
                if max_del <= 0:
                    continue
                size = random.randint(min_size, max_del)
                del seq[start:start+size]
                right_len -= size

            # After deletion in either LTR, the opposite LTR stays the same, but
            # note that for right LTR its start index shifts automatically with seq length.

    measured_pct_snp_sites = (len(snp_mutated_sites) / total_ltr_bases_initial) * 100.0
    return "".join(seq), measured_pct_snp_sites

def fmt_pct(x):
    # Nice compact formatting for percentages
    if abs(x - round(x)) < 1e-9:
        return f"{int(round(x))}%"
    return f"{x:.2f}%"

def main():
    args = parse_args()
    if args.seed is not None:
        random.seed(args.seed)

    # Sanity clamp for idel_freq
    if args.idel_freq < 0 or args.idel_freq > 1:
        print("Error: -idel_freq must be between 0 and 1 (inclusive).", file=sys.stderr)
        sys.exit(1)

    # Header pattern: expects '... LTRlen:<int> ...'
    pattern = re.compile(r'LTRlen:(\d+)')

    # Print the conceptual breakdown so it's obvious
    indel_pct = args.mutation_percent * args.idel_freq
    snp_pct = args.mutation_percent - indel_pct
    min_sz, max_sz = args.idel_size
    print(f"{fmt_pct(args.mutation_percent)} mutations: {fmt_pct(snp_pct)} SNPs & "
          f"{fmt_pct(indel_pct)} indel ({min_sz}–{max_sz} bp in length)")

    out_records = []
    stats = []

    for rec in SeqIO.parse(args.input, "fasta"):
        m = pattern.search(rec.description)
        if not m:
            print(f"Warning: no LTRlen found in header of {rec.id}, skipping.", file=sys.stderr)
            continue
        ltr_len = int(m.group(1))
        mutated_seq, meas_snp_pct = mutate_ltr(
            str(rec.seq),
            ltr_len,
            args.mutation_percent,
            args.titv_ratio,
            args.idel_size,
            args.idel_freq
        )
        rec.seq = Seq(mutated_seq)
        out_records.append(rec)
        stats.append((rec.id, meas_snp_pct))

    # write output FASTA
    SeqIO.write(out_records, args.output, "fasta")

    # print measured mutation table (SNP unique-site coverage, as in the original script)
    print("SequenceID\tMeasured_LTR_SNP_UniqueSites(%)")
    for seqid, pct in stats:
        print(f"{seqid}\t{pct:.2f}")

if __name__ == "__main__":
    main()
