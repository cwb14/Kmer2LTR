#!/usr/bin/env python3
import sys
import argparse
import statistics

def filter_monotonic_and_outliers(fname_in, fname_out=None, std_factor=2):
    """
    Reads a TSV with at least 5 columns, filters rows so that:
      1. 'start2' (column 5, zero-based index 4) is non-decreasing when sorted by 'start1' (column 3, index 2).
      2. Any row whose start1 or start2 value lies more than `std_factor` standard deviations
         away from the mean of its respective distribution (after monotonic filtering) is removed.
    Writes the header plus the remaining rows to stdout or to fname_out if provided.
    """
    # Read all lines
    with open(fname_in, 'r') as f:
        lines = f.readlines()
    if not lines:
        return

    # Extract header
    header = lines[0].rstrip("\n")
    data = []
    for line in lines[1:]:
        parts = line.rstrip("\n").split('\t')
        try:
            s1 = int(parts[2])
            s2 = int(parts[4])
            data.append((s1, s2, parts))
        except (IndexError, ValueError):
            # skip bad lines
            continue

    # 1) Sort by start1 and enforce non-decreasing start2
    data.sort(key=lambda x: x[0])
    monotonic = []
    prev_s2 = None
    for s1, s2, parts in data:
        if prev_s2 is None or s2 >= prev_s2:
            monotonic.append((s1, s2, parts))
            prev_s2 = s2

    if not monotonic:
        # nothing to output beyond header
        out = sys.stdout if fname_out is None else open(fname_out, 'w')
        print(header, file=out)
        if fname_out:
            out.close()
        return

    # 2) Compute mean and std for start1 and start2
    s1_vals = [row[0] for row in monotonic]
    s2_vals = [row[1] for row in monotonic]
    mean1, mean2 = statistics.mean(s1_vals), statistics.mean(s2_vals)
    std1 = statistics.pstdev(s1_vals)
    std2 = statistics.pstdev(s2_vals)

    # If std is zero, skip outlier filtering on that axis
    if std1 > 0:
        low1, high1 = mean1 - std_factor * std1, mean1 + std_factor * std1
    else:
        low1, high1 = float('-inf'), float('inf')
    if std2 > 0:
        low2, high2 = mean2 - std_factor * std2, mean2 + std_factor * std2
    else:
        low2, high2 = float('-inf'), float('inf')

    # 3) Remove rows where either start1 or start2 is an outlier
    final = [
        parts
        for (s1, s2, parts) in monotonic
        if low1 <= s1 <= high1 and low2 <= s2 <= high2
    ]

    # Write output
    out = sys.stdout if fname_out is None else open(fname_out, 'w')
    print(header, file=out)
    for parts in final:
        print("\t".join(parts), file=out)
    if fname_out:
        out.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Keep only rows where start2 is non-decreasing as start1 increases, "
                    "then remove outliers in start1 or start2 based on a Z-score threshold."
    )
    parser.add_argument("input", help="Input TSV file")
    parser.add_argument("-o", "--output", help="Output TSV (defaults to stdout)")
    parser.add_argument(
        "--std-factor", type=float, default=2.0,
        help="Number of standard deviations for outlier cutoff (default: 2)"
    )
    args = parser.parse_args()

    filter_monotonic_and_outliers(args.input, args.output, args.std_factor)
