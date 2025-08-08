#!/usr/bin/env python3

import argparse
import os
import shutil
import subprocess
import math
from multiprocessing import Pool
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent.resolve()

# Global args placeholder
ARGS = None


def sanitize_name(name: str) -> str:
    """
    Make a filesystem-safe folder name from a FASTA header token.
    Replaces path separators with a Unicode division slash '∕' (U+2215),
    which looks like '/' but is not a path separator.
    Also strips trailing dots/spaces which can be problematic on some FS.
    """
    # Replace forward/back slashes
    safe = name.replace("/", "∕").replace("\\", "∕")
    # Optionally, forbid leading/trailing whitespace
    safe = safe.strip()
    # Avoid names like '.' or '..'
    if safe in {".", "..", ""}:
        safe = f"header_{abs(hash(name))}"
    return safe


def run_cmd(cmd, verbose):
    """
    Run a command (list or string). If verbose, print the command.
    """
    if verbose:
        if isinstance(cmd, list):
            print("Running:", " ".join(map(str, cmd)))
        else:
            print("Running:", cmd)
    subprocess.run(cmd, shell=not isinstance(cmd, list), check=True)


def run_cmd_output(cmd, output_path, stderr=None, verbose=False):
    """
    Run a command and redirect its stdout to output_path. Prints if verbose.
    """
    if verbose:
        if isinstance(cmd, list):
            print("Running:", " ".join(map(str, cmd)), ">", output_path)
        else:
            print("Running:", cmd)
    with open(output_path, "w") as out:
        subprocess.run(cmd, shell=not isinstance(cmd, list), stdout=out, stderr=stderr, check=True)


def split_fasta(input_fasta, temp_dir):
    """
    Split a multi-sequence FASTA into one seq.fa per sequence under
    temp_dir/<sanitized_seq_id>/, preserving the original header in seq.fa.

    Writes a header_map.tsv mapping:  sanitized_name <TAB> original_header
    """
    temp_dir = Path(temp_dir)
    map_path = temp_dir / "header_map.tsv"
    with open(input_fasta) as infile, open(map_path, "w") as map_out:
        current_out = None
        current_dir = None
        for line in infile:
            if line.startswith(">"):
                header = line[1:].strip().split()[0]  # token before first whitespace
                safe = sanitize_name(header)
                dir_path = temp_dir / safe
                dir_path.mkdir(parents=True, exist_ok=True)
                # write mapping once per header encountered
                map_out.write(f"{safe}\t{header}\n")
                if current_out:
                    current_out.close()
                current_dir = dir_path
                current_out = open(dir_path / "seq.fa", "w")
                # Write original header as-is (so no loss of info for downstream)
                current_out.write(line)
            else:
                if current_out:
                    current_out.write(line)
        if current_out:
            current_out.close()


def process_dir(dir_path):
    """
    Process one sequence directory: k-mer counting, filtering, LTR extraction,
    alignment, trimming, and final wfa.
    """
    args = ARGS
    verbose = args.verbose
    out_file = args.outfile
    dist = args.dist
    kmin = args.kmin
    kmax = args.kmax
    std_factor = args.std_factor
    extension = args.extension

    dir_path = Path(dir_path)
    print(f"Processing {dir_path}")
    seq_fa = dir_path / "seq.fa"
    seq_len = sum(len(line.strip()) for line in seq_fa.open() if not line.startswith(">"))

    # K-mer counting and mapping
    for k in range(kmin, kmax + 1):
        jf_file = dir_path / "seq.jf"
        run_cmd([
            "jellyfish", "count", "-m", str(k), "-s", str(seq_len), "-t", "1",
            "-o", str(jf_file), str(seq_fa)
        ], verbose)

        run_cmd([
            "jellyfish", "dump", "-c", "-L", "2", "-U", "2",
            "-o", str(dir_path / "kmer_duet.txt"), str(jf_file)
        ], verbose)

        mapped_file = dir_path / f"kmer_duet_k{k}.mapped"
        run_cmd_output([
            "python", str(SCRIPT_DIR / "map_kmers_to_fasta.py"),
            str(dir_path / "kmer_duet.txt"), str(seq_fa),
            "-d", str(dist)
        ], mapped_file, verbose=verbose)

    # Concatenate mapped kmers
    combined_file = dir_path / f"kmer_duet_k{kmin}_{kmax}.mapped"
    with open(combined_file, "w") as outfile:
        for k in range(kmin, kmax + 1):
            part = dir_path / f"kmer_duet_k{k}.mapped"
            with open(part) as infile:
                outfile.write(infile.read())

    # Filter kmers
    filtered_file = dir_path / f"kmer_duet_k{kmin}_{kmax}.mapped.filtered"
    run_cmd([
        "python", str(SCRIPT_DIR / "filter_kmers.py"), str(combined_file),
        "--std-factor", str(std_factor),
        "-o", str(filtered_file)
    ], verbose)

    # Extract LTRs
    ltrs_fa = dir_path / "LTRs.fa"
    run_cmd([
        "python", str(SCRIPT_DIR / "extract_ltrs.py"),
        "-e", str(extension),
        str(filtered_file), str(seq_fa), str(ltrs_fa)
    ], verbose)

    # Align and trim
    aln_fa = dir_path / "LTRs.aln.fa"
    run_cmd_output([
        "mafft", "--auto", "--thread", "1", str(ltrs_fa)
    ], aln_fa, stderr=subprocess.DEVNULL, verbose=verbose)

    run_cmd([
        "trimal", "-automated1", "-in", str(aln_fa), "-out",
        str(dir_path / "LTRs.aln.clean")
    ], verbose)

    # Final alignment and append results
    cmd_str = (
        f"{SCRIPT_DIR / 'wfa'} -E 10000 -e 10000 -o 10000 -O 10000 -u {args.mutation_rate} {dir_path / 'LTRs.aln.clean'}"
        "| cut -f1,5- | sed -e 's/-0\\.000000/0.000000/g'"
    )

    if verbose:
        print("Running:", cmd_str)
    with open(out_file, "a") as ofile:
        subprocess.run(cmd_str, shell=True, check=True, stdout=ofile)

def write_summary(main_out_path: str):
    """
    Read the main outfile (one row per element), sum length/transition/transversion,
    and write pooled distances to '<outfile>.summary'.

    Columns (no header):
    1 LTR-RT  2 LTR_LEN  3 substitutions  4 transitions  5 transversions
    6 p-dist  7 p-time   8 JC69-dist      9 JC69-time     10 K2P-dist   11 K2P-time
    """
    total_len = 0
    total_ts = 0.0
    total_tv = 0.0

    with open(main_out_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 11:
                # Skip malformed lines
                continue
            # Parse needed fields
            try:
                L = int(parts[1])
                ts = float(parts[3])
                tv = float(parts[4])
            except ValueError:
                continue
            total_len += L
            total_ts += ts
            total_tv += tv

    # Guard against empty/zero totals
    if total_len == 0:
        summary_path = main_out_path + ".summary"
        with open(summary_path, "w") as out:
            out.write("total_length\t0\n")
            out.write("total_transitions\t0\n")
            out.write("total_transversions\t0\n")
            out.write("raw_d\tNA\nJC69_d\tNA\nK2P_d\tNA\n")
        return

    # Pooled proportions
    P = total_ts / total_len  # transitions / site
    Q = total_tv / total_len  # transversions / site
    raw_d = P + Q

    # JC69: d = -3/4 * ln(1 - 4/3 * p)
    jc_term = 1.0 - (4.0/3.0) * raw_d
    JC69_d = None if jc_term <= 0 else -0.75 * math.log(jc_term)

    # K2P: d = -1/2 ln(1 - 2P - Q) - 1/4 ln(1 - 2Q)
    k2p_t1 = 1.0 - 2.0 * P - Q
    k2p_t2 = 1.0 - 2.0 * Q
    K2P_d = None if (k2p_t1 <= 0 or k2p_t2 <= 0) else (-0.5 * math.log(k2p_t1) - 0.25 * math.log(k2p_t2))

    # Write summary
    summary_path = main_out_path + ".summary"
    with open(summary_path, "w") as out:
        out.write(f"total_length\t{total_len}\n")
        out.write(f"total_transitions\t{total_ts:.0f}\n")
        out.write(f"total_transversions\t{total_tv:.0f}\n")
        out.write(f"raw_d\t{raw_d:.6f}\n")
        out.write(f"JC69_d\t{('NA' if JC69_d is None else f'{JC69_d:.6f}')}\n")
        out.write(f"K2P_d\t{('NA' if K2P_d is None else f'{K2P_d:.6f}')}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process multi-seq LTR-RT FASTA to extract and align LTRs."
    )
    parser.add_argument(
        "-k", action="store_true", dest="keep_temp",
        help="Keep temp directory after processing."
    )
    parser.add_argument(
        "-v", action="store_true", dest="verbose",
        help="Verbose mode; print each command before executing."
    )
    parser.add_argument(
        "-d", type=int, default=80, dest="dist",
        help="Minimum distance between kmer pairs (default: 80)."
    )
    parser.add_argument(
        "-l", type=int, default=8, dest="kmin",
        help="Minimum kmer length (default: 8)."
    )
    parser.add_argument(
        "-U", type=int, default=12, dest="kmax",
        help="Maximum kmer length (default: 12)."
    )
    parser.add_argument(
        "-u", type=float, default=3e-8, dest="mutation_rate",
        help="Mutation rate μ (default: 3e-8)."
    )
    parser.add_argument(
        "-f", type=float, default=2.0, dest="std_factor",
        help="Standard deviation factor for kmer filtering."
    )
    parser.add_argument(
        "-e", type=int, default=65, dest="extension",
        help="Extension length for LTR extraction (default: 65)."
    )
    parser.add_argument(
        "-t", default="./temp", dest="temp_dir",
        help="Temporary directory name (default: ./temp)."
    )
    parser.add_argument(
        "-o", default="./LTRs.alns.results", dest="outfile",
        help="Output filename (default: ./LTRs.alns.results)."
    )
    parser.add_argument(
        "-p", type=int, default=20, dest="threads",
        help="Number of parallel threads (default: 20)."
    )
    parser.add_argument(
        "input_fasta",
        help="Path to multi-sequence LTR-RT FASTA file."
    )

    args = parser.parse_args()

    # Prepare workspace
    if os.path.exists(args.temp_dir):
        shutil.rmtree(args.temp_dir)
    os.makedirs(args.temp_dir, exist_ok=True)
    open(args.outfile, "w").close()

    # Split FASTA into separate seq.fa files (safe dir names; original headers preserved)
    split_fasta(args.input_fasta, args.temp_dir)

    # Set global arguments for process_dir
    ARGS = args

    # Discover sequence directories
    dirs = [str(p) for p in Path(args.temp_dir).iterdir() if p.is_dir()]

    # Process in parallel
    with Pool(processes=args.threads) as pool:
        pool.map(process_dir, dirs)

    write_summary(args.outfile)

    # Cleanup
    if args.keep_temp:
        print(f"Keeping temporary directory {args.temp_dir}")
    else:
        print(f"Removing temporary directory {args.temp_dir}")
        shutil.rmtree(args.temp_dir)

    print(f"All sequences processed. Output in {args.outfile}")
