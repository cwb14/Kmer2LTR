#!/usr/bin/env python3

import argparse
import os
import shutil
import subprocess
from multiprocessing import Pool
from pathlib import Path

# Global args placeholder\ARGS = None

def run_cmd(cmd, verbose):
    """
    Run a command (list or string). If verbose, print the command.
    """
    if verbose:
        if isinstance(cmd, list):
            print("Running:", " ".join(cmd))
        else:
            print("Running:", cmd)
    subprocess.run(cmd, shell=not isinstance(cmd, list), check=True)


def run_cmd_output(cmd, output_path, stderr=None, verbose=False):
    """
    Run a command and redirect its stdout to output_path. Prints if verbose.
    """
    if verbose:
        if isinstance(cmd, list):
            print("Running:", " ".join(cmd), ">", output_path)
        else:
            print("Running:", cmd)
    with open(output_path, "w") as out:
        subprocess.run(cmd, shell=not isinstance(cmd, list), stdout=out, stderr=stderr, check=True)


def split_fasta(input_fasta, temp_dir):
    """
    Split a multi-sequence FASTA into one seq.fa per sequence under temp_dir/<seq_id>/
    """
    current_out = None
    with open(input_fasta) as infile:
        for line in infile:
            if line.startswith(">"):
                header = line[1:].strip().split()[0]
                dir_path = Path(temp_dir) / header
                dir_path.mkdir(parents=True, exist_ok=True)
                if current_out:
                    current_out.close()
                current_out = open(dir_path / "seq.fa", "w")
                current_out.write(line)
            else:
                if current_out:
                    current_out.write(line)
    if current_out:
        current_out.close()


def process_dir(dir_path):
    """
    Process one sequence directory: k-mer counting, filtering, LTR extraction, alignment, trimming, and final wfa.
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
            "python", "map_kmers_to_fasta.py",
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
        "python", "filter_kmers.py", str(combined_file),
        "--std-factor", str(std_factor),
        "-o", str(filtered_file)
    ], verbose)

    # Extract LTRs
    ltrs_fa = dir_path / "LTRs.fa"
    run_cmd([
        "python", "extract_ltrs.py",
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
        f"./wfa -E 10000 -e 10000 -o 10000 -O 10000 -u {args.mutation_rate} {dir_path / 'LTRs.aln.clean'}"
        " | cut -f1,5-"
    )

    if verbose:
        print("Running:", cmd_str)
    with open(out_file, "a") as ofile:
        subprocess.run(cmd_str, shell=True, check=True, stdout=ofile)


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

    # Split FASTA into separate seq.fa files
    split_fasta(args.input_fasta, args.temp_dir)

    # Set global arguments for process_dir
    global ARGS
    ARGS = args

    # Discover sequence directories
    dirs = [str(p) for p in Path(args.temp_dir).iterdir() if p.is_dir()]

    # Process in parallel
    with Pool(processes=args.threads) as pool:
        pool.map(process_dir, dirs)

    # Cleanup
    if args.keep_temp:
        print(f"Keeping temporary directory {args.temp_dir}")
    else:
        print(f"Removing temporary directory {args.temp_dir}")
        shutil.rmtree(args.temp_dir)

    print(f"All sequences processed. Output in {args.outfile}")
