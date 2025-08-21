#!/usr/bin/env python3

import argparse
import os
import shutil
import subprocess
import math
from multiprocessing import Pool
from pathlib import Path
import re

SCRIPT_DIR = Path(__file__).parent.resolve()

# Global args/domains placeholders
ARGS = None
DOMAINS = None  # dict[str,int] or None


def sanitize_name(name: str) -> str:
    """
    Make a filesystem-safe folder name from a FASTA header token.
    Replaces path separators with a Unicode division slash '∕' (U+2215).
    Also strips trailing dots/spaces which can be problematic on some FS.
    """
    safe = name.replace("/", "∕").replace("\\", "∕")
    safe = safe.strip()
    if safe in {".", "..", ""}:
        safe = f"header_{abs(hash(name))}"
    return safe

def sum_ungapped_bp(fasta_path: Path) -> int:
    """
    Count total basepairs across all sequences in a FASTA, excluding gaps ('-').
    Headers ('>') and whitespace are ignored.
    """
    total = 0
    with open(fasta_path) as fh:
        for line in fh:
            if not line or line.startswith(">"):
                continue
            s = line.strip()
            # count all non-gap characters
            total += sum(1 for c in s if c != "-")
    return total

def has_data_rows(path: Path) -> bool:
    """
    True if file exists and contains at least one non-empty, non-header line.
    Treats a header-only file as empty.
    """
    try:
        with open(path) as fh:
            for line in fh:
                s = line.strip()
                if not s:
                    continue
                if s.lower().startswith("kmer"):
                    continue
                return True
    except FileNotFoundError:
        return False
    return False


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
            print("Running:", cmd, ">", output_path)
    with open(output_path, "w") as out:
        subprocess.run(cmd, shell=not isinstance(cmd, list), stdout=out, stderr=stderr, check=True)

def read_domains(domains_tsv: str):
    """
    Read domains TSV: <name> <TAB> <ltr_len>. Ignores blank/comment lines.
    Returns dict[name] = int(ltr_len)
    """
    if not domains_tsv:
        return None
    mapping = {}
    with open(domains_tsv) as fh:
        for ln, line in enumerate(fh, 1):
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split("\t")
            if len(parts) < 2:
                print(f"[DOMAINS] Skipping malformed line {ln}: {s}")
                continue
            name, ltr = parts[0], parts[1]
            try:
                ltr_len = int(ltr)
            except ValueError:
                print(f"[DOMAINS] Skipping non-integer LTR length at line {ln}: {s}")
                continue
            if ltr_len <= 0:
                print(f"[DOMAINS] Skipping non-positive LTR length at line {ln}: {s}")
                continue
            mapping[name] = ltr_len
    if not mapping:
        print("[DOMAINS] No valid entries parsed; ignoring domains file.")
        return None
    print(f"[DOMAINS] Loaded {len(mapping)} entries from {domains_tsv}")
    return mapping


def split_fasta(input_fasta, temp_dir):
    """
    Split a multi-sequence FASTA into one seq.fa per sequence under
    temp_dir/<sanitized_seq_id>[/__dupN]/, preserving the original header in seq.fa.

    Writes a header_map.tsv mapping: sanitized_name <TAB> original_header_token

    If duplicate header tokens occur, create unique dirs with a __dupN suffix and
    write a DUPLICATE_NAME marker file in *all* such directories so we can safely
    fall back to the full pipeline for those.
    """
    temp_dir = Path(temp_dir)
    map_path = temp_dir / "header_map.tsv"
    counts = {}
    with open(input_fasta) as infile, open(map_path, "w") as map_out:
        current_out = None
        current_dir = None
        for line in infile:
            if line.startswith(">"):
                header_line = line[1:].strip()
                header_token = header_line.split()[0]  # token before first whitespace
                counts[header_token] = counts.get(header_token, 0) + 1

                base_safe = sanitize_name(header_token)
                # First occurrence: base name; duplicates get a suffix
                if counts[header_token] == 1:
                    safe = base_safe
                else:
                    safe = f"{base_safe}__dup{counts[header_token]}"

                dir_path = temp_dir / safe
                dir_path.mkdir(parents=True, exist_ok=True)

                # If we just learned this token is a duplicate (count==2),
                # go back and mark the first directory as duplicate as well.
                if counts[header_token] == 2:
                    first_dir = temp_dir / base_safe
                    try:
                        (first_dir / "DUPLICATE_NAME").write_text("")
                    except Exception:
                        pass

                # Always mark duplicates (count>1)
                if counts[header_token] > 1:
                    try:
                        (dir_path / "DUPLICATE_NAME").write_text("")
                    except Exception:
                        pass

                # mapping: sanitized dir -> original header token
                map_out.write(f"{safe}\t{header_token}\n")

                if current_out:
                    current_out.close()
                current_dir = dir_path
                current_out = open(dir_path / "seq.fa", "w")
                current_out.write(">" + header_line + "\n")
            else:
                if current_out:
                    current_out.write(line)
        if current_out:
            current_out.close()


def read_seq_string(seq_fa: Path) -> str:
    """Read a single FASTA sequence (concatenate all non-header lines)."""
    seq_chunks = []
    with open(seq_fa) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            seq_chunks.append(line.strip())
    return "".join(seq_chunks)


def try_fast_path(dir_path: Path, header_token: str, seq_len: int) -> bool:
    """
    Fast-path when DOMAINS TSV provides LTR length:
      - Still runs MAFFT and trimAl (we only skip k-mer discovery/extraction).
      - Extract 5′ and 3′ segments of length min(ltr_len + extension, floor(seq_len/2)).
      - Apply the same trimmed-vs-raw safety short-circuit before WFA.
    Returns True if fast-path executed (even if it ends up skipping due to checks),
    False if we should run the full pipeline.
    """
    global ARGS, DOMAINS
    if not DOMAINS:
        return False

    # Duplicate name? -> uncertain mapping -> full pipeline
    if (dir_path / "DUPLICATE_NAME").exists():
        if ARGS.verbose:
            print(f"[FAST-PATH SKIP] Duplicate header name detected for '{header_token}'.")
        return False

    if header_token not in DOMAINS:
        return False

    ltr_len = DOMAINS[header_token]

    # Safety check: length_of_element > 2*ltr_len + 50
    if not (seq_len > (2 * ltr_len + 50)):
        if ARGS.verbose:
            print(f"[FAST-PATH SKIP] Safety check failed for '{header_token}': "
                  f"seq_len={seq_len} <= 2*{ltr_len}+50")
        return False

    # Extract first/last (ltr_len + extension), capped at seq_len/2
    seq = read_seq_string(dir_path / "seq.fa")
    ext_len = min(ltr_len + ARGS.extension, seq_len // 2)
    five_p = seq[:ext_len]
    three_p = seq[-ext_len:]

    # Write LTRs.fa to feed into MAFFT -> trimAl
    ltrs_fa = dir_path / "LTRs.fa"
    with open(ltrs_fa, "w") as out:
        out.write(f">{header_token}_5prime_LTR_len{ltr_len}_ext{ext_len}\n")
        out.write(five_p + "\n")
        out.write(f">{header_token}_3prime_LTR_len{ltr_len}_ext{ext_len}\n")
        out.write(three_p + "\n")

    if ARGS.verbose:
        print(f"[FAST-PATH] {dir_path.name}: extracted 5′/3′ LTRs (len={ltr_len}, used={ext_len}) "
              f"from domains TSV; running MAFFT + trimAl.")

    # Align and trim (same as heavy pipeline)
    aln_fa = dir_path / "LTRs.aln.fa"
    run_cmd_output([
        "mafft", "--auto", "--thread", "1", str(ltrs_fa)
    ], aln_fa, stderr=subprocess.DEVNULL, verbose=ARGS.verbose)

    run_cmd([
        "trimal",
        "-automated1",
        "-keepheader",
        "-in", str(aln_fa),
        "-out", str(dir_path / "LTRs.aln.clean")
    ], ARGS.verbose)

    # Apply the same short-circuit as the heavy pipeline
    try:
        raw_bp = sum_ungapped_bp(aln_fa)
        clean_bp = sum_ungapped_bp(dir_path / "LTRs.aln.clean")
    except FileNotFoundError as e:
        print(f"[SKIP] Missing alignment file for {dir_path.name}: {e}")
        try:
            (dir_path / "SKIPPED.missing_alignment").write_text(str(e))
        except Exception:
            pass
        return True

    threshold = 0.6 * raw_bp
    if (clean_bp + ARGS.extension) <= threshold:
        print(
            f"[SKIP] Trimmed alignment too short for {dir_path.name}: "
            f"clean_bp({clean_bp}) + ext({ARGS.extension}) <= 0.6 * raw_bp({raw_bp:.0f})"
        )
        try:
            (dir_path / "SKIPPED.too_trimmed").write_text(
                f"raw_bp={raw_bp}\nclean_bp={clean_bp}\nextension={ARGS.extension}\nthreshold={threshold}\n"
            )
        except Exception:
            pass
        return True

    # Final WFA and append results (same as heavy pipeline)
    cmd_str = (
        f"{SCRIPT_DIR / 'wfa'} -E 10000 -e 10000 -o 10000 -O 10000 -u {ARGS.mutation_rate} {dir_path / 'LTRs.aln.clean'}"
        " | cut -f1,5- | sed -e 's/-0\\.000000/0.000000/g' | "
        "awk -F'\t' '{if($1~/_(5prime|3prime)/) sub(/_(5prime|3prime).*$/, \"\", $1); print}' OFS='\t' | "
        "awk '$6 <= 0.35' "
    )

    if ARGS.verbose:
        print("Running:", cmd_str)
    with open(ARGS.outfile, "a") as ofile:
        subprocess.run(cmd_str, shell=True, check=True, stdout=ofile)

    return True


def process_dir(dir_path):
    """
    Process one sequence directory:
      - If domains TSV provides a safe fast-path, extract LTRs directly (with extension cap)
        and run MAFFT -> trimAl -> WFA.
      - Otherwise run the original heavy pipeline (k-mers -> map -> filter -> extract LTRs ->
        MAFFT -> trimAl -> WFA).
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
    # Header token (the same notion used for domains lookup)
    with open(seq_fa) as fh:
        first = next((ln for ln in fh if ln.startswith(">")), None)
    header_token = first[1:].strip().split()[0] if first else dir_path.name

    seq_len = sum(len(line.strip()) for line in seq_fa.open() if not line.startswith(">"))

    # ===== FAST-PATH TRY =====
    try:
        if try_fast_path(dir_path, header_token, seq_len):
            return
    except Exception as e:
        print(f"[FAST-PATH ERROR] {header_token}: {e}. Falling back to full pipeline.")

    # ===== ORIGINAL HEAVY PIPELINE =====

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

    filtered_file = dir_path / f"kmer_duet_k{kmin}_{kmax}.mapped.filtered"
    run_cmd([
        "python", str(SCRIPT_DIR / "filter_kmers.py"), str(combined_file),
        "--std-factor", str(std_factor),
        "-o", str(filtered_file)
    ], verbose)

    # Short-circuit if empty or header-only
    if not has_data_rows(filtered_file):
        print(f"[SKIP] No k-mer pairs after filtering for {dir_path.name} "
              f"({filtered_file} is empty/header-only).")
        try:
            (dir_path / "SKIPPED.no_filtered_pairs").write_text("")
        except Exception:
            pass
        return

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
        "trimal",
        "-automated1",
        "-keepheader",  # For benchmarking LTR boundry classification, its helpful to have the LTR length in the header. Without "-keepheader", trimal modifies 'Gypsy-1_CiLe_Cimex#LTR/Ty3~LTRlen:490' to 'Gypsy-1_CiLe_Cimex#LTR/Ty3~LTRlen'.
        "-in", str(aln_fa),
        "-out", str(dir_path / "LTRs.aln.clean")
    ], verbose)

    # ----- Short-circuit if trimmed alignment is too short vs raw alignment -----
    # MAFFT and trimal dont filter spurious likely false positives.
    # This is aimed at conservatively cleaning some false positives.
    # 0.6 is used, but may need adjusted.
    try:
        raw_bp = sum_ungapped_bp(aln_fa)                 # from "LTRs.aln.fa"
        clean_bp = sum_ungapped_bp(dir_path / "LTRs.aln.clean")  # from "LTRs.aln.clean"
    except FileNotFoundError as e:
        print(f"[SKIP] Missing alignment file for {dir_path.name}: {e}")
        try:
            (dir_path / "SKIPPED.missing_alignment").write_text(str(e))
        except Exception:
            pass
        return

    # Require: sum_bp(clean) + extension > 0.6 * sum_bp(raw)
    threshold = 0.6 * raw_bp
    if (clean_bp + extension) <= threshold:
        print(
            f"[SKIP] Trimmed alignment too short for {dir_path.name}: "
            f"clean_bp({clean_bp}) + ext({extension}) <= 0.6 * raw_bp({raw_bp:.0f})"
        )
        try:
            (dir_path / "SKIPPED.too_trimmed").write_text(
                f"raw_bp={raw_bp}\nclean_bp={clean_bp}\nextension={extension}\nthreshold={threshold}\n"
            )
        except Exception:
            pass
        return

    # Final WFA and append results
    cmd_str = (
        f"{SCRIPT_DIR / 'wfa'} -E 10000 -e 10000 -o 10000 -O 10000 -u {args.mutation_rate} {dir_path / 'LTRs.aln.clean'}"
        " | cut -f1,5- | sed -e 's/-0\\.000000/0.000000/g' | "
        "awk -F'\t' '{if($1~/_(5prime|3prime)/) sub(/_(5prime|3prime).*$/, \"\", $1); print}' OFS='\t' | "
        "awk '$6 <= 0.35' " # Filter out likely false positives.
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
                continue
            try:
                L = int(parts[1])
                ts = float(parts[3])
                tv = float(parts[4])
            except ValueError:
                continue
            total_len += L
            total_ts += ts
            total_tv += tv

    summary_path = main_out_path + ".summary"
    if total_len == 0:
        with open(summary_path, "w") as out:
            out.write("total_length\t0\n")
            out.write("total_transitions\t0\n")
            out.write("total_transversions\t0\n")
            out.write("raw_d\tNA\nJC69_d\tNA\nK2P_d\tNA\n")
        return

    P = total_ts / total_len
    Q = total_tv / total_len
    raw_d = P + Q

    jc_term = 1.0 - (4.0/3.0) * raw_d
    JC69_d = None if jc_term <= 0 else -0.75 * math.log(jc_term)

    k2p_t1 = 1.0 - 2.0 * P - Q
    k2p_t2 = 1.0 - 2.0 * Q
    K2P_d = None if (k2p_t1 <= 0 or k2p_t2 <= 0) else (-0.5 * math.log(k2p_t1) - 0.25 * math.log(k2p_t2))

    with open(summary_path, "w") as out:
        out.write(f"total_length\t{total_len}\n")
        out.write(f"total_transitions\t{total_ts:.0f}\n")
        out.write(f"total_transversions\t{total_tv:.0f}\n")
        out.write(f"raw_d\t{raw_d:.6f}\n")
        out.write(f"JC69_d\t{('NA' if JC69_d is None else f'{JC69_d:.6f}')}\n")
        out.write(f"K2P_d\t{('NA' if K2P_d is None else f'{K2P_d:.6f}')}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process multi-seq LTR-RT FASTA to extract and align LTRs. "
                    "Optionally fast-path known LTR lengths via a domains TSV."
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
        "-D", "--domains", dest="domains_tsv", default=None,
        help=("Optional TSV (name\\tLTR_len). If provided, sequences whose names appear here "
              "and pass safety checks will SKIP k-mer/duet discovery and directly extract the "
              "first/last min(LTR_len + extension, seq_len/2) bp, then run MAFFT + trimAl before WFA.")
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

    # Load domains mapping (if provided)
    DOMAINS = read_domains(args.domains_tsv)

    # Split FASTA into separate seq.fa files (safe/unique dir names; original headers preserved)
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
