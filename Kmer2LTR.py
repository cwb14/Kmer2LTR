#!/usr/bin/env python3

import argparse
import os
import shutil
import subprocess
import math
from multiprocessing import Pool
from pathlib import Path
import re
import shlex
import sys
import time

q = shlex.quote

SCRIPT_DIR = Path(__file__).parent.resolve()

# Global args/domains placeholders
ARGS = None
DOMAINS = None  # dict[str,int] or None
LOGFILE = None

def sanitize_name(name: str) -> str:
    """
    Make a filesystem- and shell-safe folder name from a FASTA header token.
    - Replace path separators with Unicode division slash '∕' (U+2215).
    - Replace shell metacharacters with visually similar fullwidth variants or underscores.
    - Strip trailing dots/spaces which can be problematic on some FS.
    """
    # Friendly translations for the worst offenders
    trans = str.maketrans({
        "/": "∕",
        "\\": "∕",
        ":": "꞉",   # U+A789 modifier letter colon (or use "：" U+FF1A)
        ";": "；",   # fullwidth semicolon
        "&": "＆",
        "|": "︱",
        " ": "_",
    })
    safe = name.translate(trans)

    # Replace any remaining unsafe characters with underscore
    safe = re.sub(r"[<>?*'\"`$(){}\[\]#\t\n\r]", "_", safe)

    # Trim and guard reserved names
    safe = safe.strip().rstrip(".")
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

def _format_hms(seconds: float) -> str:
    """Return 'Hh:Mm:Ss' for a nonnegative seconds value."""
    if not math.isfinite(seconds) or seconds < 0:
        return "0h:0m:0s"
    s = int(round(seconds))
    h, s = divmod(s, 3600)
    m, s = divmod(s, 60)
    return f"{h}h:{m}m:{s}s"

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

def log_msg(msg: str):
    """Write warnings/skips/etc. to the logfile instead of cluttering stderr."""
    global LOGFILE
    if LOGFILE:
        LOGFILE.write(msg.strip() + "\n")
        LOGFILE.flush()

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


def parse_proposed_ltr_len_from_aln(aln_fa: Path) -> int | None:
    """
    Look for a 5' header line like:
    >... 5p:1-2688
    Return the integer (e.g., 2688) or None if not found.
    """
    try:
        with open(aln_fa) as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                # tolerate tabs/spaces; search only the 5p tag
                m = re.search(r"5p:1-(\d+)", line)
                if m:
                    return int(m.group(1))
    except FileNotFoundError:
        return None
    return None


def try_fast_path(dir_path: Path, header_token: str, seq_len: int) -> bool:
    """
    Fast-path when domains TSV provides LTR length:
      - Safety checks 
      - Extract first/last min(ltr_len + extension, seq_len//2) bp
      - Write LTRs.fa
      - Run MAFFT -> trimal -> WFA, then append to main outfile

    Returns True if fast-path executed (regardless of downstream success/exception),
    False if we should run the full k-mer pipeline.
    """
    global ARGS, DOMAINS
    if not DOMAINS:
        return False

    # Duplicate name? -> uncertain mapping -> full pipeline unless user overrides.
    if (dir_path / "DUPLICATE_NAME").exists() and not ARGS.assume_dup_same_ltr:
        if ARGS.verbose:
            print(f"[FAST-PATH SKIP] Duplicate header name detected for '{header_token}'. "
                  f"Use --assume-duplicate-same-ltr to override (risky).")
        return False
    # If ARGS.assume_dup_same_ltr is True, proceed assuming all duplicates with this
    # header token share the same LTR length from DOMAINS.

    ltr_len = DOMAINS[header_token]

    # Safety check: length_of_element > 2*ltr_len + 50 (unchanged)
    if not (seq_len > (2 * ltr_len + 50)):
        if ARGS.verbose:
            print(f"[FAST-PATH SKIP] Safety check failed for '{header_token}': "
                  f"seq_len={seq_len} <= 2*{ltr_len}+50")
        return False

    # Determine extraction length: min(ltr_len + extension, seq_len//2)
    ext_len = ltr_len + ARGS.extension
    cap = seq_len // 2
    if ext_len > cap:
        ext_len = cap

    # Perform direct extraction into LTRs.fa (not aln.clean)
    seq = read_seq_string(dir_path / "seq.fa")
    five_p = seq[:ext_len]
    three_p = seq[-ext_len:]

    ltrs_fa = dir_path / "LTRs.fa"
    with open(ltrs_fa, "w") as out:
        out.write(f">{header_token}_5prime_LTR_len{ltr_len}_ext{ext_len}\n")
        out.write(five_p + "\n")
        out.write(f">{header_token}_3prime_LTR_len{ltr_len}_ext{ext_len}\n")
        out.write(three_p + "\n")

    if ARGS.verbose:
        print(f"[FAST-PATH] {dir_path.name}: extracted 5′/3′ ({ext_len} bp each) using domains TSV (LTR={ltr_len}, ext={ARGS.extension}).")

    # Align and trim (MAFFT -> trimal)
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

    # Short-circuit if trimmed alignment is too short vs raw alignment
    try:
        raw_bp = sum_ungapped_bp(aln_fa)
        clean_bp = sum_ungapped_bp(dir_path / "LTRs.aln.clean")
    except FileNotFoundError as e:
        print(f"[SKIP] Missing alignment file for {dir_path.name}: {e}")
        try:
            (dir_path / "SKIPPED.missing_alignment").write_text(str(e))
        except Exception:
            pass
        return True  # fast-path executed, even if skipped downstream


    threshold = args.min_retained_fraction * raw_bp
    if (clean_bp + ARGS.extension) <= threshold:
        log_msg(
        #print(
            f"[SKIP] Trimmed alignment too short for {dir_path.name}: "
            f"clean_bp({clean_bp}) + ext({extension}) <= "
            f"{args.min_retained_fraction:.3f} * raw_bp({raw_bp:.0f})"
        )
        try:
            (dir_path / "SKIPPED.too_trimmed").write_text(
                f"raw_bp={raw_bp}\nclean_bp={clean_bp}\nextension={extension}\n"
                f"min_retained_fraction={args.min_retained_fraction}\nthreshold={threshold}\n"
            )
        except Exception:
            pass
        return True  # fast-path executed, but result skipped


    # Fast-path: use DOMAINS_TSV value directly, do NOT subtract extension
    reported_len = ltr_len

    # Final WFA and append results (same filtering as heavy path)
    cmd_str = (
        f"{q(str(SCRIPT_DIR / 'wfa'))} -E 10000 -e 10000 -o 10000 -O 10000 -u {ARGS.mutation_rate} "
        f"{q(str(dir_path / 'LTRs.aln.clean'))}"
        " | cut -f1,5- | sed -e 's/-0\\.000000/0.000000/g' | "
        "awk -F'\t' '{if($1~/_(5prime|3prime)/) sub(/_(5prime|3prime).*$/, \"\", $1); print}' OFS='\t' | "
        "awk '$6 <= 0.35' | "
        f"awk -v P={reported_len} 'BEGIN{{OFS=\"\\t\"}} {{for(i=NF;i>=2;i--) $(i+1)=$(i); $2=P; print}}' "
    )

    if ARGS.verbose:
        print("Running:", cmd_str)
    with open(ARGS.outfile, "a") as ofile:
        subprocess.run(cmd_str, shell=True, check=True, stdout=ofile)

    return True


def process_dir(dir_path):
    """
    Process one sequence directory:
      - If domains TSV provides a safe fast-path, extract LTRs directly with extension,
        then run MAFFT -> trimal -> WFA.
      - Otherwise run the original heavy pipeline (k-mers -> map -> filter -> extract LTRs -> MAFFT -> trimal -> WFA).
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
    if ARGS.verbose:
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
    # Threshold is configurable via --min-retained-fraction (default 0.6).
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

    # Require: sum_bp(clean) + extension > (min_retained_fraction * sum_bp(raw))
    threshold = args.min_retained_fraction * raw_bp
    
    if (clean_bp + extension) <= threshold:
        log_msg(
        #print(
            f"[SKIP] Trimmed alignment too short for {dir_path.name}: "
            f"clean_bp({clean_bp}) + ext({extension}) <= "
            f"{args.min_retained_fraction:.3f} * raw_bp({raw_bp:.0f})"
        )
        try:
            (dir_path / "SKIPPED.too_trimmed").write_text(
                f"raw_bp={raw_bp}\nclean_bp={clean_bp}\nextension={extension}\n"
                f"min_retained_fraction={args.min_retained_fraction}\nthreshold={threshold}\n"
            )
        except Exception:
            pass
        return

    # Compute proposed LTR len from LTRs.aln.fa headers and adjust by extension
    proposed_len = parse_proposed_ltr_len_from_aln(aln_fa)
    if proposed_len is None:
        proposed_adj = 0
        if ARGS.verbose:
            print(f"[WARN] Could not parse proposed LTR length for {dir_path.name}; inserting 0.")
    else:
        proposed_adj = max(proposed_len - ARGS.extension, 0)


    # Final WFA and append results
    cmd_str = (
        f"{q(str(SCRIPT_DIR / 'wfa'))} -E 10000 -e 10000 -o 10000 -O 10000 -u {args.mutation_rate} "
        f"{q(str(dir_path / 'LTRs.aln.clean'))}"
        " | cut -f1,5- | sed -e 's/-0\\.000000/0.000000/g' | "
        "awk -F'\t' '{if($1~/_(5prime|3prime)/) sub(/_(5prime|3prime).*$/, \"\", $1); print}' OFS='\t' | "
        "awk '$6 <= 0.35' | "
        f"awk -v P={proposed_adj} 'BEGIN{{OFS=\"\\t\"}} {{for(i=NF;i>=2;i--) $(i+1)=$(i); $2=P; print}}' "
    )
    if verbose:
        print("Running:", cmd_str)
    with open(out_file, "a") as ofile:
        subprocess.run(cmd_str, shell=True, check=True, stdout=ofile)

def write_summary(main_out_path: str):
    """
    Read the main outfile (one row per element), sum length/transition/transversion,
    and write pooled distances to '<outfile>.summary'.

    Columns (no header) AFTER the new insertion:
    1 LTR-RT  2 proposed_len_minus_ext  3 LTR_LEN  4 substitutions
    5 transitions  6 transversions  7 p-dist  8 p-time
    9 JC69-dist  10 JC69-time  11 K2P-dist  12 K2P-time
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
            # need at least 12 columns now
            if len(parts) < 12:
                continue
            try:
                # NOTE: indexes shifted by +1 due to new column at position 2
                L = int(parts[2])             # LTR_LEN now in column 3
                ts = float(parts[4])          # transitions now in column 5
                tv = float(parts[5])          # transversions now in column 6
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
        help="Optional TSV (name\\tLTR_len). If provided, skip LTR discovery and go stright to alignment"
    )
    parser.add_argument(
        "input_fasta",
        help="Path to multi-sequence LTR-RT FASTA file."
    )
    parser.add_argument(
        "--min-retained-fraction", type=float, default=0.6, dest="min_retained_fraction",
        help="Minimum fraction of ungapped columns retained after trimming required to proceed (default: 0.6)."
    )
    parser.add_argument(
        "--assume-duplicate-same-ltr", action="store_true", dest="assume_dup_same_ltr",
        help="Override duplicate-header safety fallback in fast path. Assumes all duplicate header tokens share the same LTR length from the domains TSV. Use with caution."
    )

    args = parser.parse_args()

    # Prepare workspace
    if os.path.exists(args.temp_dir):
        shutil.rmtree(args.temp_dir)
    os.makedirs(args.temp_dir, exist_ok=True)
    open(args.outfile, "w").close()

    # Load domains mapping (if provided)
    DOMAINS = read_domains(args.domains_tsv)
    
    LOGFILE = open(args.outfile + ".log", "w")

    # Split FASTA into separate seq.fa files (safe/unique dir names; original headers preserved)
    split_fasta(args.input_fasta, args.temp_dir)

    # Set global arguments for process_dir
    ARGS = args

    # Discover sequence directories
    dirs = [str(p) for p in Path(args.temp_dir).iterdir() if p.is_dir()]

    # Process in parallel with live progress/ETA (single updating line on stderr)
    total = len(dirs)
    start = time.monotonic()
    last_update = 0.0  # throttle screen updates
    update_interval = 0.25  # seconds between redraws (tweak as desired)

    # We iterate results so we can update the counter as work completes
    completed = 0
    print("", file=sys.stderr)  # ensure stderr stream exists

    with Pool(processes=args.threads) as pool:
        # chunksize=1 gives smoother progress; increase for less overhead
        for _ in pool.imap_unordered(process_dir, dirs, chunksize=1):
            completed += 1
            now = time.monotonic()
            if (now - last_update) >= update_interval or completed == total:
                elapsed = now - start
                pct = (completed / total) * 100 if total else 100.0
                rate = (completed / elapsed) if elapsed > 0 else 0.0
                remaining = ((total - completed) / rate) if rate > 0 else float("inf")
                eta_text = _format_hms(remaining)

                # Carriage return to redraw the same line; no newline
                msg = (f"\rProcessing {completed}/{total}. "
                    f"{pct:.2f}% complete. Estimated time remaining {eta_text}")
                print(msg, end="", file=sys.stderr, flush=True)
                last_update = now

    # finish the progress line with a newline
    print("", file=sys.stderr)

    write_summary(args.outfile)

    # Cleanup
    if args.keep_temp:
        print(f"Keeping temporary directory {args.temp_dir}")
    else:
        print(f"Removing temporary directory {args.temp_dir}")
        shutil.rmtree(args.temp_dir)

    print(f"All sequences processed. Output in {args.outfile}")
    
    if LOGFILE:
        LOGFILE.close()
