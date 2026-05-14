#!/usr/bin/env python3

import argparse
import os
import shutil
import subprocess
import math
from collections import namedtuple
from multiprocessing import Pool
from pathlib import Path
import re
import shlex
import sys
import tempfile
import time
from copy import deepcopy
import platform

q = shlex.quote
MIN_SEQ_BP = 80
SCRIPT_DIR = Path(__file__).parent.resolve()

# Problems to address later:
#  1. Indel-induced offset drift. The current algorithm bins LTR-LTR pairs at 24 bp offset resolution and keeps a single bin. Indels in the LTR-vs-LTR alignment shift downstream pairs' offsets out of that
#  bin, so middle-LTR pairs get filtered (the dumbbell pattern). Heavily diverged LTR-RTs lose enough pairs to fail qualification entirely.
#  2. Silent flank contamination. The boundary call anchors the extract at positions 0 and seq_len, assuming the input IS the LTR-RT. Flanks ≤ 100 bp pass the termini-reach test and get silently glued onto
#   the extracts, contaminating downstream MAFFT alignment and inflating divergence estimates.
#  Both share a root cause: the algorithm searches for the LTR in 1D (an offset bucket) and assumes a property of the input (clean LTR-RT, no flank) that may not hold.
#  The reframe.
#  The LTR-LTR signal is a 2D geometric feature — pairs forming a collinear line segment in (s1, s2) space with slope ≈ 1. Non-LTR noise is scattered with no collinear structure. Detect the segment, and
#  three useful numbers fall out of its geometry:
#  - LTR boundaries within the input = the segment's endpoints in (s1, s2).
#  - Inter-LTR spacing = the segment's offset (s2 − s1), now a measured property rather than a search parameter.
#  - Flank size = the gap between the segment's endpoints and the plot's corners.
#  Indels become small steps in offset along the segment — still collinear, still continuous in (s1, s2), still part of the same LTR signal. Overextension becomes a visible inset of the segment from the
#  plot's edges — the input is no longer assumed to BE the LTR-RT; it's a search space within which the LTR-RT is geometrically located.
#  Formally: line-segment detection on a sparse 2D point cloud — same problem family as Hough-transform line detection, RANSAC line fitting, or collinear seed chaining in long-read aligners (minimap2,
#  LASTZ).
#  We need a homology detector that doesn't presume where the homologous regions sit within the input. Find the longest contiguous chain of evidence that "these two regions are duplicates of each
#  other"; report the chain's endpoints as the LTR boundaries, the chain's offset(s) as the inter-LTR spacing, and the gaps between chain endpoints and input edges as flank size. Three useful numbers, one
#  geometric measurement, no input assumptions.

# ---------------- WFA binary selection with fallback ---------------- #
# Tries wfa_linux1 first, then wfa_linux2. If both fail, instruct compile.

def _test_wfa(bin_path: Path) -> bool:
    """Return True if bin_path exists and can execute a trivial --help command."""
    if not bin_path.exists():
        return False
    try:
        # Use --help or -h depending on the binary's behavior; stderr ignored.
        subprocess.run(
            [str(bin_path), "-h"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True
        )
        return True
    except Exception:
        return False

if platform.system() == "Darwin":
    WFA_BIN = SCRIPT_DIR / "wfa_mac"
    if not _test_wfa(WFA_BIN):
        print("[ERROR] wfa_mac does not execute correctly. "
              "Please compile wfa.cpp for macOS.", file=sys.stderr)
        sys.exit(1)

else:
    # Linux: try two possible binaries
    WFA1 = SCRIPT_DIR / "wfa_linux1"
    WFA2 = SCRIPT_DIR / "wfa_linux2"

    if _test_wfa(WFA1):
        WFA_BIN = WFA1
        print(f"[INFO] Using WFA binary: {WFA1}")
    elif _test_wfa(WFA2):
        WFA_BIN = WFA2
        print(f"[INFO] Primary WFA failed; using fallback: {WFA2}")
    else:
        print(
            "[FATAL] Neither wfa_linux1 nor wfa_linux2 executed successfully.\n"
            "You may need to compile Kmer2LTR/wfa.cpp locally for your system.\n"
            "To compile manually:\n"
            "    g++ -O3 -o wfa_linux1 wfa.cpp\n",
            file=sys.stderr
        )
        sys.exit(1)
# -------------------------------------------------------------------- #



# Global args/domains placeholders
ARGS = None
DOMAINS = None  # dict[str,int] or None
LOGFILE = None
LOG_PATH = None  # <- add this


# ---------------- Debug-mode helpers ---------------- #
# These are all no-ops when ARGS.debug is False, so they're safe to sprinkle
# liberally through the pipeline without affecting normal-mode performance.

def _debug_enabled() -> bool:
    """True if --debug was set. Safe to call before ARGS exists (e.g. before _pool_init)."""
    return bool(getattr(ARGS, "debug", False)) if ARGS is not None else False


def _ensure_debug_dir(dir_path) -> "Path | None":
    """Return <dir_path>/debug, creating it. Returns None when debug is off."""
    if not _debug_enabled():
        return None
    debug_dir = Path(dir_path) / "debug"
    debug_dir.mkdir(parents=True, exist_ok=True)
    return debug_dir


def _debug_write(debug_dir, filename: str, content: str):
    """Write `content` to debug_dir/filename. No-op if debug_dir is None."""
    if debug_dir is None:
        return
    path = Path(debug_dir) / filename
    if not content.endswith("\n"):
        content += "\n"
    path.write_text(content)


def _debug_log(debug_dir, line: str):
    """Append a line to debug_dir/README.txt. No-op if debug_dir is None."""
    if debug_dir is None:
        return
    with open(Path(debug_dir) / "README.txt", "a") as fh:
        fh.write(line.rstrip() + "\n")


def _debug_print(msg: str):
    """Print a [DEBUG] message to stderr if debug is on. No-op otherwise."""
    if _debug_enabled():
        print(f"[DEBUG] {msg}", file=sys.stderr, flush=True)


_DEBUG_BANNER = """\
================================================================
DEBUG MODE ENABLED  (Kmer2LTR --debug)
================================================================
What --debug does:
  * Forces --verbose (commands echoed before running).
  * Forces -k / --keep_temp (no temp cleanup).
  * Adds the wfa.cpp '-d' flag to the WFA divergence call so the
    alignment view and -k/-K boundary kmers are written to a file.
  * Writes per-sequence intermediates to <temp_dir>/<seq>/debug/
    so you can diff alignments, inspect kmer pair sets, and trace
    where divergence calls came from.
  * --debug is incompatible with --purge-subdirs (which deletes
    sub-temp dirs as we go).

Pipeline at a glance (per LTR-RT in the input FASTA):
  1. Boundary discovery
       Fast path  (only if --domains TSV provides this header's LTR_len):
         extract first/last (LTR_len + --extension) bp directly.
       Otherwise  (kmer pipeline):
         a) Enumerate kmers length --kmin..--kmax that appear in BOTH
            halves of the sequence, --dist or more bp apart.
         b) Sort by 5' position; keep the longest-increasing-subsequence
            of the 3' positions => collinear pair set.
         c) 98th-percentile of 5'-end => end of 5' LTR;
            2nd-percentile of 3'-start => start of 3' LTR;
            symmetrize lengths; extend by --extension bp on each side.
  2. Pairwise alignment of 5' vs 3' LTR (MAFFT --auto, or WFA if --wfa-align).
  3. trimal -automated1 -keepheader on the alignment.
  4. Re-run WFA on the trimmed alignment for divergence stats.
       WFA uses -k 5 -K: requires 5 consecutive matches on each side
       to enter/exit the reliable region, and EXCLUDES those k columns
       from the metrics. This is how unreliable LTR-RT termini are
       chopped before counting substitutions.
  5. (Optional) Build IUPAC consensus per LTR-RT for --ltr-cluster.

Per-sequence debug files (under <temp_dir>/<sanitized_name>/debug/):
  Fast path:
    01_fastpath_extracted.fa     5'/3' LTR FASTA cut using domains LTR_len.
  Full pipeline:
    01_kmer_pairs.tsv            All candidate kmer pairs (s1,e1,s2,e2).
    01b_diagonal_pairs.tsv       Pairs surviving the LTR-RT diagonal filter
                                 (chosen offset annotated in the header).
    02_lis_pairs.tsv             Pairs that survive the LIS collinearity filter.
    03_ltr_boundaries.txt        Final 5'/3' coordinates + intermediate maths.
    03b_dotplot.png              Self-vs-self dot plot of the kmer pairs at
                                 each filtering stage (gray=all, blue=diagonal
                                 survivors, orange=LIS survivors), with the
                                 chosen diagonal and final boundaries marked.
                                 Also written when the LTR-RT is REJECTED at
                                 the diagonal or LIS step (status_note in plot).
    04_extracted_ltrs.fa         5'/3' LTR FASTA from kmer discovery.
  Both paths:
    05_alignment.fa              Pairwise alignment (MAFFT or WFA).
    06_trimal.fa                 trimal -automated1 output.
    07_wfa_cmd.txt               Exact WFA divergence command (with -d).
    08_wfa_stdout.txt            WFA divergence raw stats line.
    09_wfa_stderr.txt            WFA -d alignment view + boundary kmers.
    10_consensus_realign.fa      Realigned ungapped LTRs (only --ltr-cluster).
    11_consensus.fa              IUPAC consensus (only --ltr-cluster).
    README.txt                   Step-by-step narrative for this LTR-RT.

Things --debug is designed to help you diagnose:
  * "Why does providing -D/--domains give a very different result?"
      Same LTR-RT processed without -D writes 02_lis_pairs.tsv +
      03_ltr_boundaries.txt and 04_extracted_ltrs.fa. With -D you get
      01_fastpath_extracted.fa instead. Diff the extracted FASTAs —
      different coordinates feed into completely different alignments.
  * "Why does --wfa-align disagree with MAFFT?"
      Diff 05_alignment.fa across the two runs. trimal + divergence
      maths are deterministic given the alignment.
  * "Is wfa.cpp -k actually doing what I think?"
      09_wfa_stderr.txt shows the boundary kmers and a fallback note
      ('-k N could not find N consecutive matches ... using full
      alignment') when applicable.
  * "Why was this LTR-RT skipped?"
      Skip markers (SKIPPED.too_short / SKIPPED.no_filtered_pairs)
      land in the sequence dir; reasons are appended to README.txt.
================================================================
"""


def print_debug_banner():
    print(_DEBUG_BANNER, file=sys.stderr, flush=True)


def _plot_dotplot(debug_dir, seq_len, all_pairs, diagonal_pairs, lis_pairs,
                  chosen_offset=None, final_end1=None, final_start2=None,
                  status_note=None):
    """
    Self-vs-self dot plot of kmer pairs at each filtering stage.
    Written to <debug_dir>/03b_dotplot.png. Silent no-op on missing matplotlib.

    Each pair (s1, s2) is one dot at (x=s1, y=s2) -- the standard self dot-plot
    convention. The true 5'LTR<->3'LTR diagonal is the line y = x + (seq_len - L).
    Layering: lightgray (all candidates) under blue (diagonal-filter survivors)
    under orange (LIS survivors) so each filter stage is visible.
    """
    if debug_dir is None:
        return
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return

    fig, ax = plt.subplots(figsize=(6, 6))

    half = seq_len // 2
    ax.axvline(half, color="0.85", linestyle=":", linewidth=0.5, zorder=0)
    ax.axhline(half, color="0.85", linestyle=":", linewidth=0.5, zorder=0)
    ax.plot([0, seq_len], [0, seq_len], color="0.9", linewidth=0.5, zorder=0)

    # Layered markers chosen so all three sets stay visible even when LIS keeps
    # nearly all diagonal pairs: gray = small filled background dot; blue = open
    # ring (no fill) so an underlying orange dot still shows through; orange =
    # small filled dot that sits inside the blue ring.
    if all_pairs:
        ax.scatter([p[0] for p in all_pairs], [p[2] for p in all_pairs],
                   c="lightgray", s=4, alpha=0.7, linewidths=0,
                   label=f"all pairs (n={len(all_pairs)})")
    if diagonal_pairs:
        ax.scatter([p[0] for p in diagonal_pairs], [p[2] for p in diagonal_pairs],
                   facecolors="none", edgecolors="#0173B2", s=20, alpha=0.85,
                   linewidths=0.6,
                   label=f"diagonal survivors (n={len(diagonal_pairs)})")
    if lis_pairs:
        ax.scatter([p[0] for p in lis_pairs], [p[2] for p in lis_pairs],
                   c="#D55E00", s=6, alpha=0.95, linewidths=0,
                   label=f"LIS survivors (n={len(lis_pairs)})")

    if chosen_offset is not None:
        x_end = seq_len - chosen_offset
        ax.plot([0, x_end], [chosen_offset, seq_len],
                color="#009E73", linewidth=0.8, alpha=0.7,
                label=f"chosen diagonal (offset={chosen_offset})")

    if final_end1 is not None:
        ax.axvline(final_end1, color="black", linestyle="--", linewidth=0.8,
                   label=f"5' LTR end = {final_end1}")
    if final_start2 is not None:
        ax.axhline(final_start2, color="black", linestyle="-.", linewidth=0.8,
                   label=f"3' LTR start = {final_start2}")

    ax.set_xlabel("position in LTR-RT (bp)")
    ax.set_ylabel("position in LTR-RT (bp)")
    ax.set_xlim(0, seq_len)
    ax.set_ylim(0, seq_len)
    ax.set_aspect("equal")
    # Lower-right is the data-free triangle: pairs only exist where
    # s1 < half_len < s2 (upper-left), so the legend never overlaps dots,
    # the LTR-LTR diagonal, or the boundary lines.
    ax.legend(loc="lower right", fontsize=7, framealpha=0.9)
    ax.tick_params(labelsize=8)

    if status_note:
        ax.text(0.98, 0.02, status_note, transform=ax.transAxes,
                ha="right", va="bottom", fontsize=8, color="0.3",
                family="monospace")

    fig.tight_layout()
    fig.savefig(Path(debug_dir) / "03b_dotplot.png", dpi=200)
    plt.close(fig)


# -------------- end debug helpers -------------- #


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
            total += len(s) - s.count("-")
    return total


def _sum_ungapped_bp_str(data: str) -> int:
    """Count total ungapped basepairs from a FASTA string."""
    total = 0
    for line in data.splitlines():
        if not line or line.startswith(">"):
            continue
        s = line.strip()
        total += len(s) - s.count("-")
    return total


def _read_aln_stats(aln_fa: Path):
    """Read proposed LTR length from header and count ungapped bp in one pass."""
    proposed_len = None
    total_bp = 0
    with open(aln_fa) as fh:
        for line in fh:
            if line.startswith(">"):
                if proposed_len is None:
                    m = re.search(r"5p:1-(\d+)", line)
                    if m:
                        proposed_len = int(m.group(1))
            else:
                s = line.strip()
                if s:
                    total_bp += len(s) - s.count("-")
    return proposed_len, total_bp


def _read_aln_stats_str(data: str):
    """Extract proposed LTR length and count ungapped bp from a FASTA string."""
    proposed_len = None
    total_bp = 0
    for line in data.splitlines():
        if line.startswith(">"):
            if proposed_len is None:
                m = re.search(r"5p:1-(\d+)", line)
                if m:
                    proposed_len = int(m.group(1))
        else:
            s = line.strip()
            if s:
                total_bp += len(s) - s.count("-")
    return proposed_len, total_bp


# ---- In-memory alignment/trimal helpers ---- #

def _run_mafft_mem(fasta_str: str, verbose: bool) -> str:
    """Align two sequences via mafft. Uses /dev/shm if available (Linux), otherwise system temp, captures aligned output in memory."""
    # Use shared memory if available (Linux performance), fallback to system temp (macOS compatibility)
    temp_dir = '/dev/shm' if os.path.exists('/dev/shm') else None
    fd, tmp_path = tempfile.mkstemp(suffix='.fa', dir=temp_dir)
    try:
        with os.fdopen(fd, 'w') as f:
            f.write(fasta_str)
        cmd = ["mafft", "--auto", "--thread", "1", tmp_path]
        if verbose:
            print("Running:", " ".join(cmd))
        result = subprocess.run(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.DEVNULL, text=True)
        if result.returncode != 0:
            raise subprocess.CalledProcessError(result.returncode, cmd)
        return result.stdout
    finally:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass


def _run_wfa_align_mem(fasta_str: str, verbose: bool) -> str:
    """Align two sequences via WFA in FASTA-output mode. Reads from stdin, returns aligned FASTA string."""
    cmd = [str(WFA_BIN), "-F", "-"]
    if verbose:
        print("Running:", " ".join(cmd), "(stdin)")
    result = subprocess.run(cmd, input=fasta_str, capture_output=True, text=True)
    if result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, cmd)
    return result.stdout


def _run_align_mem(fasta_str: str, verbose: bool) -> str:
    """Dispatch to WFA or mafft based on ARGS.wfa_align."""
    if ARGS.wfa_align:
        return _run_wfa_align_mem(fasta_str, verbose)
    return _run_mafft_mem(fasta_str, verbose)


def _run_trimal_mem(aln_str: str, verbose: bool) -> str:
    """Run trimal on an alignment string. Uses /dev/shm if available (Linux), otherwise system temp (trimal needs seekable input)."""
    # Use shared memory if available (Linux performance), fallback to system temp (macOS compatibility)
    temp_dir = '/dev/shm' if os.path.exists('/dev/shm') else None
    fd, tmp_path = tempfile.mkstemp(suffix='.fa', dir=temp_dir)
    try:
        with os.fdopen(fd, 'w') as f:
            f.write(aln_str)
        cmd = ["trimal", "-automated1", "-keepheader", "-in", tmp_path]
        if verbose:
            print("Running:", " ".join(cmd))
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return result.stdout
    finally:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass


# ---- IUPAC consensus helpers (for optional --ltr-cluster output) ---- #

# Two-base IUPAC ambiguity codes.
_IUPAC_PAIR = {
    frozenset(['A', 'G']): 'R',  # purine
    frozenset(['C', 'T']): 'Y',  # pyrimidine
    frozenset(['G', 'C']): 'S',  # strong (3 H-bonds)
    frozenset(['A', 'T']): 'W',  # weak (2 H-bonds)
    frozenset(['G', 'T']): 'K',  # keto
    frozenset(['A', 'C']): 'M',  # amino
}


def _strip_gaps_fasta(aln_str: str) -> str:
    """Remove '-' gaps from sequence lines of a FASTA string; preserve headers."""
    out_lines = []
    for line in aln_str.splitlines():
        if line.startswith('>'):
            out_lines.append(line)
        else:
            out_lines.append(line.replace('-', ''))
    return '\n'.join(out_lines) + '\n'


def _parse_fasta_records(fasta_str: str):
    """Yield (header, sequence) pairs from a FASTA string."""
    header = None
    seq_parts = []
    for line in fasta_str.splitlines():
        if not line:
            continue
        if line.startswith('>'):
            if header is not None:
                yield header, ''.join(seq_parts)
            header = line[1:].strip()
            seq_parts = []
        else:
            seq_parts.append(line.strip())
    if header is not None:
        yield header, ''.join(seq_parts)


def _build_iupac_consensus_pairwise(aligned_fasta_str: str):
    """
    Build an IUPAC consensus from a 2-sequence pairwise alignment string.

    Per-column rules (uppercased):
      - both ACGT and equal -> that base
      - both ACGT and different -> 2-way IUPAC ambiguity (R/Y/S/W/K/M)
      - exactly one '-' -> the non-gap base (or 'N' if non-ACGT)
      - both '-' -> column skipped
      - any non-ACGT, non-gap base -> 'N'

    Returns the consensus string, or None if the input cannot be parsed as
    exactly two equal-length non-empty sequences.
    """
    records = list(_parse_fasta_records(aligned_fasta_str))
    if len(records) != 2:
        return None
    s1 = records[0][1].upper()
    s2 = records[1][1].upper()
    if not s1 or len(s1) != len(s2):
        return None

    bases = {'A', 'C', 'G', 'T'}
    out = []
    for a, b in zip(s1, s2):
        if a == '-' and b == '-':
            continue
        if a == '-':
            out.append(b if b in bases else 'N')
            continue
        if b == '-':
            out.append(a if a in bases else 'N')
            continue
        if a == b:
            out.append(a if a in bases else 'N')
            continue
        if a in bases and b in bases:
            out.append(_IUPAC_PAIR.get(frozenset([a, b]), 'N'))
        else:
            out.append('N')
    return ''.join(out)


def _build_consensus_from_trimmed(clean_data: str, full_header: str, verbose: bool, debug_dir=None):
    """
    Build an IUPAC consensus FASTA record from a trimal-cleaned 2-sequence alignment.
    Pipeline: strip gaps, re-align with WFA, compute IUPAC consensus.
    Returns the FASTA record (with trailing newline), or None on failure.

    debug_dir: optional Path; when set, retain
        10_consensus_realign.fa  (WFA realignment of the ungapped trimmed LTRs)
    """
    ungapped = _strip_gaps_fasta(clean_data)
    try:
        wfa_aln = _run_wfa_align_mem(ungapped, verbose)
    except subprocess.CalledProcessError as exc:
        log_msg(f"[CONSENSUS-SKIP] WFA realignment failed for "
                f"{full_header.split()[0] if full_header else '?'}: exit {exc.returncode}")
        if debug_dir is not None:
            _debug_log(debug_dir,
                       f"[CONSENSUS] WFA realignment failed (exit {exc.returncode}); no consensus written.")
        return None
    if debug_dir is not None:
        _debug_write(debug_dir, "10_consensus_realign.fa", wfa_aln)
    cons = _build_iupac_consensus_pairwise(wfa_aln)
    if not cons:
        log_msg(f"[CONSENSUS-SKIP] Empty/invalid consensus for "
                f"{full_header.split()[0] if full_header else '?'}")
        if debug_dir is not None:
            _debug_log(debug_dir,
                       "[CONSENSUS] _build_iupac_consensus_pairwise returned None (parsing/length mismatch).")
        return None
    return f">{full_header}\n{cons}\n"


def _format_hms(seconds: float) -> str:
    """Return 'Hh:Mm:Ss' for a nonnegative seconds value."""
    if not math.isfinite(seconds) or seconds < 0:
        return "0h:0m:0s"
    s = int(round(seconds))
    h, s = divmod(s, 3600)
    m, s = divmod(s, 60)
    return f"{h}h:{m}m:{s}s"
    
def _pool_init(args_ns, domains_map, log_path):
    """Initializer for multiprocessing workers (macOS spawn-safe)."""
    global ARGS, DOMAINS, LOGFILE, LOG_PATH
    ARGS = args_ns
    DOMAINS = domains_map
    LOG_PATH = log_path
    LOGFILE = None  # don't share a handle across processes

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
    global LOGFILE, LOG_PATH
    # Prefer path-based logging to avoid sharing file handles across processes.
    if LOG_PATH:
        try:
            with open(LOG_PATH, "a") as lf:
                lf.write(msg.strip() + "\n")
        except Exception:
            pass
    elif LOGFILE:
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


FastaEntry = namedtuple("FastaEntry", ["safe_name", "header_token", "file_offset", "length_bytes", "is_duplicate"])


def index_fasta(input_fasta, temp_dir):
    """
    Scan the input FASTA in a single pass WITHOUT writing any seq.fa files.

    Returns:
        entries: list[FastaEntry] — one per sequence, recording file offset/length
                 so each sequence can be materialized later via write_batch().

    Side effects:
        Writes header_map.tsv into temp_dir (needed for --reuse-existing filtering).
    """
    temp_dir = Path(temp_dir)
    map_path = temp_dir / "header_map.tsv"
    counts = {}       # header_token -> occurrence count
    entries = []
    dup_tokens = set() # tokens that appear more than once

    with open(input_fasta, "rb") as infile, open(map_path, "w") as map_out:
        current_offset = None
        current_safe = None
        current_token = None

        for line in infile:
            if line.startswith(b">"):
                # Finalize previous entry
                if current_offset is not None:
                    length = infile.tell() - len(line) - current_offset
                    entries.append(FastaEntry(
                        safe_name=current_safe,
                        header_token=current_token,
                        file_offset=current_offset,
                        length_bytes=length,
                        is_duplicate=False,  # placeholder, resolved below
                    ))

                header_line = line[1:].strip().decode("utf-8", errors="replace")
                header_token = header_line.split()[0] if header_line else ""
                counts[header_token] = counts.get(header_token, 0) + 1

                base_safe = sanitize_name(header_token)
                if counts[header_token] == 1:
                    safe = base_safe
                else:
                    safe = f"{base_safe}__dup{counts[header_token]}"

                if counts[header_token] == 2:
                    dup_tokens.add(header_token)

                map_out.write(f"{safe}\t{header_token}\n")

                current_offset = infile.tell() - len(line)
                current_safe = safe
                current_token = header_token

        # Finalize the last entry
        if current_offset is not None:
            length = infile.tell() - current_offset
            entries.append(FastaEntry(
                safe_name=current_safe,
                header_token=current_token,
                file_offset=current_offset,
                length_bytes=length,
                is_duplicate=False,
            ))

    # Mark duplicate entries
    if dup_tokens:
        entries = [
            e._replace(is_duplicate=True) if e.header_token in dup_tokens else e
            for e in entries
        ]

    return entries


def write_batch(entries, input_fasta, temp_dir):
    """
    Given a slice of FastaEntry items, seek into the source FASTA and write
    each entry's seq.fa into its temp subdir. Also creates DUPLICATE_NAME
    marker files where needed.
    """
    temp_dir = Path(temp_dir)
    with open(input_fasta, "rb") as fh:
        for entry in entries:
            dir_path = temp_dir / entry.safe_name
            dir_path.mkdir(parents=True, exist_ok=True)

            # Write seq.fa
            fh.seek(entry.file_offset)
            data = fh.read(entry.length_bytes)
            with open(dir_path / "seq.fa", "wb") as out:
                out.write(data)

            # Write DUPLICATE_NAME marker if needed
            if entry.is_duplicate:
                try:
                    (dir_path / "DUPLICATE_NAME").write_text("")
                except Exception:
                    pass


def read_seq_string(seq_fa: Path) -> str:
    """Read a single FASTA sequence (concatenate all non-header lines)."""
    seq_chunks = []
    with open(seq_fa) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            seq_chunks.append(line.strip())
    return "".join(seq_chunks)


def _read_fasta_entry(seq_fa: Path):
    """Read header text and sequence string from a single-sequence FASTA in one pass."""
    header = ""
    seq_chunks = []
    with open(seq_fa) as fh:
        for line in fh:
            if line.startswith(">"):
                header = line[1:].strip()
            else:
                seq_chunks.append(line.strip())
    return header, "".join(seq_chunks)


def _lis_indices(values):
    """
    Longest Increasing Subsequence (non-strict) via patience sort.
    Returns the indices into `values` that form the LIS.
    O(n log n) time.
    """
    import bisect
    n = len(values)
    if n == 0:
        return []

    # tails[i] = smallest tail value for an increasing subsequence of length i+1
    tails = []
    # For each position, which LIS length it extends
    lis_len_at = [0] * n
    # For backtracking: predecessor index
    predecessor = [-1] * n
    # tail_idx[i] = index in values[] of tails[i]
    tail_idx = []

    for i in range(n):
        v = values[i]
        # bisect_right for non-strict (>=) increasing
        pos = bisect.bisect_right(tails, v)
        if pos == len(tails):
            tails.append(v)
            tail_idx.append(i)
        else:
            tails[pos] = v
            tail_idx[pos] = i
        lis_len_at[i] = pos
        predecessor[i] = tail_idx[pos - 1] if pos > 0 else -1

    # Backtrack to recover the actual subsequence
    lis_length = len(tails)
    result = []
    idx = tail_idx[lis_length - 1]
    for _ in range(lis_length):
        result.append(idx)
        idx = predecessor[idx]
    result.reverse()
    return result


# Maximum number of occurrences for a kmer to be considered informative.
# Above this, the kmer is likely low-complexity / microsatellite.
_MAX_KMER_OCCURRENCES = 10

# Minimum number of collinear pairs after LIS to trust the boundary call.
_MIN_LIS_PAIRS = 20

# Diagonal-filter defaults. bin_width must be wide enough to absorb
# divergence-induced offset jitter on a single diagonal -- a few-bp indel in
# the LTR alignment shifts s2-s1 by that many bp, so a single true diagonal
# typically spans 5-15 bp of offset. Setting bin_width too small splits the
# true diagonal across adjacent bins, neither of which then passes the
# termini-reach test individually. 24 bp gives comfortable headroom while
# still cleanly separating diagonals from spurious neighbors.
#
# The termini-reach test uses the SMALLER of (frac * seq_len, hard-bp-cap).
# For short inputs (~500 bp) the fractional gap dominates (46 bp at frac=0.10);
# for long inputs (~15 kb) the fractional gap (1.5 kb) is much larger than
# typical LTRs (100-2000 bp), letting LTR-internal cross-match diagonals pass
# as "reaching the terminus" -- the bp cap stops that. 100 bp is generous
# enough to absorb boundary kmer dropout from divergence/indels at LTR ends.
_DIAG_BIN_WIDTH = 24
_DIAG_MIN_PAIRS_PER_BIN = 5
_DIAG_TERMINUS_GAP_FRAC = 0.10
_DIAG_TERMINUS_GAP_MAX_BP = 100


def _pick_ltr_diagonal(all_mapped, seq_len,
                       bin_width=_DIAG_BIN_WIDTH,
                       min_pairs_per_bin=_DIAG_MIN_PAIRS_PER_BIN,
                       max_terminus_gap_frac=_DIAG_TERMINUS_GAP_FRAC,
                       max_terminus_gap_bp=_DIAG_TERMINUS_GAP_MAX_BP):
    """
    Pick the kmer-pair subset on the LTR-RT 5'LTR<->3'LTR dot-plot diagonal.

    Why: LIS on s2 alone finds the LONGEST monotonic kmer-pair set, which on
    internally repetitive LTR-RTs is often a spurious diagonal (internal-region
    repeat phase) rather than the true LTR<->LTR diagonal. Restricting to a
    single diagonal before LIS fixes the boundary call.

    Selection:
      1. Histogram offsets (s2 - s1) at `bin_width` resolution.
      2. Keep bins with >= `min_pairs_per_bin` pairs AND with some pair whose
         s1 is within `max_terminus_gap_frac * seq_len` of the 5' end AND some
         pair whose s2 is within the same fraction of the 3' end. This rejects
         internal-region repeat diagonals (which sit in the middle of the
         input and don't touch the termini).
      3. Reject bins whose implied LTR length exceeds seq_len/2 (the LTRs
         can't overlap -- this rejects internal-overlap diagonals that the
         downstream pipeline would have to cap to seq_len/2 anyway).
      4. Among survivors, pick the diagonal with the highest density-weighted
         score: pairs^2 / implied_LTR_len. Equivalent to (pair count) x
         (pair density along the diagonal). Balances raw count and
         concentration: a fat, sparse "diagonal" (often a shift artifact
         or LTR-internal cross-match that spans more of the input than the
         true LTR) loses to a dense, narrow LTR-LTR diagonal of the right
         width. The TRUE LTR<->LTR diagonal has pairs at every position
         within the LTR (high density) AND covers the LTR length
         (substantial count) -- it dominates this score among physically
         valid candidates.

    Returns (keep_indices, chosen_offset_bin) or None if no diagonal qualifies.
    Indices are positions in `all_mapped`.
    """
    if not all_mapped:
        return None

    by_bin = {}
    for idx, p in enumerate(all_mapped):
        s1, _, s2, _ = p
        b = ((s2 - s1) // bin_width) * bin_width
        by_bin.setdefault(b, []).append(idx)

    gap = min(max_terminus_gap_frac * seq_len, max_terminus_gap_bp)
    near_start = gap
    near_end = seq_len - gap

    candidates = []
    for b, indices in by_bin.items():
        if len(indices) < min_pairs_per_bin:
            continue
        implied_ltr = seq_len - b
        # Non-overlap: the two LTRs together must fit in the input.
        if 2 * implied_ltr > seq_len:
            continue
        s1_min = min(all_mapped[i][0] for i in indices)
        s1_max = max(all_mapped[i][1] for i in indices)
        s2_max = max(all_mapped[i][3] for i in indices)
        if s1_min > near_start or s2_max < near_end:
            continue
        n = len(indices)
        s1_span = s1_max - s1_min
        score = (n * n) / implied_ltr  # pairs^2 / implied_LTR
        candidates.append((score, n, s1_span, b, indices))

    if not candidates:
        return None

    # Highest density-weighted score wins; pair count and span as tiebreakers.
    candidates.sort(key=lambda c: (c[0], c[1], c[2]), reverse=True)
    best_bin = candidates[0][3]

    lo = best_bin
    hi = best_bin + bin_width
    keep = [i for i, p in enumerate(all_mapped) if lo <= (p[2] - p[0]) < hi]
    return keep, best_bin


def _discover_and_extract_ltrs(seq, kmin, kmax, dist, std_factor, extension, debug_dir=None):
    """
    In-process replacement for the jellyfish + map_kmers + filter_kmers + extract_ltrs
    subprocess chain.  Returns (ltr_5p, ltr_3p, end5p, start3p) or None.
    All coordinates are 1-based inclusive.

    debug_dir: optional Path; when set, write
        01_kmer_pairs.tsv   (all candidate pairs incl. kmer size)
        02_lis_pairs.tsv    (pairs that survived LIS)
        03_ltr_boundaries.txt (raw percentile picks, symmetrization, extension)
    """
    seq_upper = seq.upper()
    seq_len = len(seq_upper)
    half_len = seq_len // 2

    # --- Kmer discovery + position mapping (replaces jellyfish + map_kmers) ---
    # Kmers appearing up to _MAX_KMER_OCCURRENCES times are kept; all valid
    # cross-half pairs are generated.  LIS resolves ambiguity downstream.
    all_mapped = []  # (start1, end1, start2, end2) all 1-based
    all_mapped_k = []  # parallel list of kmer sizes for debug output only
    for k in range(kmin, kmax + 1):
        kmer_positions = {}
        for i in range(seq_len - k + 1):
            kmer = seq_upper[i:i + k]
            if kmer in kmer_positions:
                kmer_positions[kmer].append(i)
            else:
                kmer_positions[kmer] = [i]

        for positions in kmer_positions.values():
            if len(positions) < 2 or len(positions) > _MAX_KMER_OCCURRENCES:
                continue
            # Generate all valid (first-half, second-half) pairs
            first_half = [p for p in positions if p + 1 <= half_len]
            second_half = [p for p in positions if p + 1 > half_len]
            for pos1 in first_half:
                for pos2 in second_half:
                    if pos2 - pos1 < dist:
                        continue
                    s1 = pos1 + 1  # 1-based
                    s2 = pos2 + 1
                    all_mapped.append((s1, pos1 + k, s2, pos2 + k))
                    if debug_dir is not None:
                        all_mapped_k.append(k)

    if debug_dir is not None:
        rows = ["# all cross-half kmer pairs (1-based, inclusive)\n"
                "# columns: kmer_size, start_5p, end_5p, start_3p, end_3p"]
        for (s1, e1, s2, e2), k in zip(all_mapped, all_mapped_k):
            rows.append(f"{k}\t{s1}\t{e1}\t{s2}\t{e2}")
        _debug_write(debug_dir, "01_kmer_pairs.tsv", "\n".join(rows))

    if not all_mapped:
        if debug_dir is not None:
            _debug_log(debug_dir,
                       "[FULL] No kmer pairs satisfy the cross-half + --dist constraint; returning None.")
        return None

    # Snapshot the full pair set for the debug dot plot before the diagonal
    # filter overwrites all_mapped with only its survivors.
    all_pairs_snapshot = list(all_mapped) if debug_dir is not None else None

    # --- LTR-RT diagonal filter (run BEFORE LIS) ---
    # True 5'LTR<->3'LTR matches lie on a single dot-plot diagonal with
    # offset s2 - s1 ~ seq_len - LTR_len. Internally repetitive LTR-RTs can
    # produce stronger but spurious diagonals (internal-internal repeat phase,
    # LTR-internal coincidence) at SMALLER offsets. Restrict to the diagonal
    # with the largest offset that reaches both termini; LIS then picks the
    # collinear subset within that diagonal.
    diag = _pick_ltr_diagonal(all_mapped, seq_len)
    if diag is None:
        if debug_dir is not None:
            _debug_log(debug_dir,
                       f"[FULL] No kmer-pair diagonal reaches both termini "
                       f"(within min({int(_DIAG_TERMINUS_GAP_FRAC*100)}%, "
                       f"{_DIAG_TERMINUS_GAP_MAX_BP} bp) of each end) "
                       f"with >= {_DIAG_MIN_PAIRS_PER_BIN} pairs; returning None.")
            _plot_dotplot(debug_dir, seq_len, all_pairs_snapshot, None, None,
                          status_note="REJECTED: no diagonal qualified")
        return None
    keep_idx, chosen_offset = diag
    all_mapped = [all_mapped[i] for i in keep_idx]
    if debug_dir is not None:
        all_mapped_k = [all_mapped_k[i] for i in keep_idx]
        rows = [
            f"# kmer pairs surviving the LTR-RT diagonal filter\n"
            f"# chosen offset bin: {chosen_offset} (+/- {_DIAG_BIN_WIDTH}); "
            f"seq_len={seq_len}; implied LTR_len ~ seq_len - offset = "
            f"{seq_len - chosen_offset}\n"
            "# columns: kmer_size, start_5p, end_5p, start_3p, end_3p"
        ]
        for (s1, e1, s2, e2), k in zip(all_mapped, all_mapped_k):
            rows.append(f"{k}\t{s1}\t{e1}\t{s2}\t{e2}")
        _debug_write(debug_dir, "01b_diagonal_pairs.tsv", "\n".join(rows))

    # --- Longest Increasing Subsequence on s2 (replaces greedy filter) ---
    # Sort by s1, break ties by s2 so overlapping kmer sizes don't collide.
    if debug_dir is not None:
        # Sort the parallel k list along with all_mapped to preserve the mapping.
        combined = sorted(zip(all_mapped, all_mapped_k), key=lambda x: (x[0][0], x[0][2]))
        all_mapped = [c[0] for c in combined]
        all_mapped_k = [c[1] for c in combined]
    else:
        all_mapped.sort(key=lambda x: (x[0], x[2]))
    s2_values = [row[2] for row in all_mapped]
    lis_idx = _lis_indices(s2_values)
    monotonic = [all_mapped[i] for i in lis_idx]

    if debug_dir is not None:
        rows = ["# kmer pairs surviving the LIS collinearity filter\n"
                "# columns: kmer_size, start_5p, end_5p, start_3p, end_3p"]
        for i in lis_idx:
            s1, e1, s2, e2 = all_mapped[i]
            rows.append(f"{all_mapped_k[i]}\t{s1}\t{e1}\t{s2}\t{e2}")
        _debug_write(debug_dir, "02_lis_pairs.tsv", "\n".join(rows))

    if len(monotonic) < _MIN_LIS_PAIRS:
        if debug_dir is not None:
            _debug_log(debug_dir,
                       f"[FULL] LIS produced {len(monotonic)} pairs (< {_MIN_LIS_PAIRS} required); returning None.")
            _plot_dotplot(debug_dir, seq_len, all_pairs_snapshot, all_mapped, monotonic,
                          chosen_offset=chosen_offset,
                          status_note=f"REJECTED: LIS n={len(monotonic)} < {_MIN_LIS_PAIRS}")
        return None

    # --- Percentile-based LTR boundary extraction ---
    end1_vals = sorted(e1 for _, e1, _, _ in monotonic)
    start2_vals = sorted(s2 for _, _, s2, _ in monotonic)
    n = len(monotonic)

    # 98th percentile for the outer boundary of the 5' LTR
    idx_98 = min(int(n * 0.98), n - 1)
    max_end1 = end1_vals[idx_98]

    # 2nd percentile for the inner boundary of the 3' LTR
    idx_02 = max(int(n * 0.02), 0)
    min_start2 = start2_vals[idx_02]

    len5p = max_end1
    len3p = seq_len - (min_start2 - 1)
    new_end1 = max_end1
    new_start2 = min_start2

    raw_end1, raw_start2 = max_end1, min_start2  # snapshot for debug

    # Symmetrize to MAX(len5p, len3p) by growing the shorter side to match.
    # The diagonal filter above has already restricted pairs to the true
    # LTR<->LTR diagonal, so growing the shorter side doesn't risk dragging
    # in unrelated sequence -- it just adds a small amount of LTR/flank
    # symmetric to the longer raw length. The extra ~30-60 bp of flank
    # (plus the +extension below) is what MAFFT needs to anchor the
    # LTR-vs-LTR alignment via a shifted global alignment; too-tight
    # extracts force MAFFT to align internal-vs-LTR end-to-end and the
    # divergence estimate inflates.
    if len5p < len3p:
        new_end1 += len3p - len5p
        if new_end1 > seq_len:
            new_end1 = seq_len
    elif len3p < len5p:
        new_start2 -= len5p - len3p
        if new_start2 < 1:
            new_start2 = 1

    sym_end1, sym_start2 = new_end1, new_start2  # snapshot post-symmetrize

    new_end1 += extension
    if new_end1 > seq_len:
        new_end1 = seq_len
    new_start2 -= extension
    if new_start2 < 1:
        new_start2 = 1

    prop_5p = new_end1
    prop_3p = seq_len - new_start2 + 1
    max_no_overlap = seq_len // 2
    final_len = min(max_no_overlap, prop_5p, prop_3p)

    new_end1 = final_len
    new_start2 = seq_len - final_len + 1

    if debug_dir is not None:
        boundaries = (
            f"# Per-stage 5'/3' boundary decisions for this LTR-RT\n"
            f"sequence_length\t{seq_len}\n"
            f"half_length\t{half_len}\n"
            f"n_lis_pairs\t{n}\n"
            f"# 98th-percentile end of 5' (max_end1) and 2nd-percentile start of 3' (min_start2):\n"
            f"raw_end1_5p\t{raw_end1}\n"
            f"raw_start2_3p\t{raw_start2}\n"
            f"raw_len5p\t{len5p}\n"
            f"raw_len3p\t{len3p}\n"
            f"# After symmetrization to MAX(len5p, len3p) (force equal\n"
            f"# lengths by growing the shorter side; safe because the\n"
            f"# diagonal filter already restricted pairs to the true LTR\n"
            f"# diagonal, so growth doesn't drag in unrelated sequence):\n"
            f"sym_end1_5p\t{sym_end1}\n"
            f"sym_start2_3p\t{sym_start2}\n"
            f"# After --extension {extension} bp tacked on each side, and capping to seq_len//2:\n"
            f"final_end1_5p\t{new_end1}\n"
            f"final_start2_3p\t{new_start2}\n"
            f"final_5p_length\t{new_end1}\n"
            f"final_3p_length\t{seq_len - new_start2 + 1}\n"
        )
        _debug_write(debug_dir, "03_ltr_boundaries.txt", boundaries)
        _plot_dotplot(debug_dir, seq_len, all_pairs_snapshot, all_mapped, monotonic,
                      chosen_offset=chosen_offset,
                      final_end1=new_end1, final_start2=new_start2)

    return (seq[0:new_end1], seq[new_start2 - 1:], new_end1, new_start2)


def _run_wfa_and_filter(aln_clean, proposed_len, mutation_rate, max_win_overdisp, verbose, debug_dir=None):
    """
    Run WFA on a trimmed alignment and post-process in Python.
    aln_clean: Path (read from file) or str (pipe via stdin).
    debug_dir: optional Path; when set, write the WFA cmd/stdout/stderr
               and add wfa.cpp's -d flag so we capture its alignment view.
    """
    # Base flags used in both modes. -k 5 -K trims unreliable termini
    # (require 5 consecutive matches to enter/exit reliable region; exclude
    # those k columns from divergence math). With debug_dir we additionally
    # pass -d so wfa.cpp emits a 3-line alignment view + boundary kmers.
    base_flags = [
        "-E", "10000", "-e", "10000", "-o", "10000", "-O", "10000",
        "-u", str(mutation_rate), "-W", "50", "-k", "5", "-K",
    ]
    if debug_dir is not None:
        base_flags.append("-d")

    if isinstance(aln_clean, str):
        cmd = [str(WFA_BIN), *base_flags, "-"]
        if verbose:
            print("Running:", " ".join(cmd), "(stdin)")
        try:
            result = subprocess.run(cmd, input=aln_clean, capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as exc:
            print(f"[SKIP] WFA failed (exit {exc.returncode}): {exc.stderr.strip()[:200]}", flush=True)
            if debug_dir is not None:
                _debug_write(debug_dir, "07_wfa_cmd.txt",
                             "stdin=clean_data (piped)\n" + " ".join(cmd))
                _debug_write(debug_dir, "08_wfa_stdout.txt",
                             f"[ERROR] WFA exited {exc.returncode}\n")
                _debug_write(debug_dir, "09_wfa_stderr.txt",
                             (exc.stderr or "") + f"\n[ERROR] exit={exc.returncode}\n")
            return
    else:
        cmd = [str(WFA_BIN), *base_flags, str(aln_clean)]
        if verbose:
            print("Running:", " ".join(cmd))
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as exc:
            print(f"[SKIP] WFA failed (exit {exc.returncode}): {exc.stderr.strip()[:200]}", flush=True)
            if debug_dir is not None:
                _debug_write(debug_dir, "07_wfa_cmd.txt", " ".join(cmd))
                _debug_write(debug_dir, "08_wfa_stdout.txt",
                             f"[ERROR] WFA exited {exc.returncode}\n")
                _debug_write(debug_dir, "09_wfa_stderr.txt",
                             (exc.stderr or "") + f"\n[ERROR] exit={exc.returncode}\n")
            return

    # Persist the divergence call for debugging even on success.
    if debug_dir is not None:
        cmd_note = (
            "# Exact WFA divergence command:\n"
            f"{' '.join(cmd)}\n"
            "# Flag meanings:\n"
            "#   -E/-e/-o/-O 10000  : large gap penalties (we expect few gaps in true LTR pairs)\n"
            f"#   -u {mutation_rate}             : mutation rate mu, used to convert divergence -> time\n"
            "#   -W 50              : 50-bp window for win_mean / win_overdisp metrics\n"
            "#   -k 5               : require 5 consecutive matches to enter/exit reliable region\n"
            "#   -K                 : exclude those k boundary-match columns from metrics\n"
            "#   -d                 : (debug-only) emit alignment view + boundary kmers to stderr\n"
        )
        _debug_write(debug_dir, "07_wfa_cmd.txt", cmd_note)
        _debug_write(debug_dir, "08_wfa_stdout.txt", result.stdout or "")
        _debug_write(debug_dir, "09_wfa_stderr.txt", result.stderr or "")

    line = result.stdout.strip()
    if not line:
        return

    fields = line.split("\t")
    # WFA output (format=full, -W 50, -k 5 -K): 19 columns
    # 0:h1  1:h2  2:cigar  3:identity  4:pairs  5:total_subs  6:transitions
    # 7:transversions  8:raw_d  9:raw_T  10:JC69_d  11:JC69_T  12:K2P_d
    # 13:K2P_T  14:left_trim  15:right_trim  16:win_n  17:win_mean  18:win_overdisp
    if len(fields) < 14:
        return

    # Equivalent of cut -f1,5-  (keep h1 and cols 5+, i.e. fields[0] + fields[4:])
    h1 = fields[0]
    kept = list(fields[4:])  # pairs onward

    # Fix negative zeros
    kept = [f.replace("-0.000000", "0.000000") for f in kept]

    # Clean header: strip _5prime* or _3prime* suffix
    h1 = re.sub(r'_(5prime|3prime).*$', '', h1)

    # kept[4] = raw_d (p-distance).  Filter: raw_d <= 0.35
    # NaN comparisons return False, so check explicitly.
    try:
        raw_d = float(kept[4])
        if raw_d != raw_d or raw_d > 0.35:  # NaN != NaN is True
            return
    except (ValueError, IndexError):
        return

    # Optional overdispersion filter on last column (always present with -W 50)
    if max_win_overdisp is not None and len(kept) > 3:
        try:
            if float(kept[-1]) > max_win_overdisp:
                return
        except ValueError:
            pass

    # Drop last 3 columns (win_n, win_mean, win_overdisp) if present
    if len(kept) > 3:
        kept = kept[:-3]

    # Insert proposed_len as second column
    return "\t".join([h1, str(proposed_len)] + kept) + "\n"


def read_existing_ids(results_path: str) -> set[str]:
    """
    Return the set of first-column IDs already present in an existing results file.
    Lines may be blank; tolerate non-numeric fields elsewhere.
    """
    ids = set()
    try:
        with open(results_path) as fh:
            for line in fh:
                s = line.strip()
                if not s:
                    continue
                # first token is the LTR-RT name after awk header cleanup in the pipeline
                ids.add(s.split()[0])
    except FileNotFoundError:
        pass
    return ids

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


def try_fast_path(dir_path: Path, header_token: str, full_header: str, seq_len: int, seq: str):
    """
    Fast-path when domains TSV provides LTR length:
      - Safety checks
      - Extract first/last min(ltr_len + extension, seq_len//2) bp
      - Run MAFFT -> trimal -> WFA, return result line
      - Optionally also build an IUPAC consensus from the trimmed LTRs.

    Returns (attempted, (result_line, consensus_record)):
      (False, (None, None))           fast-path doesn't apply (caller falls back)
      (True,  (None, None))           fast-path ran but result was filtered/skipped
      (True,  (result, None))         result produced, consensus disabled or skipped
      (True,  (result, consensus))    result and consensus FASTA record produced
    """
    global ARGS, DOMAINS
    if not DOMAINS:
        return (False, (None, None))

    # Duplicate name? -> uncertain mapping -> full pipeline unless user overrides.
    if (dir_path / "DUPLICATE_NAME").exists() and not ARGS.assume_dup_same_ltr:
        if ARGS.verbose:
            print(f"[FAST-PATH SKIP] Duplicate header name detected for '{header_token}'. "
                  f"Use --assume-duplicate-same-ltr to override (risky).")
        return (False, (None, None))
    # If ARGS.assume_dup_same_ltr is True, proceed assuming all duplicates with this
    # header token share the same LTR length from DOMAINS.

    if header_token not in DOMAINS:
        return (False, (None, None))

    ltr_len = DOMAINS[header_token]

    # Safety check: length_of_element > 2*ltr_len + 50 (unchanged)
    if not (seq_len > (2 * ltr_len + 50)):
        if ARGS.verbose:
            print(f"[FAST-PATH SKIP] Safety check failed for '{header_token}': "
                  f"seq_len={seq_len} <= 2*{ltr_len}+50")
        return (False, (None, None))

    debug_dir = _ensure_debug_dir(dir_path)
    _debug_log(debug_dir,
               f"=== FAST PATH ===\n"
               f"header_token: {header_token}\n"
               f"seq_len: {seq_len}\n"
               f"domains_ltr_len: {ltr_len}\n"
               f"extension: {ARGS.extension}\n")

    # Determine extraction length: min(ltr_len + extension, seq_len//2)
    ext_len = ltr_len + ARGS.extension
    cap = seq_len // 2
    if ext_len > cap:
        ext_len = cap

    # Perform direct extraction
    five_p = seq[:ext_len]
    three_p = seq[-ext_len:]

    fasta_str = (f">{header_token}_5prime_LTR_len{ltr_len}_ext{ext_len}\n"
                 f"{five_p}\n"
                 f">{header_token}_3prime_LTR_len{ltr_len}_ext{ext_len}\n"
                 f"{three_p}\n")

    _debug_write(debug_dir, "01_fastpath_extracted.fa", fasta_str)
    _debug_log(debug_dir,
               f"Cut first/last {ext_len} bp (capped at seq_len//2={cap}).\n"
               f"5' coords: 1-{ext_len}\n"
               f"3' coords: {seq_len - ext_len + 1}-{seq_len}\n"
               f"Note: in the FULL pipeline the boundaries come from kmer LIS,\n"
               f"      which can pick very different coordinates. If --domains\n"
               f"      and no-domains give very different divergence, the most\n"
               f"      common cause is exactly this — the cuts disagree, and so\n"
               f"      do all downstream alignments.")

    if ARGS.verbose:
        print(f"[FAST-PATH] {dir_path.name}: extracted 5′/3′ ({ext_len} bp each) using domains TSV (LTR={ltr_len}, ext={ARGS.extension}).")

    # In-memory pipeline: MAFFT -> trimal -> quality check
    try:
        aln_data = _run_align_mem(fasta_str, ARGS.verbose)
    except subprocess.CalledProcessError:
        print(f"[SKIP] mafft failed for {dir_path.name}")
        _debug_log(debug_dir, "[FAST] alignment step failed.")
        return (True, (None, None))

    aligner_name = "WFA (--wfa-align)" if getattr(ARGS, "wfa_align", False) else "MAFFT --auto"
    _debug_write(debug_dir, "05_alignment.fa", aln_data)
    _debug_log(debug_dir, f"Pairwise alignment via {aligner_name} written to 05_alignment.fa.")

    raw_bp = _sum_ungapped_bp_str(aln_data)

    try:
        clean_data = _run_trimal_mem(aln_data, ARGS.verbose)
    except subprocess.CalledProcessError:
        print(f"[SKIP] trimal failed for {dir_path.name}")
        _debug_log(debug_dir, "[FAST] trimal step failed.")
        return (True, (None, None))

    _debug_write(debug_dir, "06_trimal.fa", clean_data)
    _debug_log(debug_dir, "trimal -automated1 -keepheader written to 06_trimal.fa.")

    clean_bp = _sum_ungapped_bp_str(clean_data)

    threshold = ARGS.min_retained_fraction * raw_bp
    if (clean_bp + ARGS.extension) <= threshold:
        log_msg(
            f"[SKIP] Trimmed alignment too short for {dir_path.name}: "
            f"clean_bp({clean_bp}) + ext({ARGS.extension}) <= "
            f"{ARGS.min_retained_fraction:.3f} * raw_bp({raw_bp:.0f})"
        )
        _debug_log(debug_dir,
                   f"[FAST] SKIP: clean_bp({clean_bp}) + ext({ARGS.extension}) "
                   f"<= {ARGS.min_retained_fraction:.3f} * raw_bp({raw_bp:.0f}).")
        return (True, (None, None))  # fast-path executed, but result skipped

    # Fast-path: use DOMAINS_TSV value directly, do NOT subtract extension
    reported_len = ltr_len

    result = _run_wfa_and_filter(
        clean_data, reported_len,
        ARGS.mutation_rate, ARGS.max_win_overdisp, ARGS.verbose,
        debug_dir=debug_dir,
    )
    if debug_dir is not None:
        _debug_log(debug_dir,
                   ("[FAST] WFA divergence produced a result (see 08_wfa_stdout.txt)."
                    if result is not None
                    else "[FAST] WFA divergence yielded nothing (filtered or error)."))

    consensus_record = None
    if result is not None and getattr(ARGS, 'ltr_cluster', False) and not getattr(ARGS, 'wfa_align', False):
        consensus_record = _build_consensus_from_trimmed(clean_data, full_header, ARGS.verbose, debug_dir=debug_dir)
        if debug_dir is not None and consensus_record:
            _debug_write(debug_dir, "11_consensus.fa", consensus_record)

    return (True, (result, consensus_record))


def process_dir(dir_path):
    """
    Process one sequence directory:
      - If domains TSV provides a safe fast-path, extract LTRs directly with extension,
        then run MAFFT -> trimal -> WFA.
      - Otherwise run the original heavy pipeline (k-mers -> map -> filter -> extract LTRs -> MAFFT -> trimal -> WFA).

    When ARGS.purge_subdirs is True, the per-sequence temp directory is removed
    at the end of processing (even if skipped or errors occur).

    Always returns a tuple (divergence_line, consensus_record); either element may be None.
    """
    dir_path = Path(dir_path)

    try:
        args = ARGS
        verbose = args.verbose
        dist = args.dist
        kmin = args.kmin
        kmax = args.kmax
        std_factor = args.std_factor
        extension = args.extension

        if ARGS.verbose:
            print(f"Processing {dir_path}")
        seq_fa = dir_path / "seq.fa"

        # Read header and sequence once (avoids repeated file reads)
        full_header, seq = _read_fasta_entry(seq_fa)
        header_token = full_header.split()[0] if full_header else dir_path.name
        seq_len = len(seq)

        if seq_len < MIN_SEQ_BP:
            log_msg(f"[SKIP] {dir_path.name} length {seq_len} < {MIN_SEQ_BP} bp; skipping as likely erroneous.")
            try:
                (dir_path / "SKIPPED.too_short").write_text(f"seq_len={seq_len}\nmin={MIN_SEQ_BP}\n")
            except Exception:
                pass
            d = _ensure_debug_dir(dir_path)
            _debug_log(d, f"=== SKIPPED.too_short ===\nseq_len={seq_len} < MIN_SEQ_BP={MIN_SEQ_BP}.")
            return (None, None)

        # ===== FAST-PATH TRY =====
        try:
            attempted, result_tuple = try_fast_path(dir_path, header_token, full_header, seq_len, seq)
            if attempted:
                return result_tuple
        except Exception as e:
            print(f"[FAST-PATH ERROR] {header_token}: {e}. Falling back to full pipeline.")

        # ===== HEAVY PIPELINE (in-process kmer discovery) =====

        debug_dir = _ensure_debug_dir(dir_path)
        _debug_log(debug_dir,
                   f"=== FULL PIPELINE ===\n"
                   f"header_token: {header_token}\n"
                   f"seq_len: {seq_len}\n"
                   f"kmin: {kmin}  kmax: {kmax}  dist: {dist}  extension: {extension}\n")

        ltr_result = _discover_and_extract_ltrs(seq, kmin, kmax, dist, std_factor, extension, debug_dir=debug_dir)
        if ltr_result is None:
            log_msg(f"[SKIP] No valid kmer pairs for {dir_path.name}.")
            try:
                (dir_path / "SKIPPED.no_filtered_pairs").write_text("")
            except Exception:
                pass
            _debug_log(debug_dir, "[FULL] No filtered kmer pairs; sequence skipped.")
            return (None, None)

        ltr_5p, ltr_3p, end5p, start3p = ltr_result
        fasta_str = (f">{full_header}\t5p:1-{end5p}\n{ltr_5p}\n"
                     f">{full_header}\t3p:{start3p}-{seq_len}\n{ltr_3p}\n")

        _debug_write(debug_dir, "04_extracted_ltrs.fa", fasta_str)
        _debug_log(debug_dir,
                   f"Extracted 5' LTR (1-{end5p}, {end5p} bp) and 3' LTR "
                   f"({start3p}-{seq_len}, {seq_len - start3p + 1} bp).")

        # In-memory pipeline: MAFFT -> trimal -> quality check -> WFA
        try:
            aln_data = _run_align_mem(fasta_str, verbose)
        except subprocess.CalledProcessError:
            print(f"[SKIP] mafft failed for {dir_path.name}")
            _debug_log(debug_dir, "[FULL] alignment step failed.")
            return (None, None)

        aligner_name = "WFA (--wfa-align)" if getattr(ARGS, "wfa_align", False) else "MAFFT --auto"
        _debug_write(debug_dir, "05_alignment.fa", aln_data)
        _debug_log(debug_dir, f"Pairwise alignment via {aligner_name} written to 05_alignment.fa.")

        proposed_len, raw_bp = _read_aln_stats_str(aln_data)

        try:
            clean_data = _run_trimal_mem(aln_data, verbose)
        except subprocess.CalledProcessError:
            print(f"[SKIP] trimal failed for {dir_path.name}")
            _debug_log(debug_dir, "[FULL] trimal step failed.")
            return (None, None)

        _debug_write(debug_dir, "06_trimal.fa", clean_data)
        _debug_log(debug_dir, "trimal -automated1 -keepheader written to 06_trimal.fa.")

        clean_bp = _sum_ungapped_bp_str(clean_data)

        # Require: sum_bp(clean) + extension > (min_retained_fraction * sum_bp(raw))
        threshold = args.min_retained_fraction * raw_bp

        if (clean_bp + extension) <= threshold:
            log_msg(
                f"[SKIP] Trimmed alignment too short for {dir_path.name}: "
                f"clean_bp({clean_bp}) + ext({extension}) <= "
                f"{args.min_retained_fraction:.3f} * raw_bp({raw_bp:.0f})"
            )
            _debug_log(debug_dir,
                       f"[FULL] SKIP: clean_bp({clean_bp}) + ext({extension}) "
                       f"<= {args.min_retained_fraction:.3f} * raw_bp({raw_bp:.0f}).")
            return (None, None)

        # Adjust proposed LTR len by extension
        if proposed_len is None:
            proposed_adj = 0
            if ARGS.verbose:
                print(f"[WARN] Could not parse proposed LTR length for {dir_path.name}; inserting 0.")
        else:
            proposed_adj = max(proposed_len - ARGS.extension, 0)

        # Final WFA (in-memory via stdin)
        result = _run_wfa_and_filter(
            clean_data, proposed_adj,
            args.mutation_rate, args.max_win_overdisp, verbose,
            debug_dir=debug_dir,
        )
        if debug_dir is not None:
            _debug_log(debug_dir,
                       ("[FULL] WFA divergence produced a result (see 08_wfa_stdout.txt)."
                        if result is not None
                        else "[FULL] WFA divergence yielded nothing (filtered or error)."))

        consensus_record = None
        if result is not None and getattr(ARGS, 'ltr_cluster', False) and not getattr(ARGS, 'wfa_align', False):
            consensus_record = _build_consensus_from_trimmed(clean_data, full_header, verbose, debug_dir=debug_dir)
            if debug_dir is not None and consensus_record:
                _debug_write(debug_dir, "11_consensus.fa", consensus_record)

        return (result, consensus_record)

    finally:
        # Always try to remove the subdir if purge_subdirs is enabled
        if getattr(ARGS, "purge_subdirs", False):
            try:
                shutil.rmtree(dir_path)
            except FileNotFoundError:
                # already removed or never created; that's fine
                pass
            except Exception as e:
                log_msg(f"[WARN] Failed to remove temp subdir {dir_path}: {e}")

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


def prefix_before_dot(path_like: str) -> str:
    """Return filename stem up to the first '.' (e.g., 'ZMHA_mod.pass.list.fa' -> 'ZMHA_mod')."""
    name = Path(path_like).name
    return name.split('.', 1)[0]


def domains_prefix(path_like: str) -> str:
    """Return domains file prefix for matching input FASTA prefixes.
    We take the part before the first '.' and strip a trailing '_domains' if present.
    e.g., 'ZMHA_mod_domains.tsv' -> 'ZMHA_mod'
    """
    stem = Path(path_like).name.split('.', 1)[0]
    return re.sub(r'_domains$', '', stem)


def _progress_update(completed, total, start, last_update, update_interval):
    """Print a progress update to stderr if enough time has elapsed. Returns new last_update."""
    now = time.monotonic()
    if (now - last_update) >= update_interval or completed == total:
        elapsed = now - start
        pct = (completed / total) * 100 if total else 100.0
        rate = (completed / elapsed) if elapsed > 0 else 0.0
        remaining = ((total - completed) / rate) if rate > 0 else float("inf")
        eta_text = _format_hms(remaining)
        msg = (f"\rProcessing {completed}/{total}. "
               f"{pct:.2f}% complete. Estimated time remaining {eta_text}")
        print(msg, end="", file=sys.stderr, flush=True)
        return now
    return last_update


def _consensus_path_for_outfile(outfile: str) -> str:
    """Derive consensus FASTA path from the divergence outfile path.
    If outfile ends with '.results', replace it with '.consensus.fa'; otherwise append '.consensus.fa'.
    """
    p = Path(outfile)
    if p.suffix == ".results":
        return str(p.with_suffix(".consensus.fa"))
    return f"{outfile}.consensus.fa"


def _run_mmseqs_cluster(consensus_fasta, threads: int, verbose: bool):
    """
    Cluster the IUPAC consensus LTR FASTA with mmseqs easy-cluster, retaining
    only the *_cluster.tsv. Auxiliary outputs (*_all_seqs.fasta,
    *_rep_seq.fasta, and the mmseqs tmp dir) are deleted.

    Current parameters cluster to lineage level (ie, lots of lumping). For clades and subclade level, just adjust --min-seq-id up (0.75 for sub-clade; 0.85 for recent bursts).

    Parameters were chosen via a large grid search benchmarked on Arabidopsis
    LTR annotations, with the goal of jointly minimizing:
      - singletons (single-copy clusters), and
      - family mixing (e.g., a Tork co-clustering with an Ogre is a bad sign;
        Tork co-clustering with Unknown is acceptable).

    Findings from that grid search:
      (1) mmseqs easy-cluster outperformed cd-hit-est and a custom wavefront
          alignment + clustering pipeline. Wavefront also doesn't scale: it is
          not feasible for species with more than ~2k LTR-RTs.
      (2) The IUPAC consensus LTR FASTA (this pipeline's output) outperformed
          both the full-length LTR-RT FASTA and a 5'-LTR-only FASTA as input
          to mmseqs.
    mmseqs is fast, so this step adds little to total runtime.

    Returns the Path to <prefix>_cluster.tsv on success, or None if the
    consensus FASTA is missing/empty, mmseqs is unavailable, or mmseqs fails.
    """
    consensus_fasta = Path(consensus_fasta)
    if not consensus_fasta.exists() or consensus_fasta.stat().st_size == 0:
        print(f"[CLUSTER] Consensus FASTA {consensus_fasta} is missing or empty; "
              f"skipping mmseqs.")
        return None

    cluster_prefix = consensus_fasta.with_suffix("")  # strips final '.fa'
    tmp_dir = Path(tempfile.mkdtemp(
        prefix=f"{consensus_fasta.stem}_mmseqs_",
        dir=str(consensus_fasta.parent),
    ))

    cmd = [
        "mmseqs", "easy-cluster",
        str(consensus_fasta), str(cluster_prefix), str(tmp_dir),
        "--min-seq-id", "0.6",
        "-c", "0.6",
        "--cov-mode", "1",
        "--cluster-mode", "1",
        "--mask", "0",
        "-s", "7.5",
        "--threads", str(threads),
    ]
    if verbose:
        print("Running:", " ".join(cmd))

    try:
        subprocess.run(
            cmd,
            stdout=(None if verbose else subprocess.DEVNULL),
            stderr=(None if verbose else subprocess.DEVNULL),
            check=True,
        )
    except FileNotFoundError:
        print("[CLUSTER] mmseqs not found in PATH; cannot cluster consensus LTRs.",
              file=sys.stderr)
        shutil.rmtree(tmp_dir, ignore_errors=True)
        return None
    except subprocess.CalledProcessError as exc:
        print(f"[CLUSTER] mmseqs easy-cluster failed (exit {exc.returncode}).",
              file=sys.stderr)
        shutil.rmtree(tmp_dir, ignore_errors=True)
        return None

    shutil.rmtree(tmp_dir, ignore_errors=True)
    for ext in ("_all_seqs.fasta", "_rep_seq.fasta"):
        aux = Path(f"{cluster_prefix}{ext}")
        try:
            aux.unlink()
        except FileNotFoundError:
            pass

    cluster_tsv = Path(f"{cluster_prefix}_cluster.tsv")
    if not cluster_tsv.exists():
        print(f"[CLUSTER] Expected cluster TSV not found: {cluster_tsv}",
              file=sys.stderr)
        return None
    return cluster_tsv


def process_one_input(args_base: argparse.Namespace, in_fasta: str, per_prefix_domains: dict, multi_mode: bool):
    """Process a single input FASTA using (possibly) matched domains mapping in multi-input mode."""
    global ARGS, DOMAINS, LOGFILE

    pref = prefix_before_dot(in_fasta)
    # Decide temp/out per mode
    temp_dir = f"{pref}_temp" if multi_mode else args_base.temp_dir
    outfile = f"{pref}.LTRs.alns.results" if multi_mode else args_base.outfile

    # Prepare workspace (per input)
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    os.makedirs(temp_dir, exist_ok=True)

    if args_base.reuse_existing and os.path.exists(outfile):
        # do not truncate; we will append only missing items
        pass
    else:
        open(outfile, "w").close()

    # Optional consensus FASTA output
    consensus_outfile = None
    if getattr(args_base, "ltr_cluster", False):
        consensus_outfile = _consensus_path_for_outfile(outfile)
        if args_base.reuse_existing and os.path.exists(consensus_outfile):
            # mirror reuse semantics: keep existing records, append new ones
            pass
        else:
            open(consensus_outfile, "w").close()

    # Pick matching domains mapping (None if absent)
    DOMAINS = per_prefix_domains.get(pref)
    if args_base.domains_tsvs and multi_mode and DOMAINS is None:
        print(f"[DOMAINS] No matching domains TSV for prefix '{pref}'. Proceeding without fast-path.")

    # Create a fresh ARGS namespace for this input
    ARGS = deepcopy(args_base)
    ARGS.temp_dir = temp_dir
    ARGS.outfile = outfile
    ARGS.consensus_outfile = consensus_outfile
    ARGS.input_fasta = in_fasta

    log_path = outfile + ".log"
    # create/clear the log file once in the parent
    with open(log_path, "w") as _lf:
        _lf.write("")  # truncate

    use_batching = args_base.purge_subdirs

    if use_batching:
        _process_batched(ARGS, DOMAINS, log_path, in_fasta, pref)
    else:
        _process_unbatched(ARGS, DOMAINS, log_path, in_fasta, pref)

    write_summary(ARGS.outfile)

    # Cleanup
    if ARGS.keep_temp:
        print(f"Keeping temporary directory {ARGS.temp_dir}")
    else:
        print(f"Removing temporary directory {ARGS.temp_dir}")
        shutil.rmtree(ARGS.temp_dir)

    print(f"All sequences processed for {pref}. Output in {ARGS.outfile}")
    if consensus_outfile:
        print(f"LTR consensus FASTA in {consensus_outfile}")
        cluster_tsv = _run_mmseqs_cluster(
            consensus_outfile, args_base.threads, args_base.verbose
        )
        if cluster_tsv:
            print(f"LTR cluster TSV in {cluster_tsv}")


def _process_unbatched(args, domains, log_path, in_fasta, pref):
    """Original flow: split all sequences upfront, then process all at once."""
    # Split FASTA into separate seq.fa files
    split_fasta(args.input_fasta, args.temp_dir)

    # Map sanitized dir -> original header token
    header_map = _read_header_map(args.temp_dir)

    # Discover sequence directories
    dirs = [p for p in Path(args.temp_dir).iterdir() if p.is_dir()]

    # Filter to only missing IDs if reusing
    if args.reuse_existing:
        already = read_existing_ids(args.outfile)
        def needs_run(p: Path) -> bool:
            token = header_map.get(p.name)
            if token is None:
                try:
                    with open(p / "seq.fa") as fh:
                        first = next((ln for ln in fh if ln.startswith(">")), None)
                    token = first[1:].strip().split()[0] if first else p.name
                except Exception:
                    token = p.name
            return token not in already

        dirs = [str(p) for p in dirs if needs_run(p)]
        if args.verbose:
            print(f"[REUSE] {len(already)} already in {args.outfile}; {len(dirs)} remaining to process.")
    else:
        dirs = [str(p) for p in dirs]

    total = len(dirs)
    start = time.monotonic()
    last_update = 0.0
    update_interval = 0.25
    completed = 0
    print("", file=sys.stderr)

    consensus_fh = open(args.consensus_outfile, "a") if args.consensus_outfile else None
    try:
        with open(args.outfile, "a") as outfh, Pool(
            processes=args.threads,
            initializer=_pool_init,
            initargs=(args, domains, log_path),
        ) as pool:
            for result_tuple in pool.imap_unordered(process_dir, dirs, chunksize=1):
                result_line, consensus_record = result_tuple
                if result_line:
                    outfh.write(result_line)
                    outfh.flush()
                if consensus_record and consensus_fh is not None:
                    consensus_fh.write(consensus_record)
                    consensus_fh.flush()
                completed += 1
                last_update = _progress_update(completed, total, start, last_update, update_interval)
    finally:
        if consensus_fh is not None:
            consensus_fh.close()

    print("", file=sys.stderr)


def _process_batched(args, domains, log_path, in_fasta, pref):
    """
    Batched flow: index the FASTA without writing seq.fa files,
    then materialize and process one batch at a time.
    """
    batch_size = args.batch_size

    # Step 1: Index — scan only, no files written (except header_map.tsv)
    entries = index_fasta(args.input_fasta, args.temp_dir)

    # Step 2: Filter for reuse
    if args.reuse_existing:
        already = read_existing_ids(args.outfile)
        entries = [e for e in entries if e.header_token not in already]
        if args.verbose:
            print(f"[REUSE] {len(already)} already in {args.outfile}; {len(entries)} remaining to process.")

    total = len(entries)
    start = time.monotonic()
    last_update = 0.0
    update_interval = 0.25
    completed = 0
    print("", file=sys.stderr)

    consensus_fh = open(args.consensus_outfile, "a") if args.consensus_outfile else None
    try:
        # Create pool once, reuse across all batches
        with open(args.outfile, "a") as outfh, Pool(
            processes=args.threads,
            initializer=_pool_init,
            initargs=(args, domains, log_path),
        ) as pool:
            # Process in batches
            for batch_start in range(0, total, batch_size):
                batch = entries[batch_start:batch_start + batch_size]

                # Step 3a: Write only this batch's seq.fa files
                write_batch(batch, args.input_fasta, args.temp_dir)

                # Step 3b: Process this batch
                batch_dirs = [str(Path(args.temp_dir) / e.safe_name) for e in batch]
                for result_tuple in pool.imap_unordered(process_dir, batch_dirs, chunksize=1):
                    result_line, consensus_record = result_tuple
                    if result_line:
                        outfh.write(result_line)
                        outfh.flush()
                    if consensus_record and consensus_fh is not None:
                        consensus_fh.write(consensus_record)
                        consensus_fh.flush()
                    completed += 1
                    last_update = _progress_update(completed, total, start, last_update, update_interval)
                # purge_subdirs is active, so process_dir already removed each subdir
    finally:
        if consensus_fh is not None:
            consensus_fh.close()

    print("", file=sys.stderr)


def _read_header_map(temp_dir):
    """Read header_map.tsv: sanitized_name -> original_header_token."""
    header_map = {}
    try:
        with open(Path(temp_dir) / "header_map.tsv") as hm:
            for line in hm:
                s = line.strip()
                if not s:
                    continue
                safe, token = s.split("\t", 1)
                header_map[safe] = token
    except FileNotFoundError:
        pass
    return header_map


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Process multi-seq LTR-RT FASTA(s) to extract and align LTRs. "
            "Optionally fast-path known LTR lengths via a domains TSV.\n\n"
            "Single-input: -t/--temp-dir and -o/--outfile apply as usual.\n"
            "Multi-input: -t/--temp-dir and -o/--outfile are ignored; per-file temp/out are inferred from the input filename prefix (text before first '.')."
        )
    )
    parser.add_argument(
        "-k", action="store_true", dest="keep_temp",
        help="Keep temp directory after processing."
    )
    parser.add_argument(
        "--purge-subdirs", action="store_true", dest="purge_subdirs",
        help=(
            "After each sequence directory finishes processing (including skips/errors), "
            "delete its temp subdirectory to keep file counts low. "
            "Mutually exclusive with -k/--keep_temp."
        )
    )
    parser.add_argument(
        "--batch-size", type=int, default=100, dest="batch_size",
        help=(
            "When --purge-subdirs is active, only materialize this many seq.fa files "
            "at a time to limit temp file count. Default: 100. "
            "Ignored when --purge-subdirs is not set."
        )
    )
    parser.add_argument(
        "-v", action="store_true", dest="verbose",
        help="Verbose mode; print each command before executing."
    )
    parser.add_argument(
        "--debug", action="store_true", dest="debug",
        help=(
            "Debug mode: implies -v (verbose) and -k (keep temp), and writes "
            "many per-sequence intermediates under <temp_dir>/<seq>/debug/ "
            "(kmer pairs, LIS pairs, boundary maths, extracted LTR FASTAs, "
            "MAFFT/WFA alignments, trimal output, full WFA divergence command "
            "+ stdout + stderr with the wfa.cpp '-d' flag enabled, IUPAC "
            "consensus). Adds a startup banner explaining the pipeline. "
            "Speed is not the goal; this is for diagnosing why two runs disagree "
            "(domains vs no-domains, --wfa-align vs MAFFT, -k boundary trimming, "
            "unexpected skips). Incompatible with --purge-subdirs."
        )
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
        help="Temporary directory name (single input only; ignored with multiple inputs)."
    )
    parser.add_argument(
        "-o", default="./LTRs.alns.results", dest="outfile",
        help="Output filename (single input only; ignored with multiple inputs)."
    )
    parser.add_argument(
        "-p", type=int, default=20, dest="threads",
        help="Number of parallel threads per input (default: 20)."
    )
    # in argparse section
    parser.add_argument(
        "--reuse-existing", action="store_true", dest="reuse_existing",
        help="If the per-input results file already exists, only process LTR-RTs missing from it and append."
    )
    parser.add_argument(
        "-D", "--domains", dest="domains_tsvs", nargs="*", default=None,
        help=(
            "Optional domains TSV file(s) (format: name\\tLTR_len). "
            "With a single input, provide one TSV. With multiple inputs, you may supply multiple TSVs; "
            "each TSV is matched to an input by prefix. Matching uses the TSV filename up to the first '.'; "
            "if it ends with '_domains', that suffix is ignored for matching."
        )
    )
    parser.add_argument(
        "--max-win-overdisp", type=float, default=6, dest="max_win_overdisp",
        help=(
            "Maximum allowed window overdispersion (win_overdisp). "
            "Lower values are more specific (aggressive filtering), "
            "higher values are more sensitive (retains more outputs). "
            "To disable this filter entirely, pass a very large value "
            "(e.g. 'inf' or 1e9). "
            "Default: 6."
        )
    )

    parser.add_argument(
        "--min-retained-fraction", type=float, default=0.6, dest="min_retained_fraction",
        help=(
            "Minimum fraction of ungapped columns retained after trimming required to proceed. "
            "Higher values are more specific (aggressive filtering), "
            "lower values are more sensitive (retains more outputs). "
            "To disable this filter entirely, pass 0. "
            "Default: 0.6."
        )
    )
    parser.add_argument(
        "--assume-duplicate-same-ltr", action="store_true", dest="assume_dup_same_ltr",
        help=(
            "Override duplicate-header safety fallback in fast path. "
            "Assumes all duplicate header tokens share the same LTR length from the domains TSV. Use with caution."
        )
    )
    parser.add_argument(
        "--no-plot", action="store_true", dest="no_plot",
        help="Disable plotting of results into kmer2ltr_density.pdf."
    )
    parser.add_argument(
        "--wfa-align", action="store_true", dest="wfa_align",
        help=(
            "Use WFA instead of mafft for pairwise LTR alignment. "
            "Much faster (~30-50x) but uses a different alignment algorithm, "
            "so divergence estimates may differ slightly. "
            "Default: use mafft."
        )
    )
    parser.add_argument(
        "--ltr-cluster", action="store_true", dest="ltr_cluster",
        help=(
            "Build a FASTA of IUPAC consensus LTRs (one per input LTR-RT) AND cluster "
            "them with mmseqs easy-cluster. "
            "Consensus pipeline: MAFFT -> trimal -> WFA realignment of trimmed LTRs -> "
            "IUPAC consensus; headers match the input LTR-RT FASTA exactly. "
            "Clustering parameters were chosen via a large grid search benchmarked on "
            "Arabidopsis LTR annotations to jointly minimize singletons and "
            "cross-family clusters. "
            "Single-input outputs: <outfile>.consensus.fa and "
            "<outfile>.consensus_cluster.tsv (with .results replaced). "
            "Multi-input outputs: <prefix>.LTRs.alns.consensus.fa and "
            "<prefix>.LTRs.alns.consensus_cluster.tsv per input. "
            "Requires mmseqs in PATH. Not compatible with --wfa-align."
        )
    )
    parser.add_argument(
        "-i", "--input-fastas",
        nargs="+", required=True,
        help="Path(s) to multi-sequence LTR-RT FASTA file(s)."
    )
#    parser.add_argument(
#        "input_fastas", nargs="+",
#        help="Path(s) to multi-sequence LTR-RT FASTA file(s)."
#    )

    args = parser.parse_args()

    # --debug overrides perf-oriented flags so we keep everything for inspection.
    if args.debug:
        if args.purge_subdirs:
            parser.error("--debug is incompatible with --purge-subdirs "
                         "(debug needs temp subdirs retained for inspection).")
        args.verbose = True
        args.keep_temp = True
        # Wire ARGS up here so the banner / any pre-pool prints know debug is on.
        ARGS = args
        print_debug_banner()

    # Enforce mutual exclusivity: keep_temp vs purge_subdirs
    if args.keep_temp and args.purge_subdirs:
        parser.error("Options -k/--keep_temp and --purge-subdirs are mutually exclusive.")

    # Enforce mutual exclusivity: ltr_cluster vs wfa_align
    if args.ltr_cluster and args.wfa_align:
        parser.error("--ltr-cluster and --wfa-align are mutually exclusive: "
                     "the consensus + clustering pipeline requires the initial alignment to be MAFFT.")

    # --ltr-cluster shells out to mmseqs; fail fast if it isn't installed.
    if args.ltr_cluster and shutil.which("mmseqs") is None:
        parser.error("--ltr-cluster requires mmseqs in PATH; "
                     "install (e.g. 'mamba install -c bioconda mmseqs2') and retry.")

    # Multi-input mode check and warnings for -t/-o
    multi_mode = len(args.input_fastas) > 1
    if multi_mode and (args.temp_dir != "./temp" or args.outfile != "./LTRs.alns.results"):
        print("[WARN] Multiple inputs provided; -t/--temp-dir and -o/--outfile are ignored.", file=sys.stderr)

    # Build per-prefix domains mapping
    per_prefix_domains: dict[str, dict[str, int] | None] = {}
    if args.domains_tsvs:
        for dpath in args.domains_tsvs:
            pref = domains_prefix(dpath)
            if pref in per_prefix_domains:
                print(f"[DOMAINS] Warning: duplicate domains TSV for prefix '{pref}'. Overwriting previous mapping.")
            mapping = read_domains(dpath)
            if mapping:
                per_prefix_domains[pref] = mapping

    # Process inputs in the order provided
    for in_fa in args.input_fastas:
        process_one_input(args, in_fa, per_prefix_domains, multi_mode)

    # === Plotting step (optional) ===
    # Collect all output results files
    if multi_mode:
        result_files = [f"{prefix_before_dot(fa)}.LTRs.alns.results" for fa in args.input_fastas]
    else:
        result_files = [args.outfile]

    # Filter to existing files (in case some were skipped entirely)
    result_files = [rf for rf in result_files if os.path.exists(rf)]
    if not args.no_plot:
        if not result_files:
            print("[PLOT] No results found to plot; skipping.", file=sys.stderr)
        else:
            plot_cmd = [
                sys.executable, str(SCRIPT_DIR / "plot_density.py"),
                "-in", *result_files,
                "-model", "K2P",
                "-miu", str(args.mutation_rate),
                "-out", "kmer2ltr_density.pdf",
                "--xmax", "0.2",
                "--bins", "auto",
                "--label-peaks"
            ]
            if args.verbose:
                print("Running:", " ".join(map(str, plot_cmd)))
            subprocess.run(plot_cmd, check=True)
            print("Density plot saved as kmer2ltr_density.pdf")
    else:
        print("Skipping plotting step (--no-plot).")
