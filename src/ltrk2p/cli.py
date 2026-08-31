"""Command-line interface."""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

from . import align
from .runner import count_data_lines, run


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="ltrk2p",
        description="Classify LTR-RT boundaries and report K2P divergence between the LTRs.",
    )
    p.add_argument("input", help="multi-FASTA of putative intact LTR-RTs (.fa or .fa.gz)")
    p.add_argument("-o", "--output", default="-", help="output TSV (default: stdout)")
    p.add_argument("--cs", action="store_true",
                   help="emit a minimap2 cs string instead of an extended CIGAR")
    p.add_argument("-t", "--threads", type=int, default=1, help="worker processes")
    p.add_argument("--resume", action="store_true",
                   help="skip records already present in the output and append")
    p.add_argument("-v", "--verbose", action="store_true", help="per-step progress")
    adv = p.add_argument_group("advanced (benchmark-calibrated defaults)")
    adv.add_argument("--flank-bits", type=float, default=None,
                     help="pin the evidence in bits required to call a flank; "
                          "default is a divergence-aware schedule keyed on the "
                          "element's own estimated divergence")
    adv.add_argument("--min-bitscore", type=float, default=None,
                     help="minimum alignment bit score to report a pair")
    adv.add_argument("--max-window", type=int, default=None,
                     help="cap on the prefix/suffix search window in bp")
    return p


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    inp = Path(args.input)
    if not inp.exists():
        print(f"ltrk2p: error: input file not found: {inp}", file=sys.stderr)
        return 2

    if args.threads < 1:
        print(f"ltrk2p: error: --threads must be >= 1, got {args.threads}", file=sys.stderr)
        return 2

    # Thread tuning values explicitly. Assigning align.T_BITS / align.W0 would
    # be a silent no-op: both are already bound as default arguments at def time.
    classify_kw = {}
    if args.flank_bits is not None:
        classify_kw["t_bits"] = args.flank_bits
    if args.max_window is not None:
        classify_kw["w0"] = args.max_window
    if args.min_bitscore is not None:
        classify_kw["min_bitscore"] = args.min_bitscore

    skip = 0
    if args.resume:
        if args.output == "-":
            print("ltrk2p: error: --resume requires -o/--output", file=sys.stderr)
            return 2
        skip = count_data_lines(args.output)
        if args.verbose:
            print(f"resuming: {skip} records already done", file=sys.stderr)

    t0 = time.time()
    print(f"ltrk2p: reading {inp}", file=sys.stderr)
    if args.output == "-":
        n = run(inp, sys.stdout, args.threads, args.cs, skip, args.verbose,
                resuming=args.resume, **classify_kw)
    else:
        # append whenever resuming -- even at skip == 0, where the header was already
        # flushed but no record finished. `resuming` is passed explicitly because
        # inferring it from `skip` writes a second header into the middle of the data.
        mode = "a" if args.resume else "w"
        with open(args.output, mode) as fh:
            n = run(inp, fh, args.threads, args.cs, skip, args.verbose,
                    resuming=args.resume, **classify_kw)
    print(f"ltrk2p: {n} records in {time.time() - t0:.1f}s", file=sys.stderr)
    return 0
