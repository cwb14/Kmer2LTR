"""Command-line interface."""
from __future__ import annotations

import argparse
import importlib.util
import sys
import time
from pathlib import Path

from . import align, cluster
from .extras import PERFECT_MODES, ExtraSpec, ExtraWriter, stream_paths
from .genome import Options
from .runner import run, scan_output

_SWEEP = f"{cluster.SWEEP[0]:.2f}-{cluster.SWEEP[-1]:.2f}"


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="Kmer2LTR",
        description="Classify LTR-RT boundaries and report K2P divergence between the LTRs.",
    )
    p.add_argument("input", help="multi-FASTA of putative intact LTR-RTs (.fa or .fa.gz)")
    p.add_argument("-o", "--output", default="-", help="output TSV (default: stdout)")
    p.add_argument("-u", "--mutation-rate", type=float, default=None, metavar="RATE",
                   help="neutral substitution rate per site per year; fills the "
                        "k2p_time column with years since insertion")
    p.add_argument("--cs", action="store_true",
                   help="emit a minimap2 cs string instead of an extended CIGAR")
    p.add_argument("-t", "--threads", type=int, default=1, help="worker processes")
    p.add_argument("--resume", action="store_true",
                   help="skip records already present in the output and append")
    p.add_argument("-v", "--verbose", action="store_true", help="per-step progress")
    p.add_argument("--genome", nargs="+", metavar="FASTA", default=None,
                   help="reference genome(s) the input was extracted from, plain "
                        "or gzipped. Fills the orientation and TSD columns for "
                        "records whose header carries a chrom:start-end locus, "
                        "and lets --trim-flanks correct the header of a record "
                        "stored reverse-complemented")

    ex = p.add_argument_group("extra outputs (all require -o, none by default)")
    ex.add_argument("--trim-flanks", action="store_true",
                    help="the input elements with their called flanks cut off "
                         "-> <out>.trimmed.fa")
    ex.add_argument("--ltr-cluster", action="store_true",
                    help="group elements into families by consensus LTR "
                         "-> <out>.consensus.fa and <out>.consensus_id<ID>_cluster.tsv")
    ex.add_argument("--internal-cluster", action="store_true",
                    help="the same over internal (between-LTR) regions "
                         "-> <out>.internal_id<ID>_cluster.tsv")
    ex.add_argument("--perfect-ltr-rt", nargs="+", choices=PERFECT_MODES,
                    metavar="MODE", default=None,
                    help="unmutated elements carrying one LTR copy at both ends, "
                         "as at insertion: 5p, 3p and/or consensus "
                         "-> <out>.perfect_<MODE>.fa")
    ex.add_argument("--plot", action="store_true",
                    help="K2P divergence density figure -> <out>.density.pdf")
    ex.add_argument("--min-seq-id", type=float, default=None, metavar="FLOAT",
                    help=f"cluster at this one identity instead of the {_SWEEP} sweep")

    adv = p.add_argument_group("advanced (benchmark-calibrated defaults)")
    adv.add_argument("--flank-bits", type=float, default=None,
                     help="pin the evidence in bits required to call a flank; "
                          "default is a divergence-aware schedule keyed on the "
                          "element's own estimated divergence")
    adv.add_argument("--flank-sensitivity", choices=sorted(align.FLANK_SENSITIVITY),
                     default="strict",
                     help="how much of a short flank's available evidence it must "
                          "present before a flank is called. A k-base flank can "
                          "supply at most k*alpha bits, so a flat threshold is "
                          "unreachable below k = t_bits/alpha and short flanks are "
                          "undetectable by construction. 'strict' keeps the flat "
                          "threshold and suits tightly-extracted input; the looser "
                          "settings suit input padded with genomic context, at a "
                          "real cost in false flanks (default: %(default)s)")
    adv.add_argument("--period-rule", choices=list(align.PERIOD_RULES),
                     default="best-score",
                     help="which located pair wins when a record offers more "
                          "than one. 'best-score' takes the highest-scoring "
                          "alignment of the two search windows, which is what "
                          "every calibrated default here was measured under. "
                          "'outermost' instead takes the pair reaching furthest "
                          "towards both termini among those that stay "
                          "significant and still leave an internal region -- use "
                          "it when LTRs carrying a tandem array are being called "
                          "too far in, since the aligner locks onto a register "
                          "shifted by whole array units (default: %(default)s)")
    adv.add_argument("--min-bitscore", type=float, default=None,
                     help="minimum alignment bit score to report a pair")
    adv.add_argument("--max-window", type=int, default=None,
                     help="cap on the prefix/suffix search window in bp")
    adv.add_argument("--tsd-anchor", type=float, default=0.0, metavar="BITS",
                     help="treat a target-site duplication at a record's own "
                          "termini as this many bits of evidence against calling "
                          "a flank there. Needs --genome. The default of 0 leaves "
                          "every boundary exactly where it would be without a "
                          "reference, which is what keeps the TSD an independent "
                          "check on the answer rather than an input to it "
                          "(default: %(default)s)")
    return p


def _base(output: str) -> str:
    """Stem the auxiliary outputs are named from. Raises ValueError if unusable."""
    return str(Path(output).with_suffix(""))


def _validate(args, spec: ExtraSpec) -> str | None:
    """Reject impossible combinations before any work is done. Returns a message."""
    inp = Path(args.input)
    if not inp.exists():
        return f"input file not found: {args.input}"
    if inp.is_dir():
        return f"input is a directory, not a FASTA: {args.input}"
    if args.threads < 1:
        return f"--threads must be >= 1, got {args.threads}"
    if args.mutation_rate is not None and not args.mutation_rate > 0:
        return f"-u/--mutation-rate must be > 0, got {args.mutation_rate}"
    for g in args.genome or ():
        if not Path(g).exists():
            return f"--genome file not found: {g}"
        if Path(g).is_dir():
            return f"--genome is a directory, not a FASTA: {g}"
    if args.tsd_anchor:
        if args.tsd_anchor < 0:
            return f"--tsd-anchor must be >= 0, got {args.tsd_anchor}"
        if not args.genome:
            return "--tsd-anchor needs --genome: the duplication it scores lies "\
                   "outside the record"

    wants_extra = bool(spec) or args.plot
    if wants_extra and args.output == "-":
        return "the extra-output flags name their files after -o/--output, so -o is required"
    if args.output != "-":
        try:
            base = _base(args.output)
        except ValueError:
            return f"-o/--output is not a usable path: {args.output!r}"
        # Every output is opened for WRITING before the input is read, so a
        # derived path that resolves to the input truncates the very file the
        # run is about to parse -- and the emptied stream is then deleted as
        # "received no records". Feeding a --trim-flanks output back in is the
        # obvious way to land here, so it is checked, not documented.
        # Every reference is checked alongside the input for the same reason:
        # these paths are opened for writing before anything is read.
        reads = {inp.resolve(): args.input}
        for g in args.genome or ():
            reads[Path(g).resolve()] = g
        for label, path in (("-o/--output", Path(args.output)),
                            *((f"the {k} output", v)
                              for k, v in stream_paths(base, spec).items())):
            if path.exists() and path.resolve() in reads:
                return (f"{label} would be written to an input file "
                        f"({reads[path.resolve()]}); it would be destroyed")
    if wants_extra and args.resume:
        # The TSV resumes on its data-line count; the auxiliary FASTAs hold only
        # the passing subset, so that count cannot tell us where they stopped.
        # Silently appending would duplicate or drop records in them.
        return "--resume cannot be combined with the extra-output flags"

    clustering = args.ltr_cluster or args.internal_cluster
    if clustering and not cluster.available():
        return ("clustering needs mmseqs on PATH "
                "(mamba install -c bioconda mmseqs2)")
    if args.min_seq_id is not None:
        if not 0.0 < args.min_seq_id <= 1.0:
            return f"--min-seq-id must be in (0, 1], got {args.min_seq_id}"
        # Identities are used to two decimals in both the mmseqs argument and
        # the file name, so a finer value would silently be clustered at
        # something other than what was asked for -- and named for the value it
        # was not run at. 0.999 -> 1.00, which is one cluster per sequence.
        if round(args.min_seq_id, 2) != args.min_seq_id:
            return (f"--min-seq-id is used to two decimal places; "
                    f"{args.min_seq_id} would run at {args.min_seq_id:.2f}")
        if not clustering:
            return "--min-seq-id has no effect without --ltr-cluster or --internal-cluster"
    if args.plot and importlib.util.find_spec("matplotlib") is None:
        return "--plot needs matplotlib (mamba install -c conda-forge matplotlib)"
    return None


def _cluster_and_report(fasta, args, drop: bool) -> bool:
    """Cluster one auxiliary FASTA, reporting what was written. True if all ran."""
    wanted = len(cluster.identities(args.min_seq_id))
    tsvs = cluster.cluster(fasta, args.threads, args.min_seq_id, args.verbose)
    complete = len(tsvs) == wanted
    if tsvs:
        # Name every table when the sweep is incomplete: an elided ".. " range
        # over a gapped set reads as a full sweep of identities that did not run.
        listed = (f"{tsvs[0].name} .. {tsvs[-1].name}" if complete and len(tsvs) > 1
                  else ", ".join(t.name for t in tsvs))
        print(f"Kmer2LTR: wrote {len(tsvs)}/{wanted} cluster table(s): {listed}",
              file=sys.stderr)
    if drop and complete:
        Path(fasta).unlink(missing_ok=True)
    elif drop:
        # Keep the input to a step that did not finish: it is the expensive half
        # of the run, and re-deriving it means re-aligning everything. Dropping
        # it here would make the missing identities unrecoverable.
        print(f"Kmer2LTR: kept {fasta} (clustering did not complete)", file=sys.stderr)
    return complete


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    spec = ExtraSpec(consensus=args.ltr_cluster, internal=args.internal_cluster,
                     trimmed=args.trim_flanks,
                     perfect=tuple(args.perfect_ltr_rt or ()))
    err = _validate(args, spec)
    if err:
        print(f"Kmer2LTR: error: {err}", file=sys.stderr)
        return 2

    # Thread tuning values explicitly. Assigning align.T_BITS / align.W0 would
    # be a silent no-op: both are already bound as default arguments at def time.
    classify_kw = {"flank_sensitivity": args.flank_sensitivity,
                   "mutation_rate": args.mutation_rate,
                   "period_rule": args.period_rule}
    gkw = {"genome": args.genome, "gopt": Options(anchor=args.tsd_anchor)}
    if args.flank_bits is not None:
        classify_kw["t_bits"] = args.flank_bits
    if args.max_window is not None:
        classify_kw["w0"] = args.max_window
    if args.min_bitscore is not None:
        classify_kw["min_bitscore"] = args.min_bitscore

    skip, has_header = 0, False
    if args.resume:
        if args.output == "-":
            print("Kmer2LTR: error: --resume requires -o/--output", file=sys.stderr)
            return 2
        skip, keep_bytes, has_header = scan_output(args.output)
        # Drop any partial final line before appending. A run killed mid-write
        # leaves one, and appending after it concatenates two records into a row
        # that no reader can parse -- while `skip` has already stepped past the
        # record that fragment belongs to, losing it.
        if has_header and Path(args.output).stat().st_size != keep_bytes:
            with open(args.output, "r+b") as fh:
                fh.truncate(keep_bytes)
            print(f"Kmer2LTR: discarded a partial final line before resuming",
                  file=sys.stderr)
        if args.verbose:
            print(f"resuming: {skip} records already done", file=sys.stderr)

    inp = Path(args.input)
    base = _base(args.output) if args.output != "-" else ""
    writer = ExtraWriter(base, spec) if spec else None
    t0 = time.time()
    print(f"Kmer2LTR: reading {inp}", file=sys.stderr)
    try:
        if args.output == "-":
            n = run(inp, sys.stdout, args.threads, args.cs, skip, args.verbose,
                    resuming=args.resume, spec=None, writer=None,
                    **gkw, **classify_kw)
        else:
            # append whenever resuming -- even at skip == 0, where the header was already
            # flushed but no record finished. The header is suppressed on the strength of
            # one actually being THERE (`has_header`), not merely of --resume being
            # passed: resuming onto a file that does not exist yet would otherwise write
            # a headerless TSV, which every csv.DictReader downstream silently mistakes
            # its first data row for.
            mode = "a" if args.resume else "w"
            with open(args.output, mode) as fh:
                n = run(inp, fh, args.threads, args.cs, skip, args.verbose,
                        resuming=has_header, spec=spec or None, writer=writer,
                        **gkw, **classify_kw)
    finally:
        if writer is not None:
            writer.close()
    print(f"Kmer2LTR: {n} records in {time.time() - t0:.1f}s", file=sys.stderr)

    rc = 0
    if writer is not None:
        for key, path in sorted(writer.paths.items()):
            print(f"Kmer2LTR: wrote {path} ({writer.counts[key]} records)", file=sys.stderr)
        for flag, key, drop in ((args.ltr_cluster, "consensus", False),
                                (args.internal_cluster, "internal", True)):
            if not flag:
                continue
            if key not in writer.paths:
                print(f"Kmer2LTR: warning: no {key} sequence to cluster", file=sys.stderr)
                rc = 1
            elif not _cluster_and_report(writer.paths[key], args, drop=drop):
                # A clustering that did not finish must be visible to `$?`, or a
                # shell pipeline carries on into tables that are not there.
                rc = 1
    if args.plot:
        from .plot import density_plot
        pdf = f"{base}.density.pdf"
        if density_plot(args.output, pdf, args.mutation_rate):
            print(f"Kmer2LTR: wrote {pdf}", file=sys.stderr)
    return rc
