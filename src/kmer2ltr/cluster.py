"""mmseqs2 clustering of the per-element sequences Kmer2LTR derives.

Clustering consensus LTRs groups elements into families; clustering internal
regions is a mostly independent check on that grouping. Both run the same
command with the same tuned parameters.
"""
from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# --min-seq-id sets the DEPTH of the clustering rather than its correctness:
# low values lump at the lineage level, high values split toward recent bursts.
# Different downstream questions want different depths, so the whole sweep is
# kept rather than one value picked. 1.00 is degenerate (one cluster per
# sequence); 0.98 is the practical near-identical tier.
SWEEP = (0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.98)

# Grid-searched on Arabidopsis LTR annotations against two objectives at once:
# minimise single-copy clusters, and minimise family mixing (a Tork clustering
# with an Ogre is wrong; a Tork with an Unknown is not). mmseqs beat cd-hit-est
# and a wavefront pipeline, and the consensus LTR beat both the full-length
# element and the 5' LTR alone as input. -c 0.4 overclusters and 0.6
# underclusters, on internal regions as well as on LTRs.
COVERAGE = 0.5
_FIXED = ("--cov-mode", "0", "--cluster-mode", "1", "--mask", "0", "-s", "7.5")


def available() -> bool:
    """True if mmseqs is on PATH."""
    return shutil.which("mmseqs") is not None


def identities(min_seq_id: float | None) -> tuple[float, ...]:
    """The identities `cluster` will attempt."""
    return (min_seq_id,) if min_seq_id is not None else SWEEP


def cluster(fasta, threads: int = 1, min_seq_id: float | None = None,
            verbose: bool = False) -> list[Path]:
    """Cluster a FASTA over the identity sweep; return the cluster TSVs written.

    One `<stem>_id<ID>_cluster.tsv` per identity, or a single file when
    `min_seq_id` pins one. mmseqs' own scratch -- the tmp directory,
    `_all_seqs.fasta` and `_rep_seq.fasta` -- is removed after every run,
    including the failing ones.

    A failed identity is reported and skipped rather than aborting the sweep:
    the other identities are independent of it, and the run has already paid
    for the alignment.
    """
    fasta = Path(fasta)
    if not fasta.exists() or fasta.stat().st_size == 0:
        print(f"Kmer2LTR: warning: nothing to cluster in {fasta}", file=sys.stderr)
        return []

    written: list[Path] = []
    for seq_id in identities(min_seq_id):
        tag = f"{seq_id:.2f}"
        base = fasta.with_suffix("")                     # elements.consensus
        prefix = base.with_name(f"{base.name}_id{tag}")  # elements.consensus_id0.70
        tsv = Path(f"{prefix}_cluster.tsv")
        # Remove any table from an earlier run FIRST. If this identity fails,
        # the stale file would otherwise survive beside six fresh ones and the
        # directory would look like a complete sweep of the current input.
        tsv.unlink(missing_ok=True)
        tmp = Path(tempfile.mkdtemp(prefix=f".{fasta.stem}_mmseqs_{tag}_",
                                    dir=str(fasta.parent)))
        cmd = ["mmseqs", "easy-cluster", str(fasta), str(prefix), str(tmp),
               "--min-seq-id", tag, "-c", f"{COVERAGE:g}", *_FIXED,
               "--threads", str(threads)]
        if verbose:
            print("  " + " ".join(cmd), file=sys.stderr)
        try:
            subprocess.run(cmd, check=True,
                           stdout=None if verbose else subprocess.DEVNULL,
                           stderr=None if verbose else subprocess.DEVNULL)
        except subprocess.CalledProcessError as exc:
            print(f"Kmer2LTR: warning: mmseqs failed at --min-seq-id {tag} "
                  f"(exit {exc.returncode}); continuing", file=sys.stderr)
            continue
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
            for aux in ("_all_seqs.fasta", "_rep_seq.fasta"):
                Path(f"{prefix}{aux}").unlink(missing_ok=True)

        if tsv.exists():
            written.append(tsv)
        else:
            print(f"Kmer2LTR: warning: mmseqs wrote no {tsv.name}", file=sys.stderr)
    return written
