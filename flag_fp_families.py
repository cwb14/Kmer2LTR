#!/usr/bin/env python3
"""Flag false-positive LTR-RT families, clean depth TSVs, and conditionally mask.

Stage A (always): two-gate FP-family filter from Kmer2LTR cluster tables, writing
  <prefix>.family_scores.tsv, <prefix>.fp_families.tsv, <prefix>.fp_LTRs.fa and a
  <prefix>.fp_structure.pdf (identical to flag_fp_families.py).
Stage B (always): for each --domains-tsv, write an FP-cleaned copy
  (<name>_ltr.tsv -> <name>_clean_ltr.tsv): FP rows removed and dangling nest
  cross-references scrubbed in both directions.
Stage C (only if FP fraction > --fp-mask-threshold): mmseqs easy-cluster the FP
  LTRs, then dc-megablast them against --genome and hard-mask hits to N.

FP fraction = elements in false_positive families / total consensus-cluster members.
In maize 424/132931 (0.3%) -> Stage C skipped; in dog 6284/8244 (76%) -> masking then
re-annotation lifts precision from ~0.17 to ~0.96.

awk '/^>/{printf "\n%s\n",$0;next}{gsub(/[^ACGTacgt]/,"");printf "%s",$0}END{print ""}' ./dog_depth*_ltr.fa > dog_all_ltr.fa
python Kmer2LTR/Kmer2LTR.py -i dog_all_ltr.fa -o dog_all_ltr --ltr-cluster --internal-cluster -p 200 --min-seq-id 0.75
python Kmer2LTR/flag_and_mask_fp.py --consensus-cluster dog_all_ltr.consensus_id0.75_cluster.tsv --internal-cluster dog_all_ltr.internal_id0.75_cluster.tsv --ltr-fasta dog_all_ltr.consensus.fa --domains-tsv dog_depth*[0-9]_ltr.tsv -o dog_fpcheck --genome dog.fa --threads 100

# Original (ie, no flag_and_mask_fp.py; dog_depth0_ltr.tsv):
Precision = 0.1747   (TP_pred 524 / predictions 2999)
Recall    = 0.6429   (TP_truth 524 / truths 815)   [= sensitivity]
F1        = 0.2748
FDR       = 0.8253   (1 - precision)

# Cleaned (flag_and_mask_fp.py used to purge those suspicious ones; dog_depth0_clean_ltr.tsv):
Precision = 0.9424   (TP_pred 524 / predictions 556)
Recall    = 0.6429   (TP_truth 524 / truths 815)   [= sensitivity]
F1        = 0.7644
FDR       = 0.0576   (1 - precision)

# Using the flag_and_mask_fp.py hardmasked genome to re-run LTR-RT detection; dog_hardmaskedFP_depth0_ltr.tsv)
Precision = 0.9586   (TP_pred 717 / predictions 748)
Recall    = 0.8798   (TP_truth 717 / truths 815)   [= sensitivity]
F1        = 0.9175
FDR       = 0.0414   (1 - precision)
"""

import argparse
import math
import multiprocessing as mp
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass, field

IGNORE_DEFAULT = frozenset({"unknown", "mixture"})

# Two-gate filter defaults.
DEFAULT_MIN_MEMBERS = 10
DEFAULT_MAX_RECON = 0.51
DEFAULT_MAX_DOMINANCE = 0.51

# -----------------------------------------------------------------------------
# Parsing & family grouping
# -----------------------------------------------------------------------------
def parse_clusters(path):
    """mmseqs `rep<TAB>member` TSV -> {member_id: rep_id}. Members are unique."""
    member2rep = {}
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                raise ValueError(
                    f"{path}:{lineno}: expected 2 tab-separated columns, got {len(parts)}"
                )
            rep, mem = parts[0], parts[1]
            if mem in member2rep:
                print(f"[WARN] duplicate member {mem} at {path}:{lineno}; keeping first",
                      file=sys.stderr)
                continue
            member2rep[mem] = rep
    if not member2rep:
        raise ValueError(f"{path}: no cluster rows parsed")
    return member2rep


def parse_label(element_id):
    """'chrom:start-end#Class/Super/Clade' -> (cls, superfamily, clade).
    Missing/empty tokens become 'unknown'."""
    tag = element_id.split("#", 1)[1] if "#" in element_id else ""
    bits = tag.split("/")
    cls = bits[0] if len(bits) > 0 and bits[0] else "unknown"
    sup = bits[1] if len(bits) > 1 and bits[1] else "unknown"
    clade = bits[2] if len(bits) > 2 and bits[2] else "unknown"
    return cls, sup, clade


def is_real_clade(sup, clade, ignore=IGNORE_DEFAULT):
    """A member carries an informative clade iff neither token is ignored."""
    return sup not in ignore and clade not in ignore


def group_families(member2rep):
    """{member: rep} -> {rep: [members]}."""
    fam = defaultdict(list)
    for mem, rep in member2rep.items():
        fam[rep].append(mem)
    return dict(fam)


# -----------------------------------------------------------------------------
# Family metrics
# -----------------------------------------------------------------------------
def reconstitution(members, internal_map):
    """Fraction of members whose internal cluster is shared by >=2 family members."""
    n = len(members)
    if n == 0:
        return 0.0
    counts = Counter(internal_map[m] for m in members)  # within-family sizes
    co_clustered = sum(c for c in counts.values() if c >= 2)
    return co_clustered / n


def real_clade_counts(members, ignore=IGNORE_DEFAULT):
    counts = Counter()
    for m in members:
        _, sup, clade = parse_label(m)
        if is_real_clade(sup, clade, ignore):
            counts[clade] += 1
    return counts


def dominance(members, ignore=IGNORE_DEFAULT):
    """(largest real-clade count) / n; 0 if no real-clade member."""
    n = len(members)
    if n == 0:
        return 0.0
    counts = real_clade_counts(members, ignore)
    return (max(counts.values()) / n) if counts else 0.0


def clade_entropy(members, ignore=IGNORE_DEFAULT):
    """Shannon entropy (bits) of the real-clade distribution."""
    counts = real_clade_counts(members, ignore)
    total = sum(counts.values())
    if total == 0:
        return 0.0
    ent = -sum((v / total) * math.log2(v / total) for v in counts.values())
    return ent if ent else 0.0  # collapse -0.0 to 0.0


def cross_superfamily(members, ignore=IGNORE_DEFAULT):
    """True iff >=2 distinct real superfamilies present with the 2nd having >=2
    members (a lone cross-superfamily contaminant does not qualify)."""
    sup_counts = Counter()
    for m in members:
        _, sup, clade = parse_label(m)
        if is_real_clade(sup, clade, ignore):
            sup_counts[sup] += 1
    if len(sup_counts) < 2:
        return False
    top2 = sorted(sup_counts.values(), reverse=True)[:2]
    return top2[1] >= 2


@dataclass
class FamilyMetrics:
    rep: str
    n: int
    reconstitution: float
    dominance: float
    entropy: float
    cross_superfamily: bool
    n_unknown: int
    n_mixture: int
    clade_composition: dict = field(default_factory=dict)
    verdict: str = "safe"  # 'safe' | 'recovered' | 'false_positive'


def compute_family_metrics(rep, members, internal_map, ignore=IGNORE_DEFAULT):
    comp = Counter()
    n_unknown = n_mixture = 0
    for m in members:
        _, sup, clade = parse_label(m)
        comp[f"{sup}/{clade}"] += 1
        if sup == "unknown" or clade == "unknown":
            n_unknown += 1
        if clade == "mixture":
            n_mixture += 1
    return FamilyMetrics(
        rep=rep,
        n=len(members),
        reconstitution=reconstitution(members, internal_map),
        dominance=dominance(members, ignore),
        entropy=clade_entropy(members, ignore),
        cross_superfamily=cross_superfamily(members, ignore),
        n_unknown=n_unknown,
        n_mixture=n_mixture,
        clade_composition=dict(comp),
    )


def classify_family(fm, max_recon=DEFAULT_MAX_RECON, max_dominance=DEFAULT_MAX_DOMINANCE):
    """Two-gate verdict for a family already past the size gate. Returns one of
    'safe', 'recovered', 'false_positive'.

      Gate 1 (dominance): a family with a dominant real clade
        (dominance >= max_dominance) has coherent labels and PASSES -> 'safe'.
      Gate 2 (reconstitution): a family that FAILS gate 1 gets a second chance --
        if its internals reconstitute (reconstitution > max_recon) it is
        'recovered'; otherwise (reconstitution <= max_recon) it is a
        'false_positive'.
    """
    if fm.dominance >= max_dominance:
        return "safe"                # gate 1 passed: coherent labels
    if fm.reconstitution > max_recon:
        return "recovered"           # gate 2 rescue: internals reconstitute
    return "false_positive"          # failed both gates


def classify_fp(fm, min_members=DEFAULT_MIN_MEMBERS, max_recon=DEFAULT_MAX_RECON,
                max_dominance=DEFAULT_MAX_DOMINANCE):
    """True iff the family is a false positive (size gate + two-gate filter)."""
    return (fm.n >= min_members
            and classify_family(fm, max_recon, max_dominance) == "false_positive")


# -----------------------------------------------------------------------------
# Report TSVs
# -----------------------------------------------------------------------------
_SCORE_HEADER = ["rep", "n", "reconstitution", "dominance", "entropy",
                 "cross_superfamily", "n_unknown", "n_mixture",
                 "verdict", "clade_composition"]


def format_composition(comp):
    """{'Gypsy/Tekay':4,...} -> 'Gypsy/Tekay:4;...' sorted by descending count."""
    return ";".join(f"{k}:{v}" for k, v in sorted(comp.items(), key=lambda x: (-x[1], x[0])))


def _row(fm, include_verdict):
    cells = [fm.rep, str(fm.n), f"{fm.reconstitution:.4f}", f"{fm.dominance:.4f}",
             f"{fm.entropy:.4f}", str(fm.cross_superfamily),
             str(fm.n_unknown), str(fm.n_mixture)]
    if include_verdict:
        cells.append(fm.verdict)  # safe | recovered | false_positive
    cells.append(format_composition(fm.clade_composition))
    return "\t".join(cells)


def write_family_scores(path, metrics):
    # false positives first, then recovered, then safe; within a tier by dominance, size.
    rank = {"false_positive": 0, "recovered": 1, "safe": 2}
    ordered = sorted(metrics, key=lambda x: (rank.get(x.verdict, 3), x.dominance, -x.n))
    with open(path, "w") as fh:
        fh.write("\t".join(_SCORE_HEADER) + "\n")
        for fm in ordered:
            fh.write(_row(fm, include_verdict=True) + "\n")


def write_fp_families(path, metrics):
    fp = sorted((fm for fm in metrics if fm.verdict == "false_positive"),
                key=lambda x: (x.dominance, -x.n))
    header = [h for h in _SCORE_HEADER if h != "verdict"]
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for fm in fp:
            fh.write(_row(fm, include_verdict=False) + "\n")


# -----------------------------------------------------------------------------
# Output FASTA (indexed extraction)
# -----------------------------------------------------------------------------
def write_fp_fasta(path, member_ids, ltr_fasta_path):
    """Extract each member's consensus LTR from ltr_fasta_path (indexed) -> FASTA.
    Returns (n_written, n_missing). Missing IDs are warned and skipped."""
    from pyfaidx import Fasta
    fa = Fasta(ltr_fasta_path)
    n_written = n_missing = 0
    with open(path, "w") as out:
        for mid in member_ids:
            if mid not in fa:
                print(f"[WARN] {mid} absent from {ltr_fasta_path}; skipping", file=sys.stderr)
                n_missing += 1
                continue
            seq = str(fa[mid])
            out.write(f">{mid}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")
            n_written += 1
    return n_written, n_missing


# -----------------------------------------------------------------------------
# Plot-data loading
# -----------------------------------------------------------------------------
def load_results(path):
    """Kmer2LTR results table -> {element_id: (ltr_len, k2p_or_None)}.
    Columns (0-based): 0=id, 1=LTR_LEN, 10=K2P-dist."""
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 16:
                continue
            try:
                ltr_len = int(float(f[1]))
            except ValueError:
                continue
            try:
                k2p = float(f[10])
            except (ValueError, IndexError):
                k2p = None
            out[f[0]] = (ltr_len, k2p)
    return out


def load_lengths(element_fasta_path):
    """{element_id: sequence_length} from the .fai (built via pyfaidx if absent).
    NB: uses .fai length, never the header coordinates."""
    fai = element_fasta_path + ".fai"
    if not os.path.exists(fai):
        from pyfaidx import Faidx
        Faidx(element_fasta_path)  # builds <path>.fai
    lengths = {}
    with open(fai) as fh:
        for line in fh:
            parts = line.split("\t")
            lengths[parts[0]] = int(parts[1])
    return lengths


# -----------------------------------------------------------------------------
# Protein-domain / nesting support (optional; from synLTR depth TSVs)
# Coordinate/collapse/strand logic mirrors synLTR/module2/ltrharvest_plot_struct.py.
# -----------------------------------------------------------------------------
# Domain box colors keyed by normalized gene name (LTR handled separately).
FEATURE_COLORS = {
    "GAG": "#fc850e", "PROT": "#f80bfb", "RT": "#0808f7", "RH": "#fc1413",
    "INT": "#05fc09", "CH": "#20b5f4", "CHD": "#20b5f4", "CHDCR": "#7a001a",
    "ARH": "#70403c", "ENDO": "#a67c52",
}
FEATURE_ORDER = ["GAG", "PROT", "RT", "RH", "INT", "CH", "CHD", "CHDCR", "ARH", "ENDO"]
DOMAIN_DEFAULT_COLOR = "#AAAAAA"

# Canonical 5'->3' internal-domain order per superfamily, for strand inference.
_CANON_RANK = {
    "copia": {"GAG": 0, "PROT": 1, "INT": 2, "RT": 3, "RH": 4},
    "gypsy": {"GAG": 0, "PROT": 1, "RT": 2, "RH": 3, "INT": 4},
}
_CANON_CORE = {"GAG": 0, "PROT": 1, "RT": 2, "RH": 3}


def normalize_protein_name(raw):
    p = re.sub(r"[^A-Z0-9]+", "", raw.strip().upper())
    if p == "INTEGRASE":
        return "INT"
    if p == "PROTEASE":
        return "PROT"
    if p == "RNASEH":
        return "RH"
    if p in ("ARH", "ARNASEH"):
        return "ARH"
    return p  # CH* and already-normalized names pass through


def parse_domains_field(field):
    """'gene|clade@gStart-gEnd;...' or '.' -> [(gene, gStart, gEnd)] genomic."""
    out = []
    if not field or field == ".":
        return out
    for tok in field.split(";"):
        m = re.match(r"^([^|@]+)\|[^@]*@(\d+)-(\d+)$", tok.strip())
        if not m:
            continue
        gs, ge = int(m.group(2)), int(m.group(3))
        if ge < gs:
            gs, ge = ge, gs
        out.append((normalize_protein_name(m.group(1)), gs, ge))
    return out


def parse_insertions(nest_raw):
    """'nest-outer:chrom:s-e;...' or '.' -> [(gStart, gEnd)] genomic (nested children)."""
    out = []
    if not nest_raw or nest_raw == ".":
        return out
    for tok in nest_raw.split(";"):
        tok = tok.strip()
        if not tok.startswith("nest-outer:"):
            continue
        m = re.match(r"^(.+):(\d+)-(\d+)$", tok[len("nest-outer:"):])
        if m:
            out.append((int(m.group(2)), int(m.group(3))))
    return out


def merge_intervals(iv):
    """Union of 1-based inclusive intervals -> sorted disjoint list."""
    if not iv:
        return []
    xs = sorted((min(s, e), max(s, e)) for s, e in iv)
    out = [list(xs[0])]
    for s, e in xs[1:]:
        if s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def infer_strand_from_domains(raw_domains, superfamily):
    """'+'/'-'/'?' from genomic order of ranked domains vs canonical 5'->3' order."""
    rank = _CANON_RANK.get((superfamily or "").lower(), _CANON_CORE)
    pos = {}
    for gene, gs, _ge in raw_domains:
        if gene in rank and (gene not in pos or gs < pos[gene]):
            pos[gene] = gs
    if len(pos) < 2:
        return "?"
    genes = list(pos)
    plus = minus = 0
    for i in range(len(genes)):
        for j in range(i + 1, len(genes)):
            gi, gj = genes[i], genes[j]
            if rank[gi] == rank[gj] or pos[gi] == pos[gj]:
                continue
            if (rank[gi] < rank[gj]) == (pos[gi] < pos[gj]):
                plus += 1
            else:
                minus += 1
    return "+" if plus > minus else ("-" if minus > plus else "?")


@dataclass
class DepthElement:
    start: int          # genomic start (from element ID)
    span: int           # genomic span (end - start), includes nested insertions
    ltr_len: int
    k2p: float
    superfamily: str
    raw_domains: list   # [(gene, gStart, gEnd)] genomic
    insertions: list    # [(gStart, gEnd)] genomic nested children


def load_depth_tsvs(paths):
    """Parse synLTR depth TSV(s) -> {element_id: DepthElement}. An element lives in
    exactly one depth bucket; first occurrence wins. Header/comment rows skipped."""
    out = {}
    for path in paths:
        with open(path) as fh:
            for line in fh:
                if not line or line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 3:
                    continue
                eid = f[0]
                if eid in out:
                    continue
                m = re.match(r"^(.+):(\d+)-(\d+)#", eid)
                if not m:
                    continue
                s0, e0 = int(m.group(2)), int(m.group(3))
                span = e0 - s0
                if span <= 0:
                    continue
                try:
                    ltr_len = int(float(f[1]))
                except (ValueError, IndexError):
                    ltr_len = 0
                try:
                    k2p = float(f[10])
                except (ValueError, IndexError):
                    k2p = None
                _, sup, _ = parse_label(eid)
                out[eid] = DepthElement(
                    start=s0, span=span, ltr_len=ltr_len, k2p=k2p, superfamily=sup,
                    raw_domains=parse_domains_field(f[-2]),
                    insertions=parse_insertions(f[-1]),
                )
    return out


@dataclass
class DrawSpec:
    length: int          # collapsed (nest-in-excised) in domain mode, else .fai length
    ltr_len: int
    k2p: float
    domains: list        # [(gene, start, end)] 1-based element-relative, collapsed
    ins_markers: list    # collapsed x positions of excised insertions


def _domain_spec(de):
    """DepthElement -> DrawSpec: element-relative, strand-normalized 5'->3',
    with nested insertions excised (collapsed) and domain coords remapped."""
    L = de.span

    def rel(g):
        return max(1, min(g - de.start + 1, L))

    prots = [(gene, rel(gs), rel(ge)) for gene, gs, ge in de.raw_domains]
    ins = [(min(rel(gs), rel(ge)), max(rel(gs), rel(ge))) for gs, ge in de.insertions]

    if infer_strand_from_domains(de.raw_domains, de.superfamily) == "-":
        prots = [(g, min(L - e + 1, L - s + 1), max(L - e + 1, L - s + 1)) for g, s, e in prots]
        ins = [(min(L - e + 1, L - s + 1), max(L - e + 1, L - s + 1)) for s, e in ins]

    merged = merge_intervals(ins)
    excised = sum(e - s + 1 for s, e in merged)

    def remap(x):
        shift = 0
        for s, e in merged:
            if e < x:
                shift += (e - s + 1)
            elif s <= x <= e:
                shift += (x - s)
        return x - shift

    clen = max(1, L - excised)
    cdom = []
    for g, s, e in prots:
        ns, ne = remap(s), remap(e)
        if ne < ns:
            ns, ne = ne, ns
        cdom.append((g, max(1, ns), max(1, ne)))
    markers = [min(remap(s), clen) for s, e in merged]
    return DrawSpec(length=clen, ltr_len=min(de.ltr_len, clen), k2p=de.k2p,
                    domains=cdom, ins_markers=markers)


def build_domain_specs(depth_elements):
    """{id: DepthElement} -> {id: DrawSpec} (collapsed, strand-normalized, with domains)."""
    return {eid: _domain_spec(de) for eid, de in depth_elements.items()}


def build_plain_specs(member_ids, results, lengths):
    """{id: DrawSpec} with no domains: length from .fai, LTR/K2P from results table."""
    specs = {}
    for m in member_ids:
        ltr_len, k2p = results.get(m, (0, None))
        specs[m] = DrawSpec(length=lengths.get(m, 0), ltr_len=ltr_len, k2p=k2p,
                            domains=[], ins_markers=[])
    return specs


def select_borderline_families(metrics, max_recon, max_dominance, n):
    """Non-FP families that sit closest to becoming a false positive, ranked by how
    small a threshold move would flip them. Returns [(gap, gate, fm)] (closest first).
    gate='dominance' -> family only just PASSED the dominance gate; raise
      --max-dominance to flag it.
    gate='recon'     -> family failed dominance but was only just RESCUED by recon;
      raise --max-recon to flag it."""
    cands = []
    for fm in metrics:
        if fm.verdict == "false_positive":
            continue
        passed_dominance = fm.dominance >= max_dominance
        rescued_by_recon = fm.reconstitution > max_recon
        if passed_dominance and fm.reconstitution <= max_recon:
            # cleared gate 1 only; would be FP if the dominance gate were higher
            cands.append((fm.dominance - max_dominance, "dominance", fm))
        elif (not passed_dominance) and rescued_by_recon:
            # failed gate 1, rescued at gate 2; would be FP if the recon bar were higher
            cands.append((fm.reconstitution - max_recon, "recon", fm))
        # cleared both comfortably -> not near the FP boundary; skip
    cands.sort(key=lambda x: (x[0], -x[2].n))
    return cands[:max(0, n)]


# -----------------------------------------------------------------------------
# Plotting (matplotlib Agg + PdfPages)
# -----------------------------------------------------------------------------
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

_LTR_COLOR = "#333333"
_ORPHAN_COLOR = "#BBBBBB"
# Okabe-Ito colorblind-safe palette for shared internal clusters.
_SHARED_COLORS = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9",
                  "#D55E00", "#F0E442", "#000000"]


def internal_color_map(members, internal_map):
    """member -> color: gray for within-family singleton (orphan), a cycled
    color per shared internal cluster (matching across the sharing members)."""
    counts = Counter(internal_map[m] for m in members)
    cluster_color, color, ci = {}, {}, 0
    for m in members:
        ic = internal_map[m]
        if counts[ic] >= 2:
            if ic not in cluster_color:
                cluster_color[ic] = _SHARED_COLORS[ci % len(_SHARED_COLORS)]
                ci += 1
            color[m] = cluster_color[ic]
        else:
            color[m] = _ORPHAN_COLOR
    return color


def draw_member(ax, y, spec, internal_color, height=0.65):
    """One element from a DrawSpec: outline bar, colored internal, dark LTR
    end-boxes, protein-domain boxes, and carets marking excised nested insertions."""
    L = max(1, spec.length)
    lt = min(max(0, spec.ltr_len), L / 2)
    ax.add_patch(plt.Rectangle((lt, y - height / 2), max(0, L - 2 * lt), height,
                               facecolor=internal_color, edgecolor="none"))
    for x0 in (0, L - lt):
        ax.add_patch(plt.Rectangle((x0, y - height / 2), lt, height,
                                   facecolor=_LTR_COLOR, edgecolor="black", linewidth=0.5))
    for gene, s, e in spec.domains:
        ax.add_patch(plt.Rectangle((max(0, s - 1), y - height / 2), max(1, e - s + 1), height,
                                   facecolor=FEATURE_COLORS.get(gene, DOMAIN_DEFAULT_COLOR),
                                   edgecolor="black", linewidth=0.4))
    for mx in spec.ins_markers:
        ax.plot([mx], [y - height / 2 - 0.05], marker="^", markersize=3,
                color="#222222", linestyle="none")
    ax.add_patch(plt.Rectangle((0, y - height / 2), L, height, fill=False,
                               edgecolor="black", linewidth=0.8))


def plot_family_page(pdf, fm, members, internal_map, specs, note=""):
    def _len(m):
        return specs[m].length if m in specs else 0
    counts = Counter(internal_map[m] for m in members)
    # shared clusters grouped together, orphans size-sorted; longest first throughout
    members = sorted(members, key=lambda m: (
        counts[internal_map[m]] < 2,
        internal_map[m] if counts[internal_map[m]] >= 2 else "",
        -_len(m)))
    color = internal_color_map(members, internal_map)
    n = len(members)
    max_len = max((_len(m) for m in members), default=1) or 1
    fig, ax = plt.subplots(figsize=(12, max(2.5, 0.22 * n + 1.7)))
    domains_present = set()
    for i, m in enumerate(members):
        y = n - i
        spec = specs.get(m, DrawSpec(0, 0, None, [], []))
        draw_member(ax, y, spec, color[m])
        for g, _s, _e in spec.domains:
            domains_present.add(g)
        _, sup, clade = parse_label(m)
        k2p_txt = "NA" if spec.k2p is None else f"{spec.k2p * 100:.1f}%"
        ax.text(-0.01 * max_len, y, f"{m.split('#', 1)[0]}  {sup}/{clade}  {k2p_txt}",
                ha="right", va="center", fontsize=6)
    prefix = f"[{note}]  " if note else ""
    ax.set_title(
        f"{prefix}{fm.rep.split('#', 1)[0]}   n={fm.n}  recon={fm.reconstitution:.2f}  "
        f"dom={fm.dominance:.2f}  entropy={fm.entropy:.2f}  xSuper={fm.cross_superfamily}\n"
        f"{format_composition(fm.clade_composition)}", fontsize=8)
    ax.set_xlim(-0.30 * max_len, max_len * 1.02)
    ax.set_ylim(0, n + 1)
    ax.set_yticks([])
    ax.set_xlabel("Position [bp]")
    handles = [Patch(facecolor=_LTR_COLOR, label="LTR"),
               Patch(facecolor=_ORPHAN_COLOR, label="internal: orphan (singleton)"),
               Patch(facecolor=_SHARED_COLORS[0], label="internal: shared cluster")]
    for g in FEATURE_ORDER:
        if g in domains_present:
            handles.append(Patch(facecolor=FEATURE_COLORS[g], label=g))
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1.0),
              frameon=False, fontsize=7)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def plot_summary_page(pdf, all_metrics, fp_metrics, thresholds):
    n_safe = sum(1 for fm in all_metrics if fm.verdict == "safe")
    n_rec = sum(1 for fm in all_metrics if fm.verdict == "recovered")
    n_fp = sum(1 for fm in all_metrics if fm.verdict == "false_positive")
    fig, ax = plt.subplots(figsize=(12, 9))
    ax.axis("off")
    lines = [
        "False-positive LTR-RT families - two-gate structure report", "",
        f"gate 1 (dominance): a family with a dominant real clade "
        f"(dominance >= {thresholds['max_dominance']}) PASSES -> kept as true positive.",
        f"gate 2 (reconstitution): a family that fails gate 1 is RECOVERED if its internals "
        f"reconstitute (reconstitution > {thresholds['max_recon']}),",
        "                         otherwise it is a FALSE POSITIVE.",
        f"(size gate: only families with >= {thresholds['min_members']} members are analyzed)", "",
        f"analyzed families: {len(all_metrics)}     safe (passed dominance): {n_safe}     "
        f"recovered (recon rescue): {n_rec}     false-positive: {n_fp}  "
        f"({sum(fm.n for fm in fp_metrics)} elements)", "",
        "false-positive families:",
        f"  {'rep_locus':38s} {'n':>4} {'recon':>6} {'dom':>5} {'ent':>5} {'xSup':>5}  composition",
    ]
    for fm in sorted(fp_metrics, key=lambda x: (x.dominance, -x.n)):
        lines.append(
            f"  {fm.rep.split('#', 1)[0][:38]:38s} {fm.n:>4} {fm.reconstitution:>6.2f} "
            f"{fm.dominance:>5.2f} {fm.entropy:>5.2f} {str(fm.cross_superfamily):>5}  "
            f"{format_composition(fm.clade_composition)[:56]}")
    ax.text(0.01, 0.99, "\n".join(lines), va="top", ha="left", family="monospace", fontsize=7)
    pdf.savefig(fig)
    plt.close(fig)


def select_tp_families(all_metrics, n_tp):
    """Cleanest true positives for reference: dominance-gate passers ('safe'),
    highest reconstitution first."""
    tps = [fm for fm in all_metrics if fm.verdict == "safe"]
    tps.sort(key=lambda x: (-x.reconstitution, -x.dominance, -x.n))
    return tps[:max(0, n_tp)]


def render_pdf(path, metrics, families, internal_map, specs, thresholds,
               tp_pages=3, borderline_pages=10):
    fp = sorted((fm for fm in metrics if fm.verdict == "false_positive"),
                key=lambda x: (x.dominance, -x.n))
    borderline = select_borderline_families(
        metrics, thresholds["max_recon"], thresholds["max_dominance"], borderline_pages)
    tp = select_tp_families(metrics, tp_pages)
    pages = 0
    with PdfPages(path) as pdf:
        plot_summary_page(pdf, metrics, fp, thresholds)
        pages += 1
        for fm in fp:
            plot_family_page(pdf, fm, families[fm.rep], internal_map, specs,
                             note="FALSE POSITIVE - failed both gates")
            pages += 1
        for _gap, gate, fm in borderline:
            note = ("NEAR-FP - passed dominance gate (raise --max-dominance to flag)"
                    if gate == "dominance" else
                    "NEAR-FP - recon-rescued (raise --max-recon to flag)")
            plot_family_page(pdf, fm, families[fm.rep], internal_map, specs, note=note)
            pages += 1
        for fm in tp:
            plot_family_page(pdf, fm, families[fm.rep], internal_map, specs,
                             note="SAFE reference (passed dominance gate)")
            pages += 1
    return pages


# -----------------------------------------------------------------------------
# Genome masking (folded in from mask_fp.py)
# -----------------------------------------------------------------------------
BLAST_TASK = "dc-megablast"                       # sensitive to ~85% id copies
OUTFMT = "6 sseqid pident sstart send qcovhsp"    # only the columns the mask needs
WRAP = 60                                         # output FASTA line width


def die(msg):
    sys.exit(f"error: {msg}")


def run(cmd, **kw):
    """Print a command, then run it; abort on non-zero exit."""
    print("  " + " ".join(cmd), file=sys.stderr)
    if subprocess.run(cmd, **kw).returncode != 0:
        die(f"command failed: {' '.join(cmd)}")


def iter_fasta(path):
    """Yield (header, sequence) one record at a time."""
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].rstrip("\n"), []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def merge(intervals):
    """Merge overlapping/adjacent 0-based half-open intervals."""
    intervals.sort()
    out = []
    for s, e in intervals:
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def coord_of(element_id):
    """'chrom:start-end#Class/Super/Clade' -> 'chrom:start-end' (nest tokens omit label)."""
    return element_id.split("#", 1)[0]


def clean_output_path(in_path):
    """Depth-TSV clean-copy path: '<name>_ltr.tsv' -> '<name>_clean_ltr.tsv', else
    insert '.clean' before the extension. Kept in the input's directory."""
    d, base = os.path.split(in_path)
    if base.endswith("_ltr.tsv"):
        base = base[:-len("_ltr.tsv")] + "_clean_ltr.tsv"
    else:
        root, ext = os.path.splitext(base)
        base = root + ".clean" + ext
    return os.path.join(d, base)


_NEST_PREFIXES = ("nest-inner:", "nest-outer:")


def scrub_nest_field(nest, fp_coords):
    """Drop any nest-inner:/nest-outer: token pointing at a deleted FP coord (both
    directions). Returns (new_field, n_dropped); '.' when nothing remains."""
    if not nest or nest == ".":
        return nest, 0
    kept, dropped = [], 0
    for tok in nest.split(";"):
        coord = None
        for p in _NEST_PREFIXES:
            if tok.startswith(p):
                coord = tok[len(p):]
                break
        if coord is not None and coord in fp_coords:
            dropped += 1
        else:
            kept.append(tok)
    return (";".join(kept) if kept else "."), dropped


def clean_depth_tsv(in_path, out_path, fp_coords):
    """Write an FP-cleaned copy of a depth TSV: drop rows whose element coord is an
    FP, and scrub dangling nest tokens (both directions) from surviving rows. Header,
    comment, and short/malformed lines pass through verbatim. Unchanged data lines are
    written byte-for-byte. Returns (rows_removed, tokens_scrubbed)."""
    removed = scrubbed = 0
    with open(in_path) as fh, open(out_path, "w") as out:
        for line in fh:
            if line.startswith("#"):
                out.write(line)
                continue
            raw = line.rstrip("\n")
            f = raw.split("\t")
            if len(f) < 2:
                out.write(line)
                continue
            if coord_of(f[0]) in fp_coords:
                removed += 1
                continue
            new_nest, n = scrub_nest_field(f[-1], fp_coords)
            if n:
                scrubbed += n
                f[-1] = new_nest
                out.write("\t".join(f) + "\n")
            else:
                out.write(line)
    return removed, scrubbed


def fp_fraction(n_fp, n_total):
    """FP elements / total consensus-cluster members; 0.0 if no members."""
    return (n_fp / n_total) if n_total else 0.0


# --- BLAST worker state (set once per process via the Pool initializer) ---
_DB = _PIDENT = _QCOV = None


def _init(db, pident, qcov):
    global _DB, _PIDENT, _QCOV
    _DB, _PIDENT, _QCOV = db, pident, qcov


def parse_blast_hits(stdout, pident, qcov):
    """Parse outfmt6 'sseqid pident sstart send qcovhsp' -> passing 0-based half-open
    (sseqid, start, end) intervals (BLAST 1-based inclusive -> 0-based half-open)."""
    hits = []
    for line in stdout.splitlines():
        if not line.strip():
            continue
        sseqid, pid, ss, se, qc = line.split("\t")
        if float(pid) >= pident and float(qc) >= qcov:
            s, e = int(ss), int(se)
            if s > e:
                s, e = e, s
            hits.append((sseqid, s - 1, e))
    return hits


def _blast_one(record):
    """BLAST one query record against the worker DB; return passing intervals."""
    header, seq = record
    cmd = ["blastn", "-query", "-", "-db", _DB,
           "-task", BLAST_TASK, "-dust", "no", "-evalue", "1e-5",
           "-max_target_seqs", "100000", "-num_threads", "1", "-outfmt", OUTFMT]
    r = subprocess.run(cmd, input=f">{header}\n{seq}\n", text=True, capture_output=True)
    if r.returncode != 0:
        raise RuntimeError(f"blastn failed on {header.split()[0]}: {r.stderr.strip()}")
    return parse_blast_hits(r.stdout, _PIDENT, _QCOV)


def mask_intervals_in_seq(seq, intervals):
    """Set each 0-based half-open interval in seq to N. Returns (masked_seq, n_bp)."""
    ba = bytearray(seq, "ascii")
    n = len(ba)
    masked = 0
    for s, e in intervals:
        s, e = max(0, s), min(n, e)
        if e > s:
            ba[s:e] = b"N" * (e - s)
            masked += e - s
    return ba.decode("ascii"), masked


def write_masked_genome(genome_path, merged, out_path):
    """Stream genome FASTA; mask each sseqid's merged intervals to N (one record in
    memory at a time). Returns total bp masked."""
    masked = 0
    with open(out_path, "w") as out:
        for header, seq in iter_fasta(genome_path):
            sid = header.split()[0]
            if sid in merged:
                seq, n = mask_intervals_in_seq(seq, merged[sid])
                masked += n
            out.write(f">{header}\n")
            for i in range(0, len(seq), WRAP):
                out.write(seq[i:i + WRAP] + "\n")
    return masked


# mmseqs easy-cluster settings (baked in; only --threads is exposed).
MMSEQS_MIN_SEQ_ID = 0.9
MMSEQS_COV = 0.95
MMSEQS_COV_MODE = 0
MMSEQS_CLUSTER_MODE = 1
MMSEQS_MASK = 0
MMSEQS_SENS = 7.5


def run_mmseqs(fp_fasta, out_prefix, threads):
    """mmseqs easy-cluster to purge redundancy among the FP LTRs. Returns the
    representative-sequence FASTA path (<out_prefix>_rep_seq.fasta)."""
    if shutil.which("mmseqs") is None:
        die("mmseqs not found on PATH")
    tmp = tempfile.mkdtemp(prefix="mmseqs.")
    try:
        run(["mmseqs", "easy-cluster", fp_fasta, out_prefix, tmp,
             "--min-seq-id", str(MMSEQS_MIN_SEQ_ID),
             "-c", str(MMSEQS_COV),
             "--cov-mode", str(MMSEQS_COV_MODE),
             "--cluster-mode", str(MMSEQS_CLUSTER_MODE),
             "--mask", str(MMSEQS_MASK),
             "-s", str(MMSEQS_SENS),
             "--threads", str(threads)])
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    rep = out_prefix + "_rep_seq.fasta"
    if not os.path.isfile(rep):
        die(f"mmseqs did not produce {rep}")
    return rep


def mask_genome(genome, query, out, pident=80.0, qcov=90.0, threads=None):
    """dc-megablast the query FASTA against the genome and hard-mask hits to N.
    One chromosome in memory at a time. Returns total bp masked."""
    for tool in ("makeblastdb", "blastn"):
        if shutil.which(tool) is None:
            die(f"{tool} not found on PATH")
    threads = threads or os.cpu_count()
    tmp = tempfile.mkdtemp(prefix="mask_fp.")
    try:
        db = os.path.join(tmp, "db")
        run(["makeblastdb", "-in", genome, "-dbtype", "nucl", "-out", db],
            stdout=subprocess.DEVNULL)
        queries = list(iter_fasta(query))
        intervals = defaultdict(list)
        with mp.Pool(threads, initializer=_init, initargs=(db, pident, qcov)) as pool:
            for hits in pool.imap_unordered(_blast_one, queries, chunksize=1):
                for sseqid, s, e in hits:
                    intervals[sseqid].append((s, e))
        merged = {sid: merge(v) for sid, v in intervals.items()}
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return write_masked_genome(genome, merged, out)


# -----------------------------------------------------------------------------
# CLI + orchestration
# -----------------------------------------------------------------------------
def parse_args(argv=None):
    ap = argparse.ArgumentParser(
        description="Flag false-positive LTR-RT families, write FP-cleaned depth TSVs, "
                    "and (only if FPs exceed --fp-mask-threshold) mmseqs-dedup the FP "
                    "LTRs and hard-mask them out of --genome.")
    ap.add_argument("--consensus-cluster", required=True,
                    help="mmseqs consensus-LTR cluster TSV (rep<TAB>member); defines families")
    ap.add_argument("--internal-cluster", required=True,
                    help="mmseqs internal-region cluster TSV (rep<TAB>member); same element IDs")
    ap.add_argument("--ltr-fasta", required=True,
                    help="consensus LTR FASTA (*.consensus.fa); source of FP output sequences")
    ap.add_argument("-o", "--out-prefix", required=True, help="output path prefix")
    ap.add_argument("--domains-tsv", nargs="+", default=None,
                    help="synLTR depth TSV(s); enables the domain overlay AND is the source "
                         "for the FP-cleaned *_clean_ltr.tsv outputs")
    ap.add_argument("--results", help="Kmer2LTR results table (plot geometry when --domains-tsv absent)")
    ap.add_argument("--element-fasta", help="full-element FASTA (plot geometry when --domains-tsv absent)")
    ap.add_argument("--min-members", type=int, default=DEFAULT_MIN_MEMBERS,
                    help="size gate: only families with >= this many members are analyzed (default 10)")
    ap.add_argument("--max-dominance", type=float, default=DEFAULT_MAX_DOMINANCE,
                    help="gate 1 (dominance): top real clade covers >= this fraction -> kept (default 0.51)")
    ap.add_argument("--max-recon", type=float, default=DEFAULT_MAX_RECON,
                    help="gate 2 (recon rescue): survives if reconstitution > this, else FP (default 0.51)")
    ap.add_argument("--ignore-clades", default="unknown,mixture",
                    help="comma-separated clade/superfamily tokens treated as uninformative")
    ap.add_argument("--no-plot", action="store_true", help="skip the structure PDF")
    ap.add_argument("--borderline-pages", type=int, default=10,
                    help="pages for non-flagged families closest to the FP boundary (default 10; 0 disables)")
    ap.add_argument("--tp-pages", type=int, default=3,
                    help="clean true-positive reference pages (default 3; 0 disables)")
    ap.add_argument("--pdf-out", default=None, help="PDF path (default <prefix>.fp_structure.pdf)")
    # --- masking stage ---
    ap.add_argument("-g", "--genome", default=None,
                    help="genome FASTA to hard-mask (only used if FP fraction > --fp-mask-threshold)")
    ap.add_argument("--masked-out", default=None,
                    help="masked-genome output (default <genome_stem>_FP_masked.fa)")
    ap.add_argument("--fp-mask-threshold", type=float, default=0.10,
                    help="mask iff (FP elements / total members) > this fraction (default 0.10; "
                         "0 = always mask, >=1 = never)")
    ap.add_argument("--pident", type=float, default=80.0,
                    help="masking: min %% identity to mask an HSP (default 80)")
    ap.add_argument("--qcov", type=float, default=90.0,
                    help="masking: min %% of query covered, qcovhsp (default 90)")
    ap.add_argument("--threads", type=int, default=os.cpu_count(),
                    help="threads for mmseqs and BLAST (default all cores)")
    ap.add_argument("-v", "--verbose", action="store_true", help="per-family metric lines")
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    ignore = frozenset(s for s in args.ignore_clades.split(",") if s)

    # ---- Stage A: flag FP families (unchanged behavior) ----
    member2rep = parse_clusters(args.consensus_cluster)
    internal_raw = parse_clusters(args.internal_cluster)

    cons_ids, int_ids = set(member2rep), set(internal_raw)
    missing = cons_ids - int_ids
    if missing:
        frac = len(missing) / len(cons_ids)
        if frac > 0.01:
            sys.exit(f"[ERROR] {len(missing)} ({frac:.1%}) consensus members absent from "
                     f"the internal TSV; aborting (expected identical member sets)")
        print(f"[WARN] {len(missing)} members absent from internal TSV; treating as orphans",
              file=sys.stderr)
    internal_map = dict(internal_raw)
    for m in missing:
        internal_map[m] = f"__orphan__{m}"

    all_families = group_families(member2rep)
    families = {rep: mem for rep, mem in all_families.items() if len(mem) >= args.min_members}
    print(f"[INFO] {len(member2rep)} elements, {len(all_families)} families, "
          f"{len(families)} with >= {args.min_members} members", file=sys.stderr)

    metrics = []
    for rep, members in families.items():
        fmet = compute_family_metrics(rep, members, internal_map, ignore)
        fmet.verdict = classify_family(fmet, args.max_recon, args.max_dominance)
        metrics.append(fmet)
        if args.verbose:
            tag = {"false_positive": "FP ", "recovered": "rec", "safe": "   "}[fmet.verdict]
            print(f"[{tag}] {rep.split('#', 1)[0]}  n={fmet.n} recon={fmet.reconstitution:.2f} "
                  f"dom={fmet.dominance:.2f} ent={fmet.entropy:.2f}  -> {fmet.verdict}", file=sys.stderr)

    fp = [fmet for fmet in metrics if fmet.verdict == "false_positive"]
    n_safe = sum(1 for fmet in metrics if fmet.verdict == "safe")
    n_rec = sum(1 for fmet in metrics if fmet.verdict == "recovered")
    print(f"[INFO] gate 1 passed (safe): {n_safe}   gate 2 rescued (recovered): {n_rec}   "
          f"false-positive: {len(fp)} families / {sum(fmet.n for fmet in fp)} elements",
          file=sys.stderr)

    write_family_scores(args.out_prefix + ".family_scores.tsv", metrics)
    write_fp_families(args.out_prefix + ".fp_families.tsv", metrics)

    fp_members = [m for fmet in fp for m in families[fmet.rep]]
    n_written, n_missing = write_fp_fasta(args.out_prefix + ".fp_LTRs.fa", fp_members, args.ltr_fasta)
    print(f"[INFO] wrote {n_written} LTR sequences to {args.out_prefix}.fp_LTRs.fa "
          f"({n_missing} missing)", file=sys.stderr)

    if not args.no_plot:
        needed = {m for mem in families.values() for m in mem}
        if args.domains_tsv:
            depth = load_depth_tsvs(args.domains_tsv)
            n_missing_dom = len(needed - set(depth))
            if n_missing_dom:
                print(f"[WARN] {n_missing_dom} family members absent from the depth TSV(s); "
                      f"drawn without domains", file=sys.stderr)
            specs = build_domain_specs({eid: de for eid, de in depth.items() if eid in needed})
        else:
            if not args.results or not args.element_fasta:
                sys.exit("[ERROR] plotting needs either --domains-tsv, or both --results and "
                         "--element-fasta (or pass --no-plot)")
            results = load_results(args.results)
            lengths = load_lengths(args.element_fasta)
            specs = build_plain_specs(needed, results, lengths)
        pdf_out = args.pdf_out or (args.out_prefix + ".fp_structure.pdf")
        thresholds = {"min_members": args.min_members, "max_recon": args.max_recon,
                      "max_dominance": args.max_dominance}
        pages = render_pdf(pdf_out, metrics, families, internal_map, specs, thresholds,
                           args.tp_pages, args.borderline_pages)
        print(f"[INFO] wrote {pdf_out} ({pages} pages)", file=sys.stderr)

    # ---- Stage B: FP-cleaned depth TSVs (always, when --domains-tsv given) ----
    fp_coords = {coord_of(m) for m in fp_members}
    if args.domains_tsv:
        for in_path in args.domains_tsv:
            out_path = clean_output_path(in_path)
            removed, scrubbed = clean_depth_tsv(in_path, out_path, fp_coords)
            print(f"[INFO] {out_path}: removed {removed} FP rows, scrubbed {scrubbed} "
                  f"dangling nest tokens", file=sys.stderr)

    # ---- Stage C: conditional mmseqs + masking ----
    frac = fp_fraction(len(fp_members), len(member2rep))
    print(f"[INFO] FP fraction: {len(fp_members)}/{len(member2rep)} = {frac:.4f} "
          f"(threshold {args.fp_mask_threshold})", file=sys.stderr)
    if frac <= args.fp_mask_threshold:
        print("[INFO] FP fraction <= threshold; FPs are not pervasive. "
              "Skipping mmseqs + masking.", file=sys.stderr)
        return 0

    fp_fa = args.out_prefix + ".fp_LTRs.fa"
    if not args.genome:
        print(
            f"\n[ACTION REQUIRED] FP fraction {frac:.1%} exceeds --fp-mask-threshold "
            f"{args.fp_mask_threshold:.1%}, so masking is recommended, but no --genome was\n"
            f"given. Stages A+B are complete. To finish masking, rerun with --genome <genome.fa>,\n"
            f"or run manually:\n"
            f"  mmseqs easy-cluster {fp_fa} {args.out_prefix}_mmseqs <tmp> --min-seq-id 0.90 "
            f"-c 0.95 --cov-mode 0 --cluster-mode 1 --mask 0 -s 7.5 --threads {args.threads}\n"
            f"  mask via {args.out_prefix}_mmseqs_rep_seq.fasta\n",
            file=sys.stderr)
        return 2

    if not os.path.isfile(args.genome):
        die(f"no such genome file: {args.genome}")
    for tool in ("mmseqs", "makeblastdb", "blastn"):
        if shutil.which(tool) is None:
            die(f"{tool} not found on PATH")

    print(f"[INFO] FP fraction > threshold; FPs are pervasive -> mmseqs dedup + masking.",
          file=sys.stderr)
    rep = run_mmseqs(fp_fa, args.out_prefix + "_mmseqs", args.threads)
    masked_out = args.masked_out or (
        os.path.splitext(os.path.basename(args.genome))[0] + "_FP_masked.fa")
    bp = mask_genome(args.genome, rep, masked_out, args.pident, args.qcov, args.threads)
    print(f"[INFO] masked {bp:,} bp -> N; wrote {masked_out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
