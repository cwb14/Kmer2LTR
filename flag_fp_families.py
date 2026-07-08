#!/usr/bin/env python3
"""Flag false-positive LTR-RT families (chimeric non-LTR repeats) from Kmer2LTR
cluster tables with a two-gate filter, and render a structure PDF.

Each analyzed family (>= --min-members) runs a two-gate pipeline:

  Gate 1 - dominance: if one real clade dominates the family
    (dominance >= --max-dominance) the labels are coherent, so the family PASSES
    and is kept as a true positive ('safe').

  Gate 2 - reconstitution: a family that FAILS gate 1 gets a second chance. If its
    internal regions reconstitute the family (reconstitution > --max-recon) it is
    'recovered' (true positive). If the internals also fail to subcluster
    (reconstitution <= --max-recon) the family is a 'false_positive'.

A true LTR-RT family clears at least one gate: coherent labels, or reconstituting
internals. A chimeric false positive fails both. If a family is a false positive
all its members are; their LTR sequences (the actual problem repeats) are written out.

# Merge the synLTR LTR-RT files. 
awk '/^>/{printf "\n%s\n",$0;next}{gsub(/[^ACGTacgt]/,"");printf "%s",$0}END{print ""}' ./B73_LTR_depth*_ltr.fa > B73_all_ltr.fa
# Process with Kmer2LTR.
python Kmer2LTR/Kmer2LTR.py -i B73_all_ltr.fa -o B73_all_ltr --ltr-cluster --internal-fasta -p 200 --min-seq-id 0.75
# Flag false positives. 
python flag_fp_families.py --consensus-cluster B73_all_ltr.consensus_id0.75_cluster.tsv --internal-cluster B73_all_ltr.internal_id0.75_cluster.tsv --ltr-fasta B73_all_ltr.consensus.fa --domains-tsv B73_LTR_depth*_ltr.tsv -o B73_fpcheck
# I could filter those FPs from the results or, if there are a lot, I can use "B73_fpcheck.fp_LTRs.fa" to hardmask the genome, then re-run synLTR on the hardmasked genome. 
# In maize, there are 424 FP identified this way from a pool of 132931 (0.3%), so not a prolific issue and not worth tampering with. 
# In dog, this approach identifies 6284 FP from a pool of 8244 (76%), so here, FPs are a big issue (Makes sense due to dog LINE and SINE). We'd need to use "Basen_fpcheck.fp_LTRs.fa" to mask the dog genome, then reannotate LTR-RTs.
"""

import argparse
import math
import os
import re
import sys
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
# CLI + orchestration
# -----------------------------------------------------------------------------
def parse_args(argv=None):
    ap = argparse.ArgumentParser(
        description="Flag false-positive LTR-RT families (chimeric non-LTR repeats) from "
                    "Kmer2LTR cluster tables with a two-gate filter: gate 1 = dominance "
                    "(coherent labels pass and are kept); gate 2 = reconstitution (rescues "
                    "gate-1 failures whose internals subcluster; the rest are false positives).")
    ap.add_argument("--consensus-cluster", required=True,
                    help="mmseqs consensus-LTR cluster TSV (rep<TAB>member); defines families")
    ap.add_argument("--internal-cluster", required=True,
                    help="mmseqs internal-region cluster TSV (rep<TAB>member); same element IDs")
    ap.add_argument("--ltr-fasta", required=True,
                    help="consensus LTR FASTA (*.consensus.fa); source of output sequences")
    ap.add_argument("-o", "--out-prefix", required=True, help="output path prefix")
    ap.add_argument("--domains-tsv", nargs="+", default=None,
                    help="synLTR depth TSV(s), e.g. B73_LTR_depth*_ltr.tsv; enables the "
                         "protein-domain overlay and supplies plot geometry (no --results / "
                         "--element-fasta needed in this mode)")
    ap.add_argument("--results", help="Kmer2LTR results table (LTR length + K2P; plot geometry when --domains-tsv absent)")
    ap.add_argument("--element-fasta", help="full-element FASTA (element lengths via .fai; plot geometry when --domains-tsv absent)")
    ap.add_argument("--min-members", type=int, default=DEFAULT_MIN_MEMBERS,
                    help="size gate: only families with >= this many members are analyzed (default 10)")
    ap.add_argument("--max-dominance", type=float, default=DEFAULT_MAX_DOMINANCE,
                    help="gate 1 (dominance): a family whose top real clade covers >= this "
                         "fraction of members passes and is kept (default 0.51)")
    ap.add_argument("--max-recon", type=float, default=DEFAULT_MAX_RECON,
                    help="gate 2 (recon rescue): a gate-1 failure survives if reconstitution "
                         "> this, otherwise it is a false positive (default 0.51)")
    ap.add_argument("--ignore-clades", default="unknown,mixture",
                    help="comma-separated clade/superfamily tokens treated as uninformative")
    ap.add_argument("--no-plot", action="store_true", help="skip the structure PDF")
    ap.add_argument("--borderline-pages", type=int, default=10,
                    help="pages for non-flagged families closest to the FP boundary "
                         "(default 10; 0 disables)")
    ap.add_argument("--tp-pages", type=int, default=3,
                    help="clean true-positive reference pages (default 3; 0 disables)")
    ap.add_argument("--pdf-out", default=None, help="PDF path (default <prefix>.fp_structure.pdf)")
    ap.add_argument("-v", "--verbose", action="store_true", help="per-family metric lines")
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    ignore = frozenset(s for s in args.ignore_clades.split(",") if s)

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
        fm = compute_family_metrics(rep, members, internal_map, ignore)
        fm.verdict = classify_family(fm, args.max_recon, args.max_dominance)
        metrics.append(fm)
        if args.verbose:
            tag = {"false_positive": "FP ", "recovered": "rec", "safe": "   "}[fm.verdict]
            print(f"[{tag}] {rep.split('#', 1)[0]}  n={fm.n} recon={fm.reconstitution:.2f} "
                  f"dom={fm.dominance:.2f} ent={fm.entropy:.2f}  -> {fm.verdict}", file=sys.stderr)

    fp = [fm for fm in metrics if fm.verdict == "false_positive"]
    n_safe = sum(1 for fm in metrics if fm.verdict == "safe")
    n_rec = sum(1 for fm in metrics if fm.verdict == "recovered")
    print(f"[INFO] gate 1 passed (safe): {n_safe}   gate 2 rescued (recovered): {n_rec}   "
          f"false-positive: {len(fp)} families / {sum(fm.n for fm in fp)} elements",
          file=sys.stderr)

    write_family_scores(args.out_prefix + ".family_scores.tsv", metrics)
    write_fp_families(args.out_prefix + ".fp_families.tsv", metrics)

    fp_members = [m for fm in fp for m in families[fm.rep]]
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

    return 0


if __name__ == "__main__":
    sys.exit(main())
