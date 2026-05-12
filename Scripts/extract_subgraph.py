#!/opt/anaconda3/envs/bachelor_env/bin/python
"""
extract_subgraph.py
===================
Extracts forward (downstream) AND backward (upstream) subgraphs for a
given gene across all hop depths from MIN_HOPS to MAX_HOPS in one run.

Output files (all saved to OUT_DIR):
  subgraph_<GENE>_forward_nodes.csv   — all forward nodes with Hop_Distance 1-MAX_HOPS
  subgraph_<GENE>_forward_edges.csv   — all forward edges (within MAX_HOPS subgraph)
  subgraph_<GENE>_backward_nodes.csv  — all backward nodes with Hop_Distance 1-MAX_HOPS
  subgraph_<GENE>_backward_edges.csv  — all backward edges (within MAX_HOPS subgraph)
  subgraph_<GENE>_growth_summary.csv  — one row per (direction, hop), showing node/edge
                                        counts at each hop level for both directions

In Excel:
  - Filter nodes by Hop_Distance ≤ N to get the hop-N subgraph.
  - Use the growth_summary CSV to plot how the subgraph grows with each hop.
  - Compare forward vs backward in the same sheet using Direction column.
"""

import csv, os, re, keyword
from collections import deque

# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURE HERE
# ══════════════════════════════════════════════════════════════════════════════
GENE     = "MYB46"
MIN_HOPS = 1
MAX_HOPS = 12

BASE     = "/Users/elton_1stboot/Documents/Bachelor Git/bachelorprojekt"
NET_FILE = BASE + "/networks_used_by_scripts/filtered_networkL_normalized.csv"
OUT_DIR  = BASE
# ══════════════════════════════════════════════════════════════════════════════

CATEGORIES_TO_KEEP = {
    "Gene / Protein",
    "Phenotype / Trait / Disease",
    "Chemical / Metabolite / Cofactor / Ligand",
    "Biological Process / Pathway / Function / Regulatory / Signaling Mechanism",
}
ACT_REL = "Activation / Induction / Causation / Result"
SUP_REL = "Repression / Inhibition / Negative Regulation"
CAT_MAP = {
    "Gene / Protein": "Gene",
    "Phenotype / Trait / Disease": "Phenotype",
    "Chemical / Metabolite / Cofactor / Ligand": "Metabolite",
    "Biological Process / Pathway / Function / Regulatory / Signaling Mechanism": "Pathway",
}
BOOLEAN_RESERVED = {"TRUE", "FALSE", "NOT", "AND", "OR"}


def clean_name(name):
    if not name or not isinstance(name, str): return "unknown"
    name = re.sub(r"_+", "_", re.sub(r"[^a-zA-Z0-9_]", "_", name)).strip("_")
    if not name: return "unknown"
    if name[0].isdigit(): return "g" + name
    if name.upper() in BOOLEAN_RESERVED or keyword.iskeyword(name): return "gene_" + name
    return name


def short_cat(full_cat):
    return CAT_MAP.get(full_cat, full_cat or "?")


# ── Load network ──────────────────────────────────────────────────────────────
print(f"Loading network...")
edges    = []
adj_fwd  = {}
adj_bwd  = {}
node_cat = {}
node_orig = {}

with open(NET_FILE, newline="", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        if row["source_category"] not in CATEGORIES_TO_KEEP: continue
        if row["target_category"] not in CATEGORIES_TO_KEEP: continue
        rel = row["relationship_category"]
        if rel not in (ACT_REL, SUP_REL): continue
        s = clean_name(row["source"].strip())
        t = clean_name(row["target"].strip())
        sc, tc = row["source_category"], row["target_category"]
        node_cat[s] = sc; node_cat[t] = tc
        node_orig.setdefault(s, row["source"].strip())
        node_orig.setdefault(t, row["target"].strip())
        edges.append((s, t, rel, sc, tc))
        adj_fwd.setdefault(s, []).append(t)
        adj_bwd.setdefault(t, []).append(s)

print(f"  {len(node_cat):,} nodes  |  {len(edges):,} edges")

# Validate gene
gene_clean = clean_name(GENE)
if gene_clean not in node_cat:
    match = next((n for n in node_cat if n.upper() == gene_clean.upper()), None)
    if match:
        print(f"  Name corrected: '{GENE}' → '{match}'")
        gene_clean = match
    else:
        print(f"ERROR: '{GENE}' not found. Exiting."); exit(1)
print(f"  Gene: '{gene_clean}'  [{short_cat(node_cat[gene_clean])}]")


# ── BFS helper ────────────────────────────────────────────────────────────────
def bfs(start, adj, max_hops):
    """Return dict {node: hop_distance} for all nodes reachable within max_hops."""
    dist  = {start: 0}
    queue = deque([start])
    while queue:
        node = queue.popleft()
        if dist[node] >= max_hops:
            continue
        for nbr in adj.get(node, []):
            if nbr not in dist:
                dist[nbr] = dist[node] + 1
                queue.append(nbr)
    return dist


# ── Run BFS for both directions ───────────────────────────────────────────────
print(f"\nRunning BFS (hops 1–{MAX_HOPS}) for both directions...")
hop_fwd = bfs(gene_clean, adj_fwd, MAX_HOPS)   # forward  (downstream)
hop_bwd = bfs(gene_clean, adj_bwd, MAX_HOPS)   # backward (upstream)
print(f"  Forward : {len(hop_fwd):,} nodes reachable within {MAX_HOPS} hops")
print(f"  Backward: {len(hop_bwd):,} nodes reachable within {MAX_HOPS} hops")


# ── Write node CSVs ───────────────────────────────────────────────────────────
def write_nodes(hop_dict, direction, filename):
    with open(filename, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["Direction", "Gene_Cleaned", "Original_Name",
                    "Category", "Category_Short", "Hop_Distance", "Is_Centre_Gene"])
        for node in sorted(hop_dict, key=lambda n: (hop_dict[n], n)):
            w.writerow([
                direction,
                node,
                node_orig.get(node, node),
                node_cat.get(node, ""),
                short_cat(node_cat.get(node, "")),
                hop_dict[node],
                "Yes" if node == gene_clean else "No",
            ])
    print(f"  {os.path.basename(filename)}  ({len(hop_dict):,} rows)")


def write_edges(hop_dict, direction, filename):
    seen = set()
    rows = []
    for s, t, rel, sc, tc in edges:
        if s not in hop_dict or t not in hop_dict: continue
        key = (s, t, rel)
        if key in seen: continue
        seen.add(key)
        rows.append([
            direction, s, t,
            "Activation" if rel == ACT_REL else "Repression",
            short_cat(sc), short_cat(tc),
            hop_dict[s], hop_dict[t],
        ])
    with open(filename, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["Direction", "Source", "Target", "Relationship",
                    "Source_Category", "Target_Category",
                    "Source_Hop", "Target_Hop"])
        for row in rows:
            w.writerow(row)
    print(f"  {os.path.basename(filename)}  ({len(rows):,} rows)")


print("\nWriting node and edge CSVs...")
write_nodes(hop_fwd, "Forward",  os.path.join(OUT_DIR, f"subgraph_{gene_clean}_forward_nodes.csv"))
write_nodes(hop_bwd, "Backward", os.path.join(OUT_DIR, f"subgraph_{gene_clean}_backward_nodes.csv"))
write_edges(hop_fwd, "Forward",  os.path.join(OUT_DIR, f"subgraph_{gene_clean}_forward_edges.csv"))
write_edges(hop_bwd, "Backward", os.path.join(OUT_DIR, f"subgraph_{gene_clean}_backward_edges.csv"))


# ── Growth summary CSV ────────────────────────────────────────────────────────
# One row per (direction, hop): cumulative counts of nodes and edges at each level.
print("\nBuilding growth summary...")

def edge_counts_at_hop(hop_dict, max_hop):
    nodes_h = {n for n, d in hop_dict.items() if d <= max_hop}
    act = rep = 0
    seen = set()
    for s, t, rel, _, _ in edges:
        if s not in nodes_h or t not in nodes_h: continue
        key = (s, t, rel)
        if key in seen: continue
        seen.add(key)
        if rel == ACT_REL: act += 1
        else:              rep += 1
    return act, rep

summary_file = os.path.join(OUT_DIR, f"subgraph_{gene_clean}_growth_summary.csv")
with open(summary_file, "w", newline="", encoding="utf-8-sig") as f:
    w = csv.writer(f)
    w.writerow(["Direction", "Hop",
                "Cumulative_Nodes", "Cumulative_Edges",
                "Gene_Nodes", "Pathway_Nodes", "Phenotype_Nodes", "Metabolite_Nodes",
                "Activation_Edges", "Repression_Edges",
                "Pct_of_Full_Network_Nodes", "Pct_of_Full_Network_Edges"])

    total_nodes = len(node_cat)
    total_edges = len(set((s, t, rel) for s, t, rel, _, _ in edges))

    for direction, hop_dict in [("Forward", hop_fwd), ("Backward", hop_bwd)]:
        for hop in range(MIN_HOPS, MAX_HOPS + 1):
            nodes_h = {n for n, d in hop_dict.items() if d <= hop}
            cats = {}
            for n in nodes_h:
                c = short_cat(node_cat.get(n, ""))
                cats[c] = cats.get(c, 0) + 1
            act, rep = edge_counts_at_hop(hop_dict, hop)
            total_e = act + rep
            w.writerow([
                direction, hop,
                len(nodes_h), total_e,
                cats.get("Gene", 0),
                cats.get("Pathway", 0),
                cats.get("Phenotype", 0),
                cats.get("Metabolite", 0),
                act, rep,
                round(100 * len(nodes_h) / total_nodes, 2),
                round(100 * total_e     / total_edges,  2),
            ])
            print(f"  [{direction[:2]}] hop {hop:2d}  nodes={len(nodes_h):6,}  edges={total_e:6,}")

print(f"\n  {os.path.basename(summary_file)}  written")

# ── Quick comparison summary ──────────────────────────────────────────────────
print(f"\n{'='*65}")
print(f"  Full network:  {total_nodes:,} nodes  |  {total_edges:,} edges")
print(f"  Subgraph at hop {MAX_HOPS}:")
print(f"    Forward  — {len(hop_fwd):,} nodes  "
      f"({100*len(hop_fwd)/total_nodes:.1f}% of network)")
print(f"    Backward — {len(hop_bwd):,} nodes  "
      f"({100*len(hop_bwd)/total_nodes:.1f}% of network)")
print(f"\n  Filter nodes CSV by Hop_Distance ≤ N to get the hop-N subgraph.")
print(f"  Plot Cumulative_Nodes vs Hop in growth_summary for both directions.")
print(f"{'='*65}")
