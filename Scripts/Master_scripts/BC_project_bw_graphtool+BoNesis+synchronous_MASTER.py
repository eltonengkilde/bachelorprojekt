#!/opt/anaconda3/envs/bachelor_env/bin/python
import csv, os, re, sys, time, keyword
from collections import Counter
import graph_tool.all as gt

BASE_DIR   = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")

class _Tee:
    def __init__(self, *files): self.files = files
    def write(self, data):
        for f in self.files: f.write(data)
    def flush(self):
        for f in self.files: f.flush()

MAX_SIGN_ITER = 20   # max iterations for sign-propagation convergence

# ===========================================================================
# NETWORK CONFIGURATION - change only this block to switch networks
# ===========================================================================
NETWORK_FILE = "filtered_GT_normalized.csv"
#              other options:
#                filtered_GT.csv
#                filtered_networkL_normalized.csv
#                filtered_networkL.csv

# categories for filtered_GT_normalized.csv / filtered_GT.csv:
CATEGORIES_TO_KEEP = {
    "gene",
    "protein",
    "mutant",       # overexpression / transgenic lines (e.g. miR395c overexpression)
    "metabolite",
    "process",
    "phenotype",
}
CAT_MAP = {
    "gene":       "Gene",
    "protein":    "Protein",
    "mutant":     "Mutant",
    "metabolite": "Metabolite",
    "process":    "Process",
    "phenotype":  "Phenotype",
}
# gene beats protein beats metabolite beats mutant beats process beats phenotype
CAT_PRIORITY   = {"gene": 0, "protein": 1, "metabolite": 2, "mutant": 3, "process": 4, "phenotype": 5}
MOLECULAR_CATS = {"gene", "protein", "metabolite"}

# categories for filtered_networkL_normalized.csv / filtered_networkL.csv
# (replace the block above with this when switching):
#
# CATEGORIES_TO_KEEP = {
#     "Gene / Protein",
#     "Phenotype / Trait / Disease",
#     "Chemical / Metabolite / Cofactor / Ligand",
#     "Biological Process / Pathway / Function / Regulatory / Signaling Mechanism",
# }
# CAT_MAP = {
#     "Gene / Protein": "Gene",
#     "Phenotype / Trait / Disease": "Phenotype",
#     "Chemical / Metabolite / Cofactor / Ligand": "Metabolite",
#     "Biological Process / Pathway / Function / Regulatory / Signaling Mechanism": "Pathway",
# }
# ===========================================================================

BOOLEAN_RESERVED = {"TRUE", "FALSE", "NOT", "AND", "OR"}  # keeps node names from clashing with Python keywords
ACT_REL = "Activation / Induction / Causation / Result"
SUP_REL = "Repression / Inhibition / Negative Regulation"

# make sure all node names are valid identifiers before storing them
def clean_name(name):
    if not name or not isinstance(name, str): return "unknown"
    name = re.sub(r"_+", "_", re.sub(r"[^a-zA-Z0-9_]", "_", name)).strip("_")
    if not name: return "unknown"
    if name[0].isdigit(): return "n" + name
    if name.upper() in BOOLEAN_RESERVED or keyword.iskeyword(name): return "node_" + name
    return name

# read the network CSV, filter by category, build activator and suppressor dicts
# out_act and out_sup are the forward direction - what each node regulates
activators, suppressors = {}, {}
out_act, out_sup        = {}, {}
all_nodes, gene_category = set(), {}
_seen_act, _seen_sup    = set(), set()  # sets catch duplicate edges from the CSV

network_path = os.path.join(BASE_DIR, "networks_used_by_scripts", NETWORK_FILE)
print(f"Network: {NETWORK_FILE}")

with open(network_path, newline="", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        if row["source_category"] not in CATEGORIES_TO_KEEP: continue
        if row["target_category"] not in CATEGORIES_TO_KEEP: continue
        rel = row["relationship_category"]
        if rel not in (ACT_REL, SUP_REL): continue
        s, t = clean_name(row["source"]), clean_name(row["target"])
        all_nodes.update((s, t))
        for _node, _col in ((s, "source_category"), (t, "target_category")):
            _new_c = row[_col]
            if _node not in gene_category or CAT_PRIORITY.get(_new_c, 99) < CAT_PRIORITY.get(gene_category[_node], 99):
                gene_category[_node] = _new_c
        if rel == ACT_REL:
            activators.setdefault(t, set()).add(s)
            out_act.setdefault(s, set()).add(t)
            _seen_act.add((s, t))
        else:
            suppressors.setdefault(t, set()).add(s)
            out_sup.setdefault(s, set()).add(t)
            _seen_sup.add((s, t))

# convert deduplication sets to sorted lists for the graph builder
edges_act = sorted(_seen_act)
edges_sup = sorted(_seen_sup)
del _seen_act, _seen_sup

print(f"Loaded {len(edges_act)+len(edges_sup)} unique edges  ({len(all_nodes)} nodes)")

# shorthand for entity type, then print how many nodes fall into each category
def cat(g): c = gene_category.get(g); return CAT_MAP.get(c, c or "?")
# marks experimental perturbation nodes that are not endogenous gene products
def mol_flag(g): return " !" if gene_category.get(g) not in MOLECULAR_CATS else ""
cat_counts = Counter(gene_category.values())
print(f"\nEntity categories ({len(cat_counts)}):")
for c, n in sorted(cat_counts.items(), key=lambda x: x[1], reverse=True):
    print(f"  {n:>6}  {c}")

# build the graph-tool graph, both act and sup edges go in as plain directed edges
node_list = sorted(all_nodes)
node_idx  = {n: i for i, n in enumerate(node_list)}
t_total_start = time.perf_counter()

t0 = time.perf_counter()
g_full = gt.Graph(directed=True)
g_full.add_vertex(len(node_list))
for s, t in edges_act: g_full.add_edge(node_idx[s], node_idx[t])
for s, t in edges_sup: g_full.add_edge(node_idx[s], node_idx[t])
print(f"\ngraph-tool: {g_full.num_vertices()} vertices, {g_full.num_edges()} edges  "
      f"({time.perf_counter()-t0:.2f}s)")

# phase 1 - characterise the full network, find hubs with PageRank, check for feedback loops
print(f"\n{'='*70}")
print(f"PHASE 1  nodes: {g_full.num_vertices()}  |  edges: {g_full.num_edges()} "
      f"({len(edges_act)} act / {len(edges_sup)} sup)")
t0 = time.perf_counter()
pr_map = gt.pagerank(g_full)
print(f"  PageRank - structural centrality, unsigned graph ({time.perf_counter()-t0:.2f}s), top 10:")
for i in sorted(range(g_full.num_vertices()), key=lambda i: pr_map[i], reverse=True)[:10]:
    print(f"    {node_list[i]:40s}  PR: {pr_map[i]:.6f}  [{cat(node_list[i])}]")
_, hist = gt.label_components(g_full)
print(f"  SCCs: {hist.shape[0]}  |  non-trivial: {int((hist>1).sum())}  |  "
      f"largest: {int(hist.max())} nodes")

# ask which node to analyse and how many hops upstream
while True:
    _raw = input("\nTarget node for upstream regulatory analysis: ")
    target_gene = clean_name(_raw)          # apply same normalisation as CSV loading
    if target_gene not in node_idx:
        matches = sorted(n for n in all_nodes if _raw.upper() in n.upper())
        print(f"WARNING: '{_raw}' (normalised: '{target_gene}') not found. Try again.")
        if matches: print(f"  Closest matches: {matches[:5]}")
    else:
        print(f"OK: '{target_gene}'  [{cat(target_gene)}]"); break

while True:
    try:
        MAX_HOPS = int(input("\nNumber of hops to explore upstream: "))
        if MAX_HOPS < 1: print("Must be a whole number larger than 0")
        else: print(f"OK: {MAX_HOPS} hops upstream from '{target_gene}'"); break
    except ValueError:
        print("Please enter a whole number.")

if input('\nType "yolo" to run: ') != "yolo":
    print("Invalid input, exiting."); exit()

os.makedirs(OUTPUT_DIR, exist_ok=True)
_safe_gene = re.sub(r'[^a-zA-Z0-9_-]', '_', target_gene)
_outpath = os.path.join(OUTPUT_DIR,
    f"{_safe_gene}_{MAX_HOPS}hop_bw_{time.strftime('%Y%m%d_%H%M%S')}.txt")
_logfile = open(_outpath, "w", encoding="utf-8")

try:
    sys.stdout = _Tee(sys.__stdout__, _logfile)
    print(f"Output -> {_outpath}")

    # reverse all edges and BFS from the target - every node we reach is upstream
    # distance = how many hops to get back to the target; keep everything within MAX_HOPS
    # the BFS ignores edge signs, both act and sup count equally for distance
    g_rev = gt.GraphView(g_full, reversed=True)
    dist  = gt.shortest_distance(g_rev, source=g_full.vertex(node_idx[target_gene]),
                                 directed=True)
    hop_of = {node_list[i]: int(dist[i])
              for i in range(g_full.num_vertices()) if int(dist[i]) <= MAX_HOPS}

    # shell_at[h] = nodes whose shortest path back to the target is exactly h hops
    shell_at = {h: sorted(n for n, d in hop_of.items() if d == h and n != target_gene)
                for h in range(1, MAX_HOPS + 1)}

    total_upstream = sum(len(shell_at[h]) for h in range(1, MAX_HOPS + 1))

    print(f"\n{'='*70}")
    print(f"RUN METADATA")
    print(f"  Network file      : {NETWORK_FILE}")
    print(f"  Categories kept   : {', '.join(sorted(CATEGORIES_TO_KEEP))}")
    print(f"  Target node       : {target_gene}  [{cat(target_gene)}]")
    print(f"  Hops upstream     : {MAX_HOPS}")
    print(f"  Subgraph size     : {total_upstream + 1} nodes  "
          f"(target + {total_upstream} upstream: "
          f"{', '.join(f'hop {h}: {len(shell_at[h])}' for h in range(1, MAX_HOPS+1))})")
    print(f"  Method            : signed path tracing (iterative convergence over subgraph)")
    print(f"  Date/time         : {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"\n  METHOD NOTE")
    print(f"  The upstream subgraph is defined by BFS from the target on the reversed")
    print(f"  network; every node within MAX_HOPS is included.  Signs are then")
    print(f"  propagated iteratively over all edges within that subgraph until no node")
    print(f"  changes classification.  Activation preserves sign; suppression inverts")
    print(f"  it (double-negative = activation).  Lateral edges between nodes at the")
    print(f"  same hop distance are included.  Convergence is guaranteed because sign")
    print(f"  sets only grow.  Mixed means a node reaches the target through both")
    print(f"  activating and suppressing paths.")
    print(f"  ! marks experimental perturbation nodes (mutant / process / phenotype).")
    print(f"    These are not endogenous gene products - interpret with caution.")

    # direct regulators are just nodes with a 1-hop edge straight into the target
    direct_act = sorted(g for g in activators.get(target_gene, set()) if g in hop_of)
    direct_sup = sorted(g for g in suppressors.get(target_gene, set()) if g in hop_of)
    print(f"\n{'='*70}")
    print(f"DIRECT REGULATORS of '{target_gene}'  (hop 1, structural edges)")
    print(f"  Activators  ({len(direct_act)}): {', '.join(direct_act) or '(none)'}")
    print(f"  Suppressors ({len(direct_sup)}): {', '.join(direct_sup) or '(none)'}")

    # characterise the upstream subgraph topology before sign propagation
    # non-trivial SCCs mean feedback loops are present; sign propagation may produce mixed results
    subgraph_nodes = set(hop_of.keys())
    _sub_node_list = sorted(subgraph_nodes)
    _sub_node_idx  = {name: i for i, name in enumerate(_sub_node_list)}
    g_sub = gt.Graph(directed=True)
    g_sub.add_vertex(len(_sub_node_list))
    for s, t in edges_act + edges_sup:
        if s in _sub_node_idx and t in _sub_node_idx:
            g_sub.add_edge(_sub_node_idx[s], _sub_node_idx[t])
    _sub_comp, _sub_hist = gt.label_components(g_sub)
    _large_sub = int((_sub_hist > 1).sum())
    print(f"\n{'='*70}")
    print(f"UPSTREAM SUBGRAPH TOPOLOGY")
    print(f"  SCCs: {_sub_hist.shape[0]}  |  non-trivial: {_large_sub}")
    if _large_sub:
        _scc_members = [_sub_node_list[i] for i in range(g_sub.num_vertices())
                        if int(_sub_comp[i]) == int(_sub_hist.argmax())]
        print(f"  Largest SCC ({int(_sub_hist.max())} nodes) — feedback loops present; "
              f"sign propagation may yield mixed results for these nodes:")
        for m in sorted(_scc_members)[:10]: print(f"    {m:40s}  [{cat(m)}]{mol_flag(m)}")
        if len(_scc_members) > 10: print(f"    ... and {len(_scc_members)-10} more")

    # propagate signs iteratively over all edges within the subgraph until no node changes
    # lateral edges between same-hop nodes are included on every iteration
    # activation keeps the sign, suppression flips it (double negative = activating)
    # frozenset({+1}) = activator, frozenset({-1}) = suppressor, {+1,-1} = mixed, empty = no path found
    # convergence is guaranteed because sign sets only grow (monotone lattice)
    net_sign       = {target_gene: frozenset({+1})}
    shell_conns    = {}   # shell_conns[N] = [(downstream node, '+' or '-'), ...]
    for _iter in range(MAX_SIGN_ITER):
        _changed = False
        for h in range(1, MAX_HOPS + 1):
            for node in shell_at[h]:
                signs = set()
                conns = []
                for dn in out_act.get(node, set()):
                    if dn in subgraph_nodes and net_sign.get(dn):
                        signs.update(net_sign[dn])               # activation keeps the sign
                        conns.append((dn, '+'))
                for dn in out_sup.get(node, set()):
                    if dn in subgraph_nodes and net_sign.get(dn):
                        signs.update({-s for s in net_sign[dn]}) # suppression flips it
                        conns.append((dn, '-'))
                new_sign = frozenset(signs)
                if new_sign != net_sign.get(node, frozenset()):
                    _changed = True
                net_sign[node]    = new_sign
                shell_conns[node] = sorted(conns, key=lambda x: (hop_of.get(x[0], 999), x[0]))
        if not _changed:
            print(f"  Sign propagation : converged in {_iter + 1} iteration(s)")
            break
    else:
        print(f"  WARNING: sign propagation did not converge after {MAX_SIGN_ITER} iterations")

    def classify(node):
        s = net_sign.get(node, frozenset())
        if   s == frozenset({+1}):  return "activator"
        elif s == frozenset({-1}):  return "suppressor"
        elif {+1, -1} <= s:         return "mixed"
        else:                       return "no_path"

    # show which downstream nodes each regulator connects through and what the net effect is
    def fmt_conns(node, max_show=4):
        conns = shell_conns.get(node, [])
        if not conns: return ""
        parts = []
        for dn, sign in conns[:max_show]:
            dn_signs = net_sign.get(dn, frozenset())
            net_via  = frozenset({-s for s in dn_signs}) if sign == '-' else dn_signs
            if   net_via == frozenset({+1}): net_lbl = "net +"
            elif net_via == frozenset({-1}): net_lbl = "net -"
            elif {+1, -1} <= net_via:        net_lbl = "net ~"
            else:                            net_lbl = "net ?"
            arrow = "->" if sign == '+' else "-|"
            parts.append(f"{arrow} {dn} [{net_lbl}]")
        if len(conns) > max_show:
            parts.append(f"(+{len(conns)-max_show} more)")
        return "  [via: " + ", ".join(parts) + "]"

    # print each hop shell - activators first, then suppressors, then mixed
    print(f"\n{'='*70}")
    print(f"UPSTREAM REGULATORY ANALYSIS  target: '{target_gene}'  |  {MAX_HOPS} hop(s)")
    print(f"  + = net activator  (all shortest paths are activating)")
    print(f"  - = net suppressor (all shortest paths are suppressing)")
    print(f"  ~ = mixed effect   (competing activating and suppressing paths)")

    summary_act = {}   # gene -> hop
    summary_sup = {}
    summary_mix = {}

    for h in range(1, MAX_HOPS + 1):
        shell  = shell_at[h]
        acts   = sorted(n for n in shell if classify(n) == "activator")
        sups   = sorted(n for n in shell if classify(n) == "suppressor")
        mixed  = sorted(n for n in shell if classify(n) == "mixed")
        nopath = sorted(n for n in shell if classify(n) == "no_path")

        print(f"\n{'='*70}")
        print(f"  HOP {h} SHELL  |  {len(shell)} new upstream nodes")

        print(f"\n  Net activators ({len(acts)}):")
        for g in acts:
            print(f"    + {g:50s}  [{cat(g)}]{mol_flag(g)}{fmt_conns(g)}")
            summary_act[g] = h
        if not acts: print("    (none)")

        print(f"\n  Net suppressors ({len(sups)}):")
        for g in sups:
            print(f"    - {g:50s}  [{cat(g)}]{mol_flag(g)}{fmt_conns(g)}")
            summary_sup[g] = h
        if not sups: print("    (none)")

        if mixed:
            print(f"\n  Mixed effect ({len(mixed)}) - competing paths with opposite signs:")
            for g in mixed:
                print(f"    ~ {g:50s}  [{cat(g)}]{mol_flag(g)}{fmt_conns(g)}")
                summary_mix[g] = h

        # no_path means the node was reachable by the unsigned BFS but has no incoming
        # signed paths within the subgraph - should be rare after iterative propagation
        if nopath:
            print(f"\n  No signed path found ({len(nopath)}):")
            for g in nopath[:10]:
                print(f"    ? {g:50s}  [{cat(g)}]{mol_flag(g)}")
            if len(nopath) > 10: print(f"    ... and {len(nopath)-10} more")

    # summary across all hops
    print(f"\n{'='*70}")
    print(f"SUMMARY  upstream regulatory analysis of '{target_gene}'  "
          f"({MAX_HOPS} hops)  |  network: {NETWORK_FILE}")
    print(f"  Date/time: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    print(f"\n  Net activators ({len(summary_act)}):")
    for g, h in sorted(summary_act.items(), key=lambda x: (x[1], x[0])):
        print(f"    + {g:50s}  [{cat(g)}]{mol_flag(g)}  hop {h}")
    if not summary_act: print("    (none)")

    print(f"\n  Net suppressors ({len(summary_sup)}):")
    for g, h in sorted(summary_sup.items(), key=lambda x: (x[1], x[0])):
        print(f"    - {g:50s}  [{cat(g)}]{mol_flag(g)}  hop {h}")
    if not summary_sup: print("    (none)")

    if summary_mix:
        print(f"\n  Mixed effect ({len(summary_mix)}):")
        for g, h in sorted(summary_mix.items(), key=lambda x: (x[1], x[0])):
            print(f"    ~ {g:50s}  [{cat(g)}]{mol_flag(g)}  hop {h}")

    print(f"\n  Total analysis time: {time.perf_counter() - t_total_start:.2f}s")
    print(f"\n{'='*70}")

finally:
    sys.stdout = sys.__stdout__
    _logfile.close()

print(f"\nOutput saved to: {_outpath}")
