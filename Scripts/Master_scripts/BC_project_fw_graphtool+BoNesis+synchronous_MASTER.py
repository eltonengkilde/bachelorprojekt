#!/opt/anaconda3/envs/bachelor_env/bin/python
import csv, os, re, sys, time, keyword, threading, graphlib
import numpy as np
from collections import Counter
import graph_tool.all as gt
import bonesis

BASE_DIR             = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUTPUT_DIR           = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")

class _Tee:
    def __init__(self, *files): self.files = files
    def write(self, data):
        for f in self.files: f.write(data)
    def flush(self):
        for f in self.files: f.flush()

class _AttAgg:
    """Per-gene aggregate over a set of attractor states (numpy-backed).
    _stars counts states where the gene has '*' (oscillating in a cycle attractor)."""
    __slots__ = ("_n", "_ones", "_stars")
    def __init__(self, n, ones, stars=None): self._n = n; self._ones = ones; self._stars = stars or {}
    def stable_on(self, g):  return g in self._ones and self._n > 0 and self._ones[g] == self._n
    def stable_off(self, g): return g in self._ones and self._n > 0 and self._ones[g] == 0 and self._stars.get(g, 0) == 0
    def any_on(self, g):     return g in self._ones and (self._ones[g] > 0 or self._stars.get(g, 0) > 0)
    def any_off(self, g):    return g in self._ones and self._ones[g] < self._n
    def count_on(self, g):   return self._ones.get(g, 0)
    def count_off(self, g):  return self._n - self._ones[g] if g in self._ones else 0
    def __bool__(self):      return self._n > 0
    def __len__(self):       return self._n

BONESIS_TIMEOUT         = 60     # seconds; BoNesis runs first, simulation fallback if stalled
MAX_SIM_STEPS           = 1000   # max synchronous simulation steps before non-convergence warning
SINK_RECOVERY_THRESHOLD = 10000

# NETWORK CONFIGURATION - change only this block to switch networks
NETWORK_FILE = "filtered_GT_normalized.csv"
#              other options:
#                filtered_networkL_normalized.csv

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

BOOLEAN_RESERVED = {"TRUE", "FALSE", "NOT", "AND", "OR"}
ACT_REL = "Activation / Induction / Causation / Result"
SUP_REL = "Repression / Inhibition / Negative Regulation"

# make sure all node names are valid Python identifiers before building rules
def clean_name(name):
    if not name or not isinstance(name, str): return "unknown"
    name = re.sub(r"_+", "_", re.sub(r"[^a-zA-Z0-9_]", "_", name)).strip("_")
    if not name: return "unknown"
    if name[0].isdigit(): return "n" + name
    if name.upper() in BOOLEAN_RESERVED or keyword.iskeyword(name): return "node_" + name
    return name

def _to_py(r): return r.replace("!", " not ").replace("|", " or ").replace("&", " and ")

# read the network CSV, filter by entity category, build activator and suppressor dicts
activators, suppressors, all_nodes, gene_category = {}, {}, set(), {}
_seen_act, _seen_sup = set(), set()

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
            _seen_act.add((s, t))
        else:
            suppressors.setdefault(t, set()).add(s)
            _seen_sup.add((s, t))

edges_act = sorted(_seen_act); edges_sup = sorted(_seen_sup)
del _seen_act, _seen_sup

print(f"Loaded {len(edges_act)+len(edges_sup)} unique edges  ({len(all_nodes)} nodes)")

# shorthand for entity type, then print how many nodes fall into each category
def cat(g): c = gene_category.get(g); return CAT_MAP.get(c, c or "?")
def mol_flag(g): return " !" if gene_category.get(g) not in MOLECULAR_CATS else ""

cat_counts = Counter(gene_category.values())
print(f"\nEntity categories ({len(cat_counts)}):")
for c, n in sorted(cat_counts.items(), key=lambda x: x[1], reverse=True): print(f"  {n:>6}  {c}")

# build the graph-tool graph, both act and sup edges go in as plain directed edges
t0 = time.perf_counter()
node_list = sorted(all_nodes)
node_idx  = {name: i for i, name in enumerate(node_list)}
g_full = gt.Graph(directed=True)
g_full.add_vertex(len(node_list))
for s, t in edges_act: g_full.add_edge(node_idx[s], node_idx[t])
for s, t in edges_sup: g_full.add_edge(node_idx[s], node_idx[t])
print(f"\ngraph-tool: {g_full.num_vertices()} vertices, {g_full.num_edges()} edges  ({time.perf_counter()-t0:.2f}s)")

# phase 1 - characterise the full network, find hubs with PageRank, check for feedback loops
print(f"\n{'='*70}")
print(f"PHASE 1  nodes: {len(all_nodes)}  |  edges: {len(edges_act)+len(edges_sup)} ({len(edges_act)} act / {len(edges_sup)} sup)")
t0 = time.perf_counter()
pr = gt.pagerank(g_full)
print(f"  PageRank ({time.perf_counter()-t0:.2f}s), top 10:")
for i in sorted(range(g_full.num_vertices()), key=lambda i: pr[i], reverse=True)[:10]:
    print(f"    {node_list[i]:40s}  PR: {pr[i]:.6f}  [{cat(node_list[i])}]")
_, hist = gt.label_components(g_full)
print(f"  SCCs: {hist.shape[0]}  |  non-trivial: {int((hist>1).sum())}  |  largest: {int(hist.max())} nodes")

t_total_start = time.perf_counter()

# ask which node to analyse and how many hops downstream
while True:
    _raw = input("\nSource node to investigate: ")
    source_gene = clean_name(_raw)          # apply same normalisation as CSV loading
    if source_gene not in node_idx:
        matches = sorted(n for n in all_nodes if _raw.upper() in n.upper())
        print(f"WARNING: '{_raw}' (normalised: '{source_gene}') not found. Try again.")
        if matches: print(f"  Closest matches: {matches[:5]}")
    else:
        print(f"OK: '{source_gene}'  [{cat(source_gene)}]"); break

while True:
    try:
        MAX_HOPS = int(input("\nNumber of hops to explore downstream: "))
        if MAX_HOPS < 1: print("Must be a whole number larger than 0")
        else: print(f"OK: {MAX_HOPS} hop(s)"); break
    except ValueError:
        print("Please enter a whole number.")

if input('\nType "yolo" to run: ') != "yolo":
    print("Exiting."); exit()

os.makedirs(OUTPUT_DIR, exist_ok=True)
_safe_gene = re.sub(r'[^a-zA-Z0-9_-]', '_', source_gene)
_outpath = os.path.join(OUTPUT_DIR, f"{_safe_gene}_{MAX_HOPS}hop_{time.strftime('%Y%m%d_%H%M%S')}.txt")
_logfile = open(_outpath, "w", encoding="utf-8")

try:
    sys.stdout = _Tee(sys.__stdout__, _logfile)
    print(f"Output -> {_outpath}")

    # BFS from source node, keep every node reachable within MAX_HOPS directed steps
    dist = gt.shortest_distance(g_full, source=g_full.vertex(node_idx[source_gene]), directed=True)
    subgraph_nodes = {node_list[i] for i in range(g_full.num_vertices()) if dist[i] <= MAX_HOPS}
    sub_node_list  = sorted(subgraph_nodes)
    sub_node_idx   = {name: i for i, name in enumerate(sub_node_list)}

    _sub_edges = sum(1 for s, t in edges_act + edges_sup if s in subgraph_nodes and t in subgraph_nodes)
    print(f"\n{'='*70}")
    print(f"PHASE 2  {source_gene}, {MAX_HOPS} hop(s)  |  {len(subgraph_nodes)} nodes, {_sub_edges} edges")
    for c, n in sorted(Counter(cat(n) for n in subgraph_nodes).items(), key=lambda x: x[1], reverse=True):
        print(f"    {n:>4}  {c}")

    # characterise the subgraph - SCCs
    print(f"\n{'='*70}")
    print(f"PHASE 2.5  Subgraph topology (graph-tool)")
    g_sub = gt.Graph(directed=True)
    g_sub.add_vertex(len(sub_node_list))
    for s, t in edges_act + edges_sup:
        if s in sub_node_idx and t in sub_node_idx:
            g_sub.add_edge(sub_node_idx[s], sub_node_idx[t])

    sub_comp, sub_hist = gt.label_components(g_sub)
    large_sub = int((sub_hist > 1).sum())
    print(f"  SCCs (Strongly Connected Components): {sub_hist.shape[0]}  |  non-trivial: {large_sub}")
    if large_sub:
        scc_members = [sub_node_list[i] for i in range(g_sub.num_vertices()) if int(sub_comp[i]) == int(sub_hist.argmax())]
        print(f"  Largest SCC ({int(sub_hist.max())} nodes):")
        for m in sorted(scc_members)[:10]: print(f"    {m:40s}  [{cat(m)}]{mol_flag(m)}")
        if len(scc_members) > 10: print(f"    ... and {len(scc_members)-10} more")

    # build one Boolean rule per regulated node: activators OR'd, suppressors AND NOT'd
    # rule form: (act1 | act2 | ...) & !sup1 & !sup2 & ...
    bn_dict = {}
    for target in set(activators) | set(suppressors):
        if target not in subgraph_nodes: continue
        act = " | ".join(s for s in activators.get(target, []) if s in subgraph_nodes)
        sup = " & ".join("!" + s for s in suppressors.get(target, []) if s in subgraph_nodes)
        if act and sup: bn_dict[target] = "(" + act + ") & " + sup
        elif act: bn_dict[target] = act
        elif sup: bn_dict[target] = sup

    # sinks are regulated but nothing else depends on them, separate them out before simulation
    # they get excluded from the BN and evaluated afterwards using their own rules
    referenced = set()
    for rule in bn_dict.values(): referenced |= set(re.findall(r'\b[a-zA-Z_]\w*\b', str(rule)))
    sink_nodes = {g for g in bn_dict if g not in referenced and g != source_gene}
    bn_dict_pruned = {g: f for g, f in bn_dict.items() if g not in sink_nodes}
    n_pruned = len(bn_dict_pruned)
    print(f"  Regulated: {len(bn_dict)}  |  pruned: {n_pruned}  |  sinks: {len(sink_nodes)}")

    # two starting conditions: all nodes OFF (dark) or all ON (permissive), source node locked
    bn_resting_dict   = dict(bn_dict_pruned); bn_resting_dict[source_gene]   = "0"
    bn_perturbed_dict = dict(bn_dict_pruned); bn_perturbed_dict[source_gene] = "1"

    # synchronous simulation - update all nodes at once each step until a state repeats
    def simulate(rules, start_state, locked, val, max_steps=MAX_SIM_STEPS):
        genes = sorted(start_state.keys())
        idx   = {g: i for i, g in enumerate(genes)}
        li    = idx[locked]
        cc    = [(idx[g], compile(_to_py(r).strip(), "<string>", "eval"))
                 for g, r in rules.items() if g != locked]
        state = [start_state.get(g, 0) for g in genes]; state[li] = val
        history, seen = [], {}
        for i in range(max_steps):
            key = tuple(state)
            if key in seen:
                return [{g: s[idx[g]] for g in genes} for s in history[seen[key]:]], i, True
            seen[key] = i; history.append(list(state))
            ns = {g: bool(state[i]) for g, i in idx.items()}
            nw = list(state); nw[li] = val
            for gi, cd in cc: nw[gi] = 1 if eval(cd, {"__builtins__": {}}, ns) else 0
            state = nw
        return [{g: state[idx[g]] for g in genes}], max_steps, False

    def _sim_label(states, conv):
        return ("fixed point" if len(states) == 1 else f"cycle/{len(states)}") + (" (conv)" if conv else " (max steps)")

    # try BoNesis first to get the full attractor landscape, fall back to simulation on timeout
    bn_resting   = bonesis.BooleanNetwork(bn_resting_dict)
    bn_perturbed = bonesis.BooleanNetwork(bn_perturbed_dict)
    dark_start = {g: 0 for g in bn_resting_dict}
    perm_start = {g: 1 for g in bn_resting_dict}
    print(f"  Attempting BoNesis ({n_pruned} nodes, timeout: {BONESIS_TIMEOUT}s)...")
    _bonesis_ok = False; t0 = time.perf_counter(); _box = []
    # daemon thread: if BoNesis stalls past BONESIS_TIMEOUT the thread is left running in background
    # it keeps running until the process exits, there is no clean way to stop a running solver
    _t = threading.Thread(daemon=True, target=lambda: _box.append((
        list(bn_resting.attractors(reachable_from=dark_start)),
        list(bn_perturbed.attractors(reachable_from=dark_start)),
        list(bn_resting.attractors(reachable_from=perm_start)),
        list(bn_perturbed.attractors(reachable_from=perm_start)),
    )))
    _t.start(); _t.join(BONESIS_TIMEOUT)
    if _box:
        dark_resting_att, dark_perturbed_att, perm_resting_att, perm_perturbed_att = _box[0]
        print(f"  BoNesis: {time.perf_counter()-t0:.2f}s  |  attractors: "
              f"dark-rest={len(dark_resting_att)}  dark-pert={len(dark_perturbed_att)}  "
              f"perm-rest={len(perm_resting_att)}  perm-pert={len(perm_perturbed_att)}", flush=True)
        _bonesis_ok = True
    else:
        print(f"  BoNesis timed out after {BONESIS_TIMEOUT}s, falling back to synchronous simulation")
        t0 = time.perf_counter()
        dr_states, _, dr_conv = simulate(bn_dict_pruned, dark_start, source_gene, 0)
        dp_states, _, dp_conv = simulate(bn_dict_pruned, dark_start, source_gene, 1)
        pr_states, _, pr_conv = simulate(bn_dict_pruned, perm_start, source_gene, 0)
        pp_states, _, pp_conv = simulate(bn_dict_pruned, perm_start, source_gene, 1)
        print(f"  dark  rest: {_sim_label(dr_states, dr_conv)}  pert: {_sim_label(dp_states, dp_conv)}")
        print(f"  perm  rest: {_sim_label(pr_states, pr_conv)}  pert: {_sim_label(pp_states, pp_conv)}")
        print(f"  Simulation: {time.perf_counter()-t0:.3f}s")
        if not all([dr_conv, dp_conv, pr_conv, pp_conv]):
            print(f"  WARNING: not all runs converged within {MAX_SIM_STEPS} steps")
        print(f"  NOTE: simulation finds ONE attractor per condition, use BoNesis at fewer hops for the full landscape")
        dark_resting_att, dark_perturbed_att, perm_resting_att, perm_perturbed_att = dr_states, dp_states, pr_states, pp_states

    using_simulation = not _bonesis_ok

    # sinks were excluded from the BN; recover their states using two strategies:
    # _recover_sinks_numpy: vectorised numpy for large attractor sets → returns _AttAgg
    # _recover_sinks_small: per-state eval with pre-compiled rules for small sets (KO, necessity)

    def _topo_sort_sinks(sink_rules):
        sink_set = set(sink_rules)
        deps = {g: set(re.findall(r'\b[a-zA-Z_]\w*\b', str(sink_rules[g]))) & sink_set
                for g in sink_set}
        return list(graphlib.TopologicalSorter(deps).static_order())

    def _recover_sinks_numpy(att_states, sink_order, src, src_val, extra_pins=None):
        """Evaluate all sink rules as boolean column ops on a uint8 matrix.
        Returns _AttAgg aggregate instead of individual states."""
        n = len(att_states)
        if n == 0:
            return _AttAgg(0, {})
        bn_nodes = list(att_states[0].keys())
        node_col = {g: i for i, g in enumerate(bn_nodes)}
        n_bn = len(bn_nodes)
        M = np.empty((n, n_bn + len(sink_order)), dtype=np.uint8)
        # '*' = oscillating gene in a cycle attractor; treat as 0 for rule evaluation
        raw = [[att[g] for g in bn_nodes] for att in att_states]
        M[:, :n_bn] = np.array([[0 if v == '*' else int(v) for v in row] for row in raw], dtype=np.uint8)
        S = np.array([[v == '*' for v in row] for row in raw], dtype=np.uint8)  # star mask
        if src in node_col:
            M[:, node_col[src]] = src_val
            S[:, node_col[src]] = 0
        if extra_pins:
            for g, v in extra_pins.items():
                if g in node_col:
                    M[:, node_col[g]] = v
                    S[:, node_col[g]] = 0
        for j, sink in enumerate(sink_order):
            col = n_bn + j
            node_col[sink] = col
            act_cols = [node_col[a] for a in activators.get(sink, ()) if a in node_col and node_col[a] != col]
            sup_cols = [node_col[s] for s in suppressors.get(sink, ()) if s in node_col and node_col[s] != col]
            if act_cols and sup_cols:
                M[:, col] = M[:, act_cols].any(axis=1) & ~M[:, sup_cols].any(axis=1)
            elif act_cols:
                M[:, col] = M[:, act_cols].any(axis=1)
            elif sup_cols:
                M[:, col] = ~M[:, sup_cols].any(axis=1)
            else:
                M[:, col] = 0
        ones  = {g: int(M[:, c].sum()) for g, c in node_col.items()}
        stars = {g: int(S[:, c].sum()) for g, c in node_col.items() if c < n_bn}
        return _AttAgg(n, ones, stars)

    def _recover_sinks_small(att_states, sink_compiled, src, src_val, extra_pins=None):
        """Per-state sink evaluation with pre-compiled rules; for small attractor sets."""
        _NB = {"__builtins__": {}}
        result = []
        for att in att_states:
            ns = {g: (False if v == '*' else bool(v)) for g, v in att.items()}
            ns[src] = bool(src_val)
            if extra_pins:
                for k, v in extra_pins.items(): ns[k] = bool(v)
            ext = dict(att); ext[src] = src_val
            if extra_pins: ext.update(extra_pins)
            for sink, cd in sink_compiled:
                try: val = 1 if eval(cd, _NB, ns) else 0
                except: val = 0
                ext[sink] = val; ns[sink] = bool(val)
            result.append(ext)
        return result

    sink_rules = {g: bn_dict[g] for g in sink_nodes}
    _n_sinks_total = len(sink_nodes)

    # Pre-filter: evaluate each sink's rule once with a representative resting vs perturbed
    # attractor state. Only recover sinks whose rule evaluates differently (i.e. actually change
    # state when source_gene is perturbed). Topological order handles sink-on-sink dependencies.
    _sink_pre_order = _topo_sort_sinks(sink_rules)
    _ns_r = {g: bool(v) for g, v in (dark_resting_att[0]   if dark_resting_att   else {}).items()}
    _ns_p = {g: bool(v) for g, v in (dark_perturbed_att[0] if dark_perturbed_att else {}).items()}
    _ns_r[source_gene] = False; _ns_p[source_gene] = True
    _NB = {"__builtins__": {}}
    _sink_relevant = set()
    for _g in _sink_pre_order:
        _expr = _to_py(sink_rules[_g])
        try:    _rv = eval(_expr, _NB, _ns_r)
        except: _rv = None
        try:    _pv = eval(_expr, _NB, _ns_p)
        except: _pv = None
        _ns_r[_g] = bool(_rv) if _rv is not None else False
        _ns_p[_g] = bool(_pv) if _pv is not None else False
        if _rv is None or _pv is None or _rv != _pv:
            _sink_relevant.add(_g)
    sink_nodes = _sink_relevant
    sink_rules = {g: sink_rules[g] for g in sink_nodes}
    print(f"  Sink pre-filter: {len(sink_nodes)} / {_n_sinks_total} sinks change state between resting and perturbed", flush=True)

    if len(sink_nodes) <= SINK_RECOVERY_THRESHOLD:
        _n_att_total = sum(len(a) for a in (dark_resting_att, dark_perturbed_att, perm_resting_att, perm_perturbed_att))
        print(f"  Sink recovery: {len(sink_nodes)} sinks × {_n_att_total} attractor states...", flush=True)
        sink_order    = _topo_sort_sinks(sink_rules)
        sink_compiled = [(g, compile(_to_py(sink_rules[g]).strip(), "<string>", "eval")) for g in sink_order]
    else:
        sink_order    = []
        sink_compiled = []
        print(f"  Sink recovery skipped ({len(sink_nodes)} sinks > {SINK_RECOVERY_THRESHOLD}), sinks excluded from classification")

    _t_sr = time.perf_counter()
    dark_resting_full   = _recover_sinks_numpy(dark_resting_att,   sink_order, source_gene, 0)
    dark_perturbed_full = _recover_sinks_numpy(dark_perturbed_att, sink_order, source_gene, 1)
    perm_resting_full   = _recover_sinks_numpy(perm_resting_att,   sink_order, source_gene, 0)
    perm_perturbed_full = _recover_sinks_numpy(perm_perturbed_att, sink_order, source_gene, 1)
    if sink_order:
        print(f"  Sink recovery: {time.perf_counter()-_t_sr:.3f}s", flush=True)

    # direct targets: nodes whose activator or suppressor set contains source_gene
    direct_targets = sorted(
        {t for t in bn_resting_dict if source_gene in activators.get(t, set()) or source_gene in suppressors.get(t, set())}
        | {t for t in sink_nodes    if source_gene in activators.get(t, set()) or source_gene in suppressors.get(t, set())}
    )

    # classify each node: activated means it goes from OFF to ON when source is perturbed
    def stable_on(g, atts):
        if isinstance(atts, _AttAgg): return atts.stable_on(g)
        return bool(atts) and all(a.get(g, 0) == 1 for a in atts)
    def stable_off(g, atts):
        if isinstance(atts, _AttAgg): return atts.stable_off(g)
        return bool(atts) and all(a.get(g, 0) == 0 for a in atts)
    def slabel(g, atts): return "ON" if stable_on(g, atts) else "OFF" if stable_off(g, atts) else "var"
    def tags(g): return "[" + ", ".join([cat(g)] + (["sink"] if g in sink_nodes else [])) + "]" + mol_flag(g)

    G = bn_dict.keys() | {source_gene}
    dark_activated  = sorted({g for g in G if stable_off(g, dark_resting_full)  and stable_on(g,  dark_perturbed_full)} - {source_gene})
    dark_suppressed = sorted({g for g in G if stable_on(g,  dark_resting_full)  and stable_off(g, dark_perturbed_full)} - {source_gene})
    perm_activated  = sorted({g for g in G if stable_off(g, perm_resting_full)  and stable_on(g,  perm_perturbed_full)} - {source_gene})
    perm_suppressed = sorted({g for g in G if stable_on(g,  perm_resting_full)  and stable_off(g, perm_perturbed_full)} - {source_gene})
    robust_activated  = sorted(set(dark_activated)  & set(perm_activated))
    robust_suppressed = sorted(set(dark_suppressed) & set(perm_suppressed))

    conditional_derepressed = sorted(
        (g, perm_resting_full.count_off(g))
        for g in G if g not in perm_activated and g != source_gene
        and perm_resting_full.any_off(g)
        and perm_perturbed_full.any_on(g)
    )
    conditional_suppressed = sorted(
        (g, perm_resting_full.count_on(g))
        for g in G if g not in perm_suppressed and g != source_gene
        and perm_resting_full.any_on(g)
        and perm_perturbed_full.any_off(g)
    )

    # knockout from the active fixed point
    # fixed point: run KO once from the single settled state
    # cycle: run KO from every cycle state; take intersection so only robust effects are reported
    ko_all_conv   = True
    ko_cycle_note = ""

    if len(dark_perturbed_att) > 1:
        ko_cycle_note = f" (cycle/{len(dark_perturbed_att)}: intersection across all {len(dark_perturbed_att)} starting states)"
    cand_maintained = {g for g in G if stable_on(g,  dark_perturbed_full) and g != source_gene}
    cand_suppressed = {g for g in G if stable_off(g, dark_perturbed_full) and g != source_gene}
    for cycle_state in dark_perturbed_att:
        if not cand_maintained and not cand_suppressed: break
        ko_att, _, _conv = simulate(bn_dict_pruned, cycle_state, source_gene, 0)
        if not _conv: ko_all_conv = False
        ko_full = _recover_sinks_small(ko_att, sink_compiled, source_gene, 0)
        cand_maintained = {g for g in cand_maintained if not stable_on(g, ko_full)}
        cand_suppressed = {g for g in cand_suppressed if     stable_on(g, ko_full)}
    ko_maintained = sorted(cand_maintained)
    ko_suppressed = sorted(cand_suppressed)

    # KO each direct target and check if any perm-activated node loses its stable ON state
    necessary, dispensable = [], []
    if perm_activated:
        # use the same solver as the attractor analysis to stay consistent
        # BoNesis checks all attractors, simulation follows one trajectory from all-ones
        nec_method = "BoNesis" if not using_simulation else "synchronous simulation"
        print(f"\n  Necessity test ({len(direct_targets)} direct target(s), method: {nec_method})...", flush=True)
        t0 = time.perf_counter()
        for candidate in direct_targets:
            if not using_simulation:
                ko = dict(bn_perturbed_dict); ko[candidate] = "0"
                _bn_ko   = bonesis.BooleanNetwork(ko)
                _nec_box = []
                _nt = threading.Thread(daemon=True, target=lambda _b=_bn_ko, _nb=_nec_box, _k=ko:
                    _nb.append(list(_b.attractors(reachable_from={g: 1 for g in _k}))))
                _nt.start(); _nt.join(BONESIS_TIMEOUT)
                if _nec_box:
                    ko_states = _nec_box[0]
                else:
                    _ko_fb = dict(bn_dict_pruned)
                    if candidate in _ko_fb: _ko_fb[candidate] = "0"
                    ko_states, _, _ = simulate(_ko_fb, perm_start, source_gene, 1)
                    print(f"    {candidate}: necessity BoNesis timed out, used simulation fallback", flush=True)
            else:
                ko = dict(bn_dict_pruned)
                if candidate in bn_dict_pruned:
                    ko[candidate] = "0"
                ko_states, _, _ = simulate(ko, perm_start, source_gene, 1)
            ko_full = _recover_sinks_small(ko_states, sink_compiled, source_gene, 1,
                                           extra_pins={candidate: 0})
            lost = sorted(g for g in perm_activated if g != candidate and not stable_on(g, ko_full))
            if lost: necessary.append((candidate, lost))
            else:    dispensable.append(candidate)
        print(f"  Necessity: {time.perf_counter()-t0:.3f}s", flush=True)
    else:
        nec_method = "N/A"
        print(f"\n  Necessity test skipped, no nodes stably activated in permissive baseline.")

    # print run metadata
    mode = "synchronous simulation" if using_simulation else "BoNesis (complete attractor landscape)"
    print(f"\n" + "="*70)
    print(f"  RUN METADATA")
    print(f"  Network          : {NETWORK_FILE}")
    print(f"  Node             : {source_gene}  [{cat(source_gene)}]")
    print(f"  Hops             : {MAX_HOPS}  (downstream)")
    print(f"  Attractor solver : {mode}  |  BoNesis timeout: {BONESIS_TIMEOUT}s")
    print(f"  Necessity solver : {nec_method}")
    print(f"  Update scheme    : synchronous (all nodes updated simultaneously per step)")
    print(f"  Boolean rules    : activators combined with OR; suppressors combined with AND NOT")
    print(f"  Subgraph         : {len(subgraph_nodes)} nodes  |  pruned BN: {n_pruned}  |  sinks: {len(sink_nodes)}")
    print(f"  Run time         : {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"  ! marks experimental perturbation nodes (mutant / process / phenotype)")
    print(f"    These are not endogenous gene products - interpret with caution.")
    print(f"\n" + "="*70)
    print(f"  {source_gene} -- {MAX_HOPS} hop(s)  |  mode: {mode}")

    print(f"\n  Direct targets: {len(direct_targets)}")
    for g in direct_targets: print(f"    {g:40s}  {tags(g)}")

    # print the results
    if (dark_activated == perm_activated and dark_suppressed == perm_suppressed
            and not conditional_derepressed and not conditional_suppressed):
        print(f"\n{'='*70}")
        print(f"  EXP A = EXP B, dark and permissive backgrounds produce identical effects")
        print(f"  (dark: {len(dark_resting_att)} rest/{len(dark_perturbed_att)} pert att  |  perm: {len(perm_resting_att)} rest/{len(perm_perturbed_att)} pert att)")
        print(f"\n  Activated (OFF->ON): {len(perm_activated)}")
        for g in perm_activated: print(f"    + {g:40s}  {tags(g)}")
        print(f"\n  Suppressed (ON->OFF): {len(perm_suppressed)}")
        for g in perm_suppressed: print(f"    - {g:40s}  {tags(g)}")
    else:
        print(f"\n{'='*70}")
        print(f"  EXP A (dark background)  |  {len(dark_resting_att)} resting / {len(dark_perturbed_att)} perturbed attractor(s)")
        print(f"\n  Activated (OFF->ON): {len(dark_activated)}")
        for g in dark_activated: print(f"    + {g:40s}  {tags(g)}")
        print(f"\n  Suppressed (ON->OFF): {len(dark_suppressed)}")
        for g in dark_suppressed: print(f"    - {g:40s}  {tags(g)}")
        print(f"\n{'='*70}")
        print(f"  EXP B (permissive background)  |  {len(perm_resting_att)} resting / {len(perm_perturbed_att)} perturbed attractor(s)")
        print(f"\n  Activated (OFF->ON): {len(perm_activated)}")
        for g in perm_activated: print(f"    + {g:40s}  {tags(g)}")
        print(f"\n  Suppressed (ON->OFF): {len(perm_suppressed)}")
        for g in perm_suppressed: print(f"    - {g:40s}  {tags(g)}")
        if conditional_derepressed or conditional_suppressed:
            print(f"\n  Conditional effects:")
            for g, count in conditional_derepressed: print(f"    + ({count}/{len(perm_resting_att)} resting OFF) {g:36s}  {tags(g)}")
            for g, count in conditional_suppressed:  print(f"    - ({count}/{len(perm_resting_att)} resting ON)  {g:36s}  {tags(g)}")

    print(f"\n{'='*70}")
    print(f"  ROBUST EFFECTS")
    print(f"\n  Robustly activated: {len(robust_activated)}")
    for g in robust_activated: print(f"    + {g:40s}  {tags(g)}")
    print(f"\n  Robustly suppressed: {len(robust_suppressed)}")
    for g in robust_suppressed: print(f"    - {g:40s}  {tags(g)}")

    print(f"\n{'='*70}")
    print(f"  KNOCKOUT FROM ACTIVE STATE  (source removed from settled network){ko_cycle_note}")
    if not ko_all_conv: print(f"  WARNING: KO simulation did not converge within {MAX_SIM_STEPS} steps")
    print(f"\n  Activated (ON in active, falls OFF after KO): {len(ko_maintained)}")
    for g in ko_maintained: print(f"    + {g:40s}  {tags(g)}")
    print(f"\n  Suppressed (OFF in active, turns ON after KO): {len(ko_suppressed)}")
    for g in ko_suppressed: print(f"    - {g:40s}  {tags(g)}")

    if perm_activated:
        print(f"\n{'='*70}")
        print(f"  NECESSITY TEST  (method: {nec_method} from permissive background)")
        print(f"  Definition: a direct target is necessary if knocking it out (fixing to OFF)")
        print(f"  causes at least one permissive-activated node to lose stable ON status.")
        print(f"    Necessary  ({len(necessary)}):")
        for g, lost in necessary:
            print(f"    ! {g:38s}  {tags(g)}")
            print(f"        loss of stable activation: {', '.join(lost)}")
        print(f"    Dispensable ({len(dispensable)}):")
        for g in dispensable: print(f"      {g:40s}  {tags(g)}")

    print(f"\n{'='*70}")
    print(f"  Sink nodes: {len(sink_nodes)} of {_n_sinks_total} change state when {source_gene} is perturbed")
    for g in sorted(sink_nodes):
        r = slabel(g, perm_resting_full)
        p = slabel(g, perm_perturbed_full)
        print(f"    {g:40s}  resting: {r:3s}  perturbed: {p:3s}  {tags(g)}")

    print(f"\n  Total: {time.perf_counter()-t_total_start:.2f}s")
    print(f"\n{'='*70}")

finally:
    sys.stdout = sys.__stdout__
    _logfile.close()

print(f"\nOutput saved to: {_outpath}")