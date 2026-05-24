#!/opt/anaconda3/envs/bachelor_env/bin/python
import csv, os, re, time, keyword, contextlib, threading, graphlib
import numpy as np
from collections import Counter
import graph_tool.all as gt
import bonesis

GENE                    = "MYB46"
MAX_TOTAL_SECONDS       = 1800   # 0.5-hour budget; a new hop only starts if time remains
MAX_HOPS                =  30     # safety cap — time limit is the primary stop condition
MAX_SIM_STEPS           = 1000
SINK_RECOVERY_THRESHOLD = 10000
SCRIPT_NAME = os.path.basename(__file__).replace('.py', '')
OUT_DIR     = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
BASE_DIR    = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.makedirs(OUT_DIR, exist_ok=True)

NETWORK_FILE = "filtered_GT_normalized.csv"

CATEGORIES_TO_KEEP = {"gene", "protein", "mutant", "metabolite", "process", "phenotype"}
BOOLEAN_RESERVED   = {"TRUE", "FALSE", "NOT", "AND", "OR"}
ACT_REL = "Activation / Induction / Causation / Result"
SUP_REL = "Repression / Inhibition / Negative Regulation"
CAT_MAP = {
    "gene":       "Gene",
    "protein":    "Protein",
    "mutant":     "Mutant",
    "metabolite": "Metabolite",
    "process":    "Process",
    "phenotype":  "Phenotype",
}
CAT_PRIORITY   = {"gene": 0, "protein": 1, "metabolite": 2, "mutant": 3, "process": 4, "phenotype": 5}
MOLECULAR_CATS = {"gene", "protein", "metabolite"}


class _AttAgg:
    """Per-gene aggregate over attractor states. _stars counts '*' (oscillating) states."""
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


def clean_name(name):
    if not name or not isinstance(name, str): return "unknown"
    name = re.sub(r"_+", "_", re.sub(r"[^a-zA-Z0-9_]", "_", name)).strip("_")
    if not name: return "unknown"
    if name[0].isdigit(): return "n" + name
    if name.upper() in BOOLEAN_RESERVED or keyword.iskeyword(name): return "node_" + name
    return name

def _to_py(r): return r.replace("!", " not ").replace("|", " or ").replace("&", " and ")

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

def eval_rule_simple(rule, state):
    try: return 1 if eval(_to_py(rule), {"__builtins__": {}},
                          {g: (False if v == '*' else bool(v)) for g, v in state.items()}) else 0
    except: return 0

def _topo_sort_sinks(sink_rules):
    sink_set = set(sink_rules)
    deps = {g: set(re.findall(r'\b[a-zA-Z_]\w*\b', str(sink_rules[g]))) & sink_set
            for g in sink_set}
    return list(graphlib.TopologicalSorter(deps).static_order())

def recover_sinks(attractors, sink_rules, src, src_val, sink_order, extra_pins=None):
    """Evaluate sink rules in topological order; extra_pins forces additional node values."""
    result = []
    for att in attractors:
        ext = {**att, src: src_val}
        if extra_pins:
            ext.update(extra_pins)
        for sink in sink_order:
            ext[sink] = eval_rule_simple(sink_rules[sink], ext)
        result.append(ext)
    return result

def _recover_sinks_numpy(att_states, sink_order, src, src_val, extra_pins=None):
    """Vectorised numpy sink evaluation over all attractor states. Returns _AttAgg.
    Handles '*' (oscillating gene in a cycle attractor): treats as 0 in rule evaluation,
    tracks star counts separately so stable_off is not falsely asserted for oscillating genes."""
    n = len(att_states)
    if n == 0:
        return _AttAgg(0, {})
    bn_nodes = list(att_states[0].keys())
    node_col = {g: i for i, g in enumerate(bn_nodes)}
    n_bn = len(bn_nodes)
    M = np.empty((n, n_bn + len(sink_order)), dtype=np.uint8)
    raw = [[att[g] for g in bn_nodes] for att in att_states]
    M[:, :n_bn] = np.array([[0 if v == '*' else int(v) for v in row] for row in raw], dtype=np.uint8)
    S = np.array([[v == '*' for v in row] for row in raw], dtype=np.uint8)
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


# ── ONE-TIME DATA LOAD ─────────────────────────────────────────────────────────
print(f"Loading: {SCRIPT_NAME}  |  network: {NETWORK_FILE}")
activators, suppressors, all_nodes, gene_category = {}, {}, set(), {}
_seen_act, _seen_sup = set(), set()
with open(os.path.join(BASE_DIR, "networks_used_by_scripts", NETWORK_FILE), newline="", encoding="utf-8") as f:
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
        if rel == ACT_REL: activators.setdefault(t, set()).add(s); _seen_act.add((s, t))
        else:              suppressors.setdefault(t, set()).add(s); _seen_sup.add((s, t))
edges_act = sorted(_seen_act); edges_sup = sorted(_seen_sup)
del _seen_act, _seen_sup
print(f"Loaded {len(edges_act)+len(edges_sup)} unique edges  ({len(all_nodes)} nodes)")

def cat(g): c = gene_category.get(g); return CAT_MAP.get(c, c or "?")
def mol_flag(g): return " !" if gene_category.get(g) not in MOLECULAR_CATS else ""

node_list = sorted(all_nodes)
node_idx  = {n: i for i, n in enumerate(node_list)}

t0 = time.perf_counter()
g_full = gt.Graph(directed=True)
g_full.add_vertex(len(node_list))
for s, t in edges_act: g_full.add_edge(node_idx[s], node_idx[t])
for s, t in edges_sup: g_full.add_edge(node_idx[s], node_idx[t])
print(f"graph-tool: {g_full.num_vertices()} vertices, {g_full.num_edges()} edges  ({time.perf_counter()-t0:.2f}s)")

pr_map  = gt.pagerank(g_full)
_, hist = gt.label_components(g_full)
print(f"Phase 1: PageRank top hub: {node_list[max(range(g_full.num_vertices()), key=lambda i: pr_map[i])]}, SCCs: {hist.shape[0]}, largest: {int(hist.max())}")

if GENE not in node_idx:
    print(f"ERROR: '{GENE}' not found in network. Exiting."); exit(1)
print(f"Gene '{GENE}' confirmed in network [{cat(GENE)}].")
print(f"Benchmark: time limit {MAX_TOTAL_SECONDS}s ({MAX_TOTAL_SECONDS//3600}h) | Mode: BoNesis (budget-limited)\n")

class _BudgetExhausted(Exception): pass

# ── BENCHMARK LOOP ─────────────────────────────────────────────────────────────
t_benchmark_start = time.perf_counter()
timings = []
hops = 0
while True:
    hops += 1
    if hops > MAX_HOPS:
        print(f"Safety cap of {MAX_HOPS} hops reached. Stopping."); break
    remaining = MAX_TOTAL_SECONDS - (time.perf_counter() - t_benchmark_start)
    if remaining <= 0:
        print(f"Time limit ({MAX_TOTAL_SECONDS}s / {MAX_TOTAL_SECONDS//3600}h) reached before hop {hops}. Stopping."); break
    print(f"  Hop {hops}... (remaining: {remaining:.0f}s)", end=" ", flush=True)
    t_start = time.perf_counter()
    source_gene = GENE
    out_file = os.path.join(OUT_DIR, f"{SCRIPT_NAME}_MYB46_hops{hops}.txt")

    try:
        with open(out_file, 'w') as fout, contextlib.redirect_stdout(fout):
            print(f"# Benchmark: {SCRIPT_NAME}.py")
            print(f"# Gene: {GENE}  |  Hops: {hops}  |  Mode: BoNesis (complete attractor landscape)")
            print(f"# Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}")
            print(f"\n{'='*70}")
            print(f"FORWARD SIGNALING ANALYSIS — METHODOLOGY")
            print(f"  Node of interest : {GENE}")
            print(f"  Network          : {NETWORK_FILE}")
            print(f"  Hop radius       : {hops}")
            print(f"  Attractor solver : BoNesis — complete Boolean network attractor landscape.")
            print(f"                     Finds ALL attractors reachable from each initial condition.")
            print(f"")
            print(f"  Subgraph: all nodes reachable from {GENE} in ≤{hops} directed steps.")
            print(f"  Boolean rules: any activator is sufficient (OR); any suppressor dominates (AND-NOT).")
            print(f"  Rule form: (act1 | act2 | ...) & !sup1 & !sup2 & ...")
            print(f"  Update scheme: synchronous (all nodes updated simultaneously each step).")
            print(f"")
            print(f"  Experimental conditions compared:")
            print(f"    (A) Dark background      — all nodes initialised to 0 (all pathways inactive)")
            print(f"    (B) Permissive background — all nodes initialised to 1 (all pathways active)")
            print(f"  Within each background, resting ({GENE}=0) vs. perturbed ({GENE}=1) are compared.")
            print(f"  A node is classified as:")
            print(f"    Activated  : stable OFF in resting AND stable ON  in ALL perturbed attractors")
            print(f"    Suppressed : stable ON  in resting AND stable OFF in ALL perturbed attractors")
            print(f"    Robust     : activated/suppressed consistently in BOTH backgrounds")
            print(f"")
            print(f"  Sink nodes: regulated by others but regulating nothing in the subgraph.")
            print(f"  Evaluated post-hoc by applying their Boolean rule to settled attractor states.")
            print(f"  Pre-filtered: only sinks that change state between conditions are reported.")
            print(f"  ! = experimental perturbation node (mutant/process/phenotype) — not a gene product")
            print(f"{'='*70}")

            dist = gt.shortest_distance(g_full, source=g_full.vertex(node_idx[source_gene]), directed=True)
            subgraph_nodes = {node_list[i] for i in range(g_full.num_vertices()) if dist[i] <= hops}
            sub_node_list  = sorted(subgraph_nodes)
            sub_v_idx      = {name: i for i, name in enumerate(sub_node_list)}
            print(f"  [time] Subgraph extraction: {time.perf_counter()-t_start:.3f}s")

            _sub_edges = sum(1 for s, t in edges_act + edges_sup if s in subgraph_nodes and t in subgraph_nodes)
            print(f"\n{'='*70}")
            print(f"PHASE 2 — {source_gene}, {hops} hop(s)  |  {len(subgraph_nodes)} nodes, {_sub_edges} edges")
            for c, n in sorted(Counter(cat(nd) for nd in subgraph_nodes).items(), key=lambda x: x[1], reverse=True):
                print(f"    {n:>4}  {c}")

            print(f"\n{'='*70}\nPHASE 2.5 — Subgraph topology")
            g_sub = gt.Graph(directed=True)
            g_sub.add_vertex(len(sub_node_list))
            for s, t in edges_act + edges_sup:
                if s in sub_v_idx and t in sub_v_idx:
                    g_sub.add_edge(sub_v_idx[s], sub_v_idx[t])
            sub_comp, sub_hist = gt.label_components(g_sub)
            large_sub = int((sub_hist > 1).sum())
            print(f"  SCCs: {sub_hist.shape[0]}  |  non-trivial: {large_sub}")
            if large_sub:
                scc_m = [sub_node_list[i] for i in range(g_sub.num_vertices()) if int(sub_comp[i]) == int(sub_hist.argmax())]
                print(f"  Largest SCC ({int(sub_hist.max())} nodes):")
                for m in sorted(scc_m)[:10]: print(f"    {m:40s}  [{cat(m)}]")
                if len(scc_m) > 10: print(f"    ... and {len(scc_m)-10} more")

            bn_dict = {}
            for target in set(activators) | set(suppressors):
                if target not in subgraph_nodes: continue
                act = " | ".join(s for s in activators.get(target, []) if s in subgraph_nodes)
                sup = " & ".join("!" + s for s in suppressors.get(target, []) if s in subgraph_nodes)
                if act and sup: bn_dict[target] = "(" + act + ") & " + sup
                elif act: bn_dict[target] = act
                elif sup: bn_dict[target] = sup
            referenced = set()
            for rule in bn_dict.values(): referenced |= set(re.findall(r'\b[a-zA-Z_]\w*\b', str(rule)))
            sink_nodes = {g for g in bn_dict if g not in referenced and g != source_gene}
            bn_dict_pruned = {g: f for g, f in bn_dict.items() if g not in sink_nodes}
            n_pruned = len(bn_dict_pruned)
            print(f"  Regulated: {len(bn_dict)}  |  pruned: {n_pruned}  |  sinks: {len(sink_nodes)}")

            bn_resting_dict   = dict(bn_dict_pruned); bn_resting_dict[source_gene]   = "0"
            bn_perturbed_dict = dict(bn_dict_pruned); bn_perturbed_dict[source_gene] = "1"
            dark_start = {g: 0 for g in bn_resting_dict}
            perm_start = {g: 1 for g in bn_resting_dict}

            _remaining = max(0, MAX_TOTAL_SECONDS - (time.perf_counter() - t_benchmark_start))
            print(f"  Attempting BoNesis ({n_pruned} nodes — budget remaining: {_remaining:.0f}s)...")
            t_bn = time.perf_counter()
            _bn_r = bonesis.BooleanNetwork(bn_resting_dict)
            _bn_p = bonesis.BooleanNetwork(bn_perturbed_dict)
            _box = []
            _t = threading.Thread(daemon=True, target=lambda: _box.append((
                list(_bn_r.attractors(reachable_from=dark_start)),
                list(_bn_p.attractors(reachable_from=dark_start)),
                list(_bn_r.attractors(reachable_from=perm_start)),
                list(_bn_p.attractors(reachable_from=perm_start)),
            )))
            _t.start(); _t.join(_remaining)
            if _box:
                dark_resting_att, dark_perturbed_att, perm_resting_att, perm_perturbed_att = _box[0]
                print(f"  BoNesis: {time.perf_counter()-t_bn:.2f}s  |  attractors: "
                      f"dark-rest={len(dark_resting_att)}  dark-pert={len(dark_perturbed_att)}  "
                      f"perm-rest={len(perm_resting_att)}  perm-pert={len(perm_perturbed_att)}")
            else:
                print(f"  BoNesis timed out — budget exhausted after {time.perf_counter()-t_bn:.0f}s")
                raise _BudgetExhausted()

            # Sink pre-filter: evaluate each sink once per representative state in each background.
            # Topological order handles sink-on-sink dependencies. A sink is kept only if its rule
            # evaluates differently between resting and perturbed in AT LEAST ONE background
            # (dark or permissive), so sinks relevant only in the permissive background are not missed.
            sink_rules = {g: bn_dict[g] for g in sink_nodes}
            _n_sinks_total = len(sink_nodes)
            _sink_pre_order = _topo_sort_sinks(sink_rules)
            _NB = {"__builtins__": {}}
            _ns_dr = {g: (False if v == '*' else bool(v)) for g, v in (dark_resting_att[0]   if dark_resting_att   else {}).items()}
            _ns_dp = {g: (False if v == '*' else bool(v)) for g, v in (dark_perturbed_att[0] if dark_perturbed_att else {}).items()}
            _ns_pr = {g: (False if v == '*' else bool(v)) for g, v in (perm_resting_att[0]   if perm_resting_att   else {}).items()}
            _ns_pp = {g: (False if v == '*' else bool(v)) for g, v in (perm_perturbed_att[0] if perm_perturbed_att else {}).items()}
            _ns_dr[source_gene] = False; _ns_dp[source_gene] = True
            _ns_pr[source_gene] = False; _ns_pp[source_gene] = True
            _sink_relevant = set()
            for _g in _sink_pre_order:
                _expr = _to_py(sink_rules[_g])
                try:    _rv_d = eval(_expr, _NB, _ns_dr)
                except: _rv_d = None
                try:    _pv_d = eval(_expr, _NB, _ns_dp)
                except: _pv_d = None
                try:    _rv_p = eval(_expr, _NB, _ns_pr)
                except: _rv_p = None
                try:    _pv_p = eval(_expr, _NB, _ns_pp)
                except: _pv_p = None
                _ns_dr[_g] = bool(_rv_d) if _rv_d is not None else False
                _ns_dp[_g] = bool(_pv_d) if _pv_d is not None else False
                _ns_pr[_g] = bool(_rv_p) if _rv_p is not None else False
                _ns_pp[_g] = bool(_pv_p) if _pv_p is not None else False
                if (_rv_d is None or _pv_d is None or _rv_d != _pv_d or
                        _rv_p is None or _pv_p is None or _rv_p != _pv_p):
                    _sink_relevant.add(_g)
            sink_nodes = _sink_relevant
            sink_rules = {g: sink_rules[g] for g in sink_nodes}
            print(f"  Sink pre-filter: {len(sink_nodes)} / {_n_sinks_total} sinks change state in at least one background")

            if len(sink_nodes) <= SINK_RECOVERY_THRESHOLD:
                _n_att_total = sum(len(a) for a in (dark_resting_att, dark_perturbed_att, perm_resting_att, perm_perturbed_att))
                print(f"  Sink recovery: {len(sink_nodes)} sinks × {_n_att_total} attractor states...")
                sink_order = _topo_sort_sinks(sink_rules)
            else:
                sink_order = []
                print(f"  Sink recovery skipped ({len(sink_nodes)} sinks > {SINK_RECOVERY_THRESHOLD}), sinks excluded from classification")

            _t_sr = time.perf_counter()
            dark_resting_full   = _recover_sinks_numpy(dark_resting_att,   sink_order, source_gene, 0)
            dark_perturbed_full = _recover_sinks_numpy(dark_perturbed_att, sink_order, source_gene, 1)
            perm_resting_full   = _recover_sinks_numpy(perm_resting_att,   sink_order, source_gene, 0)
            perm_perturbed_full = _recover_sinks_numpy(perm_perturbed_att, sink_order, source_gene, 1)
            if sink_order:
                print(f"  Sink recovery: {time.perf_counter()-_t_sr:.3f}s")

            direct_targets = sorted(
                {t for t in bn_resting_dict if source_gene in activators.get(t, set()) or source_gene in suppressors.get(t, set())}
                | {t for t in sink_nodes    if source_gene in activators.get(t, set()) or source_gene in suppressors.get(t, set())}
            )

            def stable_on(g, atts):
                if isinstance(atts, _AttAgg): return atts.stable_on(g)
                return bool(atts) and all(a.get(g, 0) == 1 for a in atts)
            def stable_off(g, atts):
                if isinstance(atts, _AttAgg): return atts.stable_off(g)
                return bool(atts) and all(a.get(g, 0) == 0 for a in atts)
            def slabel(g, atts): return "ON" if stable_on(g, atts) else "OFF" if stable_off(g, atts) else "var"
            def tags(g): return "[" + ", ".join([cat(g)] + (["sink"] if g in sink_nodes else [])) + "]" + mol_flag(g)

            G = set(bn_dict) | {source_gene}
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

            # KO from the active (perturbed) state — unified loop handles both fixed-point and
            # multi-attractor cases; fixed-point just runs once, cycle intersects across all states
            ko_all_conv   = True
            ko_cycle_note = (f" (cycle/{len(dark_perturbed_att)}: intersection across all "
                             f"{len(dark_perturbed_att)} starting states)") if len(dark_perturbed_att) > 1 else ""
            cand_maintained = {g for g in G if stable_on(g,  dark_perturbed_full) and g != source_gene}
            cand_suppressed = {g for g in G if stable_off(g, dark_perturbed_full) and g != source_gene}
            for cs in dark_perturbed_att:
                if not cand_maintained and not cand_suppressed: break
                ko_att, _, _conv = simulate(bn_dict_pruned, cs, source_gene, 0)
                if not _conv: ko_all_conv = False
                ko_full = recover_sinks(ko_att, sink_rules, source_gene, 0, sink_order)
                cand_maintained = {g for g in cand_maintained if not stable_on(g, ko_full)}
                cand_suppressed = {g for g in cand_suppressed if     stable_on(g, ko_full)}
            ko_maintained = sorted(cand_maintained)
            ko_suppressed = sorted(cand_suppressed)

            # Variable genes: nodes with different stable states across the perturbed attractors
            # (uses raw attractor list; sinks are evaluated separately via perm_perturbed_full)
            gstates = {}
            for att in perm_perturbed_att:
                for g, v in att.items(): gstates.setdefault(g, set()).add(v)
            variable_genes = sorted(g for g, vs in gstates.items()
                                    if len(vs) > 1 and all(isinstance(v, int) for v in vs))
            decisions = {}
            for g in variable_genes:
                p = tuple(a.get(g, 0) for a in perm_perturbed_att)
                if not all(isinstance(v, int) for v in p): continue
                decisions.setdefault(p, []).append(g)

            necessary, dispensable = [], []
            if perm_activated:
                print(f"\n  Running necessity test on {len(direct_targets)} direct target(s) (BoNesis)...")
                t0_nec = time.perf_counter()
                for candidate in direct_targets:
                    _nec_remaining = MAX_TOTAL_SECONDS - (time.perf_counter() - t_benchmark_start)
                    if _nec_remaining <= 0:
                        print(f"  Necessity test budget exhausted — results partial"); break
                    ko = dict(bn_perturbed_dict); ko[candidate] = "0"
                    _bn_ko   = bonesis.BooleanNetwork(ko)
                    _nec_box = []
                    _nt = threading.Thread(daemon=True, target=lambda _b=_bn_ko, _nb=_nec_box, _k=ko:
                        _nb.append(list(_b.attractors(reachable_from={g: 1 for g in _k}))))
                    _nt.start(); _nt.join(_nec_remaining)
                    if not _nec_box:
                        print(f"    {candidate}: BoNesis timed out — skipping"); continue
                    ko_states = _nec_box[0]
                    ko_full   = recover_sinks(ko_states, sink_rules, source_gene, 1, sink_order,
                                             extra_pins={candidate: 0})
                    lost = sorted(g for g in perm_activated if g != candidate and not stable_on(g, ko_full))
                    if lost: necessary.append((candidate, lost))
                    else:    dispensable.append(candidate)
                print(f"  Necessity test completed in {time.perf_counter()-t0_nec:.3f}s")
            else:
                print(f"\n  Necessity test skipped — no genes stably activated in permissive baseline.")

            # ── OUTPUT ────────────────────────────────────────────────────────
            print(f"\n{'='*70}")
            print(f"  {source_gene} — {hops} hop(s)  |  BoNesis (complete attractor landscape)")
            print(f"  dark: {len(dark_resting_att)} rest/{len(dark_perturbed_att)} pert att  |  "
                  f"perm: {len(perm_resting_att)} rest/{len(perm_perturbed_att)} pert att")
            print(f"\n  Direct targets of {source_gene}: {len(direct_targets)}")
            for g in direct_targets: print(f"    {g:40s}  {tags(g)}")
            print(f"\n{'='*70}")
            print(f"  EXP A — dark background  |  {len(dark_resting_att)} resting / {len(dark_perturbed_att)} perturbed attractor(s)")
            print(f"\n  Activated (OFF→ON): {len(dark_activated)}")
            for g in dark_activated: print(f"    + {g:40s}  {tags(g)}")
            print(f"\n  Suppressed (ON→OFF): {len(dark_suppressed)}")
            for g in dark_suppressed: print(f"    - {g:40s}  {tags(g)}")
            print(f"\n{'='*70}")
            print(f"  EXP B — permissive background  |  {len(perm_resting_att)} resting / {len(perm_perturbed_att)} perturbed attractor(s)")
            print(f"\n  Activated (OFF→ON): {len(perm_activated)}")
            for g in perm_activated: print(f"    + {g:40s}  {tags(g)}")
            print(f"\n  Suppressed (ON→OFF): {len(perm_suppressed)}")
            for g in perm_suppressed: print(f"    - {g:40s}  {tags(g)}")
            if conditional_derepressed or conditional_suppressed:
                print(f"\n  Conditional effects:")
                for g, count in conditional_derepressed:
                    print(f"    + ({count}/{len(perm_resting_att)} resting OFF) {g:36s}  {tags(g)}")
                for g, count in conditional_suppressed:
                    print(f"    - ({count}/{len(perm_resting_att)} resting ON)  {g:36s}  {tags(g)}")
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
                print(f"  NECESSITY TEST  (method: BoNesis from permissive background)")
                print(f"  Definition: a direct target is necessary if knocking it out causes at")
                print(f"  least one permissive-activated gene to lose stable-ON status.")
                print(f"    Necessary  ({len(necessary)}):")
                for g, lost in necessary:
                    print(f"    ! {g:40s}  {tags(g)}")
                    print(f"        loss of stable activation: {', '.join(lost)}")
                print(f"    Dispensable ({len(dispensable)}):")
                for g in dispensable: print(f"      {g:40s}  {tags(g)}")
            print(f"\n{'='*70}")
            n_perturbed_att = len(perm_perturbed_att)
            print(f"  Context-dependent genes: {len(variable_genes)} gene(s) in {len(decisions)} co-varying pattern(s)")
            if variable_genes:
                print(f"  (Different stable values across {n_perturbed_att} perturbed attractors; groups share identical ON/OFF patterns.)")
            for i, (pattern, genes) in enumerate(sorted(decisions.items()), 1):
                on_in  = [j + 1 for j, v in enumerate(pattern) if v == 1]
                off_in = [j + 1 for j, v in enumerate(pattern) if v == 0]
                print(f"\n  Pattern {i} — stably ON in attractor(s) {on_in or 'none'}  |  stably OFF in {off_in or 'none'}")
                for g in sorted(genes): print(f"    {g:40s}  {tags(g)}")
            print(f"\n{'='*70}")
            print(f"  Sink nodes: {len(sink_nodes)} of {_n_sinks_total} change state when {source_gene} is perturbed")
            for g in sorted(sink_nodes):
                r = slabel(g, perm_resting_full)
                p = slabel(g, perm_perturbed_full)
                print(f"    {g:40s}  resting: {r:3s}  perturbed: {p:3s}  {tags(g)}")
            print(f"\n{'='*70}")
            print(f"  TIMING SUMMARY — {source_gene}, hop {hops}")
            print(f"    Total hop time: {time.perf_counter()-t_start:.3f}s")
            print(f"{'='*70}")

        elapsed = time.perf_counter() - t_start
        timings.append((hops, elapsed, "BoNesis"))
        print(f"{elapsed:.1f}s [BoNesis] | {len(subgraph_nodes)} nodes → {os.path.basename(out_file)}")

    except _BudgetExhausted:
        elapsed = time.perf_counter() - t_start
        timings.append((hops, elapsed, "Timeout"))
        print(f"{elapsed:.1f}s [Timeout] | {len(subgraph_nodes)} nodes → budget exhausted, stopping.")
        break

    except Exception as e:
        import traceback
        elapsed = time.perf_counter() - t_start
        print(f"ERROR after {elapsed:.1f}s: {e}")
        timings.append((hops, elapsed, "ERROR"))
        with open(out_file, 'a') as fout: fout.write(f"\n--- ERROR ---\n{traceback.format_exc()}\n")
        break

print(f"\n{'='*70}")
print(f"BENCHMARK SUMMARY: {SCRIPT_NAME}")
print(f"Gene: {GENE}  |  Mode: BoNesis (budget-limited)")
print(f"{'Hops':>6}  {'Time (s)':>10}  {'Mode':>12}  Output file")
print(f"{'-'*6}  {'-'*10}  {'-'*12}  {'-'*45}")
for hops, t, mode in timings:
    print(f"{hops:>6}  {t:>10.2f}  {mode:>12}  {SCRIPT_NAME}_MYB46_hops{hops}.txt")
