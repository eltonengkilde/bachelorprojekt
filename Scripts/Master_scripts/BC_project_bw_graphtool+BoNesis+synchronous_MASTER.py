#!/opt/anaconda3/envs/bachelor_env/bin/python
import csv, os, re, sys, time, keyword, threading
import graph_tool.all as gt
import bonesis

BASE_DIR        = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUTPUT_DIR      = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")

class _Tee:
    def __init__(self, *files): self.files = files
    def write(self, data):
        for f in self.files: f.write(data)
    def flush(self):
        for f in self.files: f.flush()
BONESIS_TIMEOUT = 60    # seconds; BoNesis necessity runs first, simulation fallback if stalled
MAX_SIM_STEPS   = 1000  # max synchronous simulation steps before non-convergence warning
CATEGORIES_TO_KEEP = {
    "Gene / Protein",
    "Phenotype / Trait / Disease",
    "Chemical / Metabolite / Cofactor / Ligand",
    "Biological Process / Pathway / Function / Regulatory / Signaling Mechanism",
}
BOOLEAN_RESERVED = {"TRUE", "FALSE", "NOT", "AND", "OR"}
ACT_REL = "Activation / Induction / Causation / Result"
SUP_REL = "Repression / Inhibition / Negative Regulation"
CAT_MAP = {
    "Gene / Protein": "Gene",
    "Phenotype / Trait / Disease": "Phenotype",
    "Chemical / Metabolite / Cofactor / Ligand": "Metabolite",
    "Biological Process / Pathway / Function / Regulatory / Signaling Mechanism": "Pathway",
}

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
activators, suppressors, edges_act, edges_sup, all_nodes, gene_category = {}, {}, [], [], set(), {}
with open(os.path.join(BASE_DIR, "networks_used_by_scripts", "filtered_networkL_normalized.csv"), newline="", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        if row["source_category"] not in CATEGORIES_TO_KEEP or row["target_category"] not in CATEGORIES_TO_KEEP: continue
        rel = row["relationship_category"]
        if rel not in (ACT_REL, SUP_REL): continue
        s, t = clean_name(row["source"]), clean_name(row["target"])
        all_nodes.update((s, t))
        gene_category[s] = row["source_category"]; gene_category[t] = row["target_category"]
        if rel == ACT_REL: activators.setdefault(t, set()).add(s); edges_act.append((s, t))
        else:              suppressors.setdefault(t, set()).add(s); edges_sup.append((s, t))

print(f"Loaded {len(edges_act)+len(edges_sup)} edges  ({len(all_nodes)} nodes)")

# shorthand for entity type, print how many nodes fall into each category
def cat(g): c = gene_category.get(g); return CAT_MAP.get(c, c or "?")
cat_counts = {}
for c in gene_category.values(): cat_counts[c] = cat_counts.get(c, 0) + 1
print(f"\nEntity categories ({len(cat_counts)}):")
for c, n in sorted(cat_counts.items(), key=lambda x: x[1], reverse=True): print(f"  {n:>6}  {c}")

# build the graph-tool graph, activation edges flagged True and suppression edges False
node_list = sorted(all_nodes)
node_idx  = {n: i for i, n in enumerate(node_list)}
t_total_start = time.perf_counter()

t0 = time.perf_counter()
g_full = gt.Graph(directed=True)
g_full.add_vertex(len(node_list))
for s, t in edges_act: g_full.add_edge(node_idx[s], node_idx[t])
for s, t in edges_sup: g_full.add_edge(node_idx[s], node_idx[t])
print(f"\ngraph-tool: {g_full.num_vertices()} vertices, {g_full.num_edges()} edges  ({time.perf_counter()-t0:.2f}s)")

# phase 1 - characterise the full network, find hubs with PageRank, check for feedback loops
print(f"\n{'='*70}")
print(f"PHASE 1  nodes: {g_full.num_vertices()}  |  edges: {g_full.num_edges()} ({len(edges_act)} act / {len(edges_sup)} sup)")
t0 = time.perf_counter()
pr_map = gt.pagerank(g_full)
print(f"  PageRank ({time.perf_counter()-t0:.2f}s), top 10:")
for i in sorted(range(g_full.num_vertices()), key=lambda i: pr_map[i], reverse=True)[:10]:
    print(f"    {node_list[i]:40s}  PR: {pr_map[i]:.6f}  [{cat(node_list[i])}]")
_, hist = gt.label_components(g_full)
print(f"  SCCs: {hist.shape[0]}  |  non-trivial: {int((hist>1).sum())}  |  largest: {int(hist.max())} nodes")

# ask which node to analyse and how many hops upstream
while True:
    target_gene = input("\nTarget node for upstream regulatory analysis: ")
    if target_gene not in node_idx:
        matches = sorted(n for n in all_nodes if target_gene.upper() in n.upper())
        print(f"WARNING: '{target_gene}' not found. Try again.")
        if matches: print(f"  Did you mean one of these? {matches[:5]}")
    else:
        print(f"OK: '{target_gene}'  [{cat(target_gene)}]"); break

while True:
    MAX_HOPS = int(input("\nNumber of hops to explore upstream: "))
    if MAX_HOPS < 1: print("Must be a whole number larger than 0")
    else: print(f"OK: {MAX_HOPS} hops upstream from '{target_gene}'"); break

if input('\nType and enter "yolo" to run the analysis: ') != "yolo":
    print("Invalid input, exiting."); exit()

os.makedirs(OUTPUT_DIR, exist_ok=True)
_safe_gene = re.sub(r'[^a-zA-Z0-9_-]', '_', target_gene)
_outpath = os.path.join(OUTPUT_DIR, f"{_safe_gene}_{MAX_HOPS}hop_bw_{time.strftime('%Y%m%d_%H%M%S')}.txt")
_logfile = open(_outpath, "w", encoding="utf-8")
sys.stdout = _Tee(sys.__stdout__, _logfile)
print(f"Output → {_outpath}")

# reverse the graph and BFS from target node, collect all nodes that can reach it within MAX_HOPS
g_rev = gt.GraphView(g_full, reversed=True)
dist  = gt.shortest_distance(g_rev, source=g_full.vertex(node_idx[target_gene]), directed=True)
subgraph_nodes = {node_list[i] for i in range(g_full.num_vertices()) if dist[i] <= MAX_HOPS}
sub_node_list  = sorted(subgraph_nodes)
sub_v_idx      = {name: i for i, name in enumerate(sub_node_list)}
hop_of = {node_list[i]: int(dist[i]) for i in range(g_full.num_vertices()) if int(dist[i]) <= MAX_HOPS}

print(f"\n{'='*70}")
print(f"PHASE 2  Upstream subgraph: {target_gene}, {MAX_HOPS} hop(s)  |  {len(subgraph_nodes)} nodes, "
      f"{sum(1 for s,t in edges_act+edges_sup if s in sub_v_idx and t in sub_v_idx)} edges")
sc = {}
for n in subgraph_nodes: sc[cat(n)] = sc.get(cat(n), 0) + 1
for c, n in sorted(sc.items(), key=lambda x: x[1], reverse=True): print(f"    {n:>4}  {c}")

# characterise the subgraph - SCCs
print(f"\n{'='*70}")
print(f"PHASE 2.5  Subgraph topology (graph-tool)")
g_sub = gt.Graph(directed=True)
g_sub.add_vertex(len(sub_node_list))
for s, t in set((s, t) for s, t in edges_act if s in sub_v_idx and t in sub_v_idx):
    g_sub.add_edge(sub_v_idx[s], sub_v_idx[t])
for s, t in set((s, t) for s, t in edges_sup if s in sub_v_idx and t in sub_v_idx):
    g_sub.add_edge(sub_v_idx[s], sub_v_idx[t])

sub_comp, sub_hist = gt.label_components(g_sub)
large_sub = int((sub_hist > 1).sum())
print(f"  SCCs: {sub_hist.shape[0]}  |  non-trivial: {large_sub}")
if large_sub:
    scc_members = [sub_node_list[i] for i in range(g_sub.num_vertices()) if int(sub_comp[i]) == int(sub_hist.argmax())]
    print(f"  Largest SCC ({int(sub_hist.max())} nodes):")
    for m in sorted(scc_members)[:10]: print(f"    {m:40s}  [{cat(m)}]")
    if len(scc_members) > 10: print(f"    ... and {len(scc_members)-10} more")

# build one Boolean rule per regulated node, activators OR'd, suppressors AND NOT'd
# any single activator is enough to turn a node ON (functional redundancy)
# any suppressor keeps it OFF regardless of activators (dominant repression)
# rule form: (act1 | act2 | ...) & !sup1 & !sup2 & ...
bn_dict = {}
for tgt in set(activators) | set(suppressors):
    if tgt not in subgraph_nodes: continue
    act = " | ".join(s for s in activators.get(tgt, []) if s in subgraph_nodes)
    sup = " & ".join("!" + s for s in suppressors.get(tgt, []) if s in subgraph_nodes)
    if act and sup: bn_dict[tgt] = "(" + act + ") & " + sup
    elif act: bn_dict[tgt] = act
    elif sup: bn_dict[tgt] = sup

# in the upstream context, sinks have no incoming edges in this subgraph
# they are candidate master regulators since nothing drives them in this model
referenced = set()
for rule in bn_dict.values(): referenced |= set(re.findall(r'\b[a-zA-Z_]\w*\b', str(rule)))
sink_nodes = {gene for gene in bn_dict if gene not in referenced and gene != target_gene}
bn_dict_pruned = {g: f for g, f in bn_dict.items() if g not in sink_nodes}
n_pruned = len(bn_dict_pruned)
print(f"  Regulated: {len(bn_dict)}  |  pruned: {n_pruned}  |  sinks: {len(sink_nodes)}")

# helper functions for output labels, simulation and attractor checks
def tags(g):
    parts = [cat(g)]
    if g in sink_nodes: parts.append("sink")
    h = hop_of.get(g)
    if h is not None and h > 1: parts.append(f"hop {h}")
    return "[" + ", ".join(parts) + "]"

def simulate(rules, start_state, locked=None, val=None, max_steps=MAX_SIM_STEPS):
    genes = sorted(start_state.keys())
    idx   = {g: i for i, g in enumerate(genes)}
    li    = idx[locked] if locked is not None and locked in idx else None
    cc    = [(idx[g], compile(_to_py(r).strip(), "<string>", "eval"))
             for g, r in rules.items() if g != locked]
    state = [start_state.get(g, 0) for g in genes]
    if li is not None: state[li] = val
    history, seen = [], {}
    for i in range(max_steps):
        key = tuple(state)
        if key in seen:
            return [{g: s[idx[g]] for g in genes} for s in history[seen[key]:]], i, True
        seen[key] = i; history.append(list(state))
        ns = {g: bool(state[i]) for g, i in idx.items()}
        nw = list(state)
        if li is not None: nw[li] = val
        for gi, cd in cc: nw[gi] = 1 if eval(cd, {"__builtins__": {}}, ns) else 0
        state = nw
    return [{g: state[idx[g]] for g in genes}], max_steps, False

def target_on_any(atts):  return any(a.get(target_gene, 0) == 1 for a in atts)
def stable_on(g, atts):  return bool(atts) and all(a.get(g, 0) == 1 for a in atts)
def stable_off(g, atts): return bool(atts) and all(a.get(g, 0) == 0 for a in atts)

# two starting conditions: all nodes OFF (dark) and all nodes ON (permissive)
all_zeros_sub = {gene: 0 for gene in subgraph_nodes}
all_ones_sub  = {gene: 1 for gene in subgraph_nodes}

# attractor trace always uses synchronous simulation; necessity method determined at runtime
trace_mode = f"synchronous simulation ({n_pruned} nodes)"
print(f"  Attractor trace: {trace_mode}")

direct_activators  = sorted(g for g in activators.get(target_gene, set()) if g in subgraph_nodes)
direct_suppressors = sorted(g for g in suppressors.get(target_gene, set()) if g in subgraph_nodes)
all_upstream       = sorted(subgraph_nodes - {target_gene})

# STEP 1: four baseline simulations — two backgrounds × two target states
# Each pair holds the initial state constant and only varies the target lock.
# This ensures any differential is caused by the target, not the global state flip.
# dark background (all-zeros): target=0 vs target=1
# perm background (all-ones):  target=0 vs target=1
t0 = time.perf_counter()
dark_att_off, _, dark_off_conv = simulate(bn_dict_pruned, all_zeros_sub, locked=target_gene, val=0)
dark_att_on,  _, dark_on_conv  = simulate(bn_dict_pruned, all_zeros_sub, locked=target_gene, val=1)
perm_att_off, _, perm_off_conv = simulate(bn_dict_pruned, all_ones_sub,  locked=target_gene, val=0)
perm_att_on,  _, perm_on_conv  = simulate(bn_dict_pruned, all_ones_sub,  locked=target_gene, val=1)
print(f"\n  Baseline simulations (4): {time.perf_counter()-t0:.3f}s")
for _lbl, _conv in [("dark/off", dark_off_conv), ("dark/on", dark_on_conv),
                    ("perm/off", perm_off_conv), ("perm/on",  perm_on_conv)]:
    if not _conv: print(f"  WARNING: {_lbl} attractor did not converge within {MAX_SIM_STEPS} steps")
print(f"  NOTE: simulation finds ONE attractor per starting condition, "
      f"use BoNesis at fewer hops for the full attractor landscape")
target_in_perm = target_on_any(perm_att_on)
_att_label = lambda atts, conv: (
    ("fixed point" if len(atts) == 1 else f"cycle/{len(atts)}") +
    (" (converged)" if conv else " (max steps)"))
print(f"  Dark/OFF (all-zeros, target=0): {_att_label(dark_att_off, dark_off_conv)}")
print(f"  Dark/ON  (all-zeros, target=1): {_att_label(dark_att_on,  dark_on_conv)}")
print(f"  Perm/OFF (all-ones,  target=0): {_att_label(perm_att_off, perm_off_conv)}")
print(f"  Perm/ON  (all-ones,  target=1): {_att_label(perm_att_on,  perm_on_conv)}")

# STEP 2: hop-by-hop attractor state trace
# a node is target-responsive if its stable state flips when the target flips,
# holding the background (all-zeros or all-ones) constant.
# robust: responsive in BOTH backgrounds — genuine feedback regardless of context.
all_hops = sorted({hop_of[g] for g in all_upstream if hop_of.get(g) is not None})
hop_dark_on  = {h: [] for h in all_hops}   # target-responsive in dark background (OFF→ON)
hop_perm_on  = {h: [] for h in all_hops}   # target-responsive in perm background (OFF→ON)
hop_robust   = {h: [] for h in all_hops}   # responsive in both backgrounds
for g in all_upstream:
    h = hop_of.get(g)
    if h is None: continue
    dark_resp = stable_off(g, dark_att_off) and stable_on(g, dark_att_on)
    perm_resp = stable_off(g, perm_att_off) and stable_on(g, perm_att_on)
    if dark_resp: hop_dark_on[h].append(g)
    if perm_resp: hop_perm_on[h].append(g)
    if dark_resp and perm_resp: hop_robust[h].append(g)

# STEP 3: filter candidates by attractor state before per-node testing
# necessity (activator): stably ON in perm_att_on — these are candidates that could be required
# necessity (suppressor): stably OFF in perm_att_on — their absence may be required for target ON
# sufficiency: genuinely flips OFF→ON when target flips 0→1 from dark background
# suppressor release: only tested when target is OFF in perm_att_on
perm_on_stable  = sorted(g for g in all_upstream if stable_on(g,  perm_att_on))
perm_off_stable = sorted(g for g in all_upstream if stable_off(g, perm_att_on))

necessity_act_candidates = perm_on_stable if target_in_perm else []
sufficiency_candidates   = sorted(
    g for g in all_upstream
    if stable_off(g, dark_att_off) and stable_on(g, dark_att_on)
) if target_in_perm else []
necessity_sup_candidates = perm_off_stable if target_in_perm else []
suppressor_release_cands = perm_off_stable if not target_in_perm else []

print(f"\n  Attractor-filtered candidate pools:")
print(f"    Sufficiency (target-responsive in dark):  {len(sufficiency_candidates)}")
print(f"    Necessity, activators (stably ON-perm):   {len(necessity_act_candidates)}")
print(f"    Necessity, suppressors (stably OFF-perm): {len(necessity_sup_candidates)}")

# STEP 4: per-node tests
GENE_TEST_BUDGET = 120   # total seconds for sufficiency + simulation fallback
t_gene_tests = time.perf_counter()
budget_exceeded = False

# sufficiency: fix candidate to ON from all-zeros; does target turn ON?
sufficient_activators = []
print(f"\n  Sufficiency test: {len(sufficiency_candidates)} candidates (budget: {GENE_TEST_BUDGET}s)...")
t0 = time.perf_counter()
for candidate in sufficiency_candidates:
    if time.perf_counter() - t_gene_tests > GENE_TEST_BUDGET:
        print(f"    Stopped: budget reached after {time.perf_counter()-t_gene_tests:.0f}s")
        budget_exceeded = True; break
    states, _, _ = simulate(bn_dict_pruned, all_zeros_sub, candidate, 1)
    if target_on_any(states) and not target_on_any(dark_att_off):
        sufficient_activators.append(candidate)
print(f"  Sufficiency: {time.perf_counter()-t0:.3f}s")

# necessity and suppressor tests - BoNesis first, simulation fallback if it times out
necessary_activators, redundant_activators, necessary_suppressors, suppressor_releases = [], [], [], []
print(f"\n  BoNesis necessity test: {len(necessity_act_candidates)} act + {len(necessity_sup_candidates)} sup candidates ({n_pruned} nodes, timeout: {BONESIS_TIMEOUT}s)...")
t_bn = time.perf_counter(); _bn_box = []
_ones_p = {g: 1 for g in bn_dict_pruned}
def _run_bonesis_necessity():
    baseline = target_on_any(list(bonesis.BooleanNetwork(bn_dict_pruned).attractors(reachable_from=_ones_p)))
    _nec, _red, _nsup, _sup = [], [], [], []
    for cand in necessity_act_candidates:
        ko = dict(bn_dict_pruned); ko[cand] = "0"
        ko_t = target_on_any(list(bonesis.BooleanNetwork(ko).attractors(reachable_from={g: 1 for g in ko})))
        cb = baseline
        if cand not in bn_dict_pruned:
            pd = dict(bn_dict_pruned); pd[cand] = "1"
            cb = target_on_any(list(bonesis.BooleanNetwork(pd).attractors(reachable_from={g: 1 for g in pd})))
        if cb and not ko_t: _nec.append(cand)
        elif cb and ko_t:   _red.append(cand)
    for cand in necessity_sup_candidates:
        forced = dict(bn_dict_pruned); forced[cand] = "1"
        f_t = target_on_any(list(bonesis.BooleanNetwork(forced).attractors(reachable_from={g: 1 for g in forced})))
        if baseline and not f_t: _nsup.append(cand)
    for cand in suppressor_release_cands:
        ko = dict(bn_dict_pruned); ko[cand] = "0"
        ko_t = target_on_any(list(bonesis.BooleanNetwork(ko).attractors(reachable_from={g: 1 for g in ko})))
        if not baseline and ko_t: _sup.append(cand)
    _bn_box.append((baseline, _nec, _red, _nsup, _sup))
# daemon thread: if BoNesis stalls past BONESIS_TIMEOUT the thread is left running in the background
# it keeps running until the process exits, there is no clean way to stop a running solver
_t = threading.Thread(daemon=True, target=_run_bonesis_necessity)
_t.start(); _t.join(BONESIS_TIMEOUT)
if _bn_box:
    using_bonesis = True
    baseline_perm_on, necessary_activators, redundant_activators, necessary_suppressors, suppressor_releases = _bn_box[0]
    print(f"  BoNesis necessity: {time.perf_counter()-t_bn:.3f}s")
else:
    using_bonesis = False
    print(f"  BoNesis timed out after {BONESIS_TIMEOUT}s, falling back to synchronous simulation")
    baseline_perm_on = target_in_perm
    t_sim = time.perf_counter()
    for candidate in necessity_act_candidates:
        if time.perf_counter() - t_gene_tests > GENE_TEST_BUDGET:
            print(f"    Stopped: budget reached after {time.perf_counter()-t_gene_tests:.0f}s")
            budget_exceeded = True; break
        ko = dict(bn_dict_pruned); ko[candidate] = "0"
        ko_st, _, _ = simulate(ko, all_ones_sub, candidate, 0)
        if baseline_perm_on and not target_on_any(ko_st): necessary_activators.append(candidate)
        elif baseline_perm_on and target_on_any(ko_st):   redundant_activators.append(candidate)
    for candidate in necessity_sup_candidates:
        if time.perf_counter() - t_gene_tests > GENE_TEST_BUDGET:
            print(f"    Stopped: budget reached after {time.perf_counter()-t_gene_tests:.0f}s")
            budget_exceeded = True; break
        forced = dict(bn_dict_pruned); forced[candidate] = "1"
        f_st, _, _ = simulate(forced, all_ones_sub, candidate, 1)
        if baseline_perm_on and not target_on_any(f_st): necessary_suppressors.append(candidate)
    for candidate in suppressor_release_cands:
        if time.perf_counter() - t_gene_tests > GENE_TEST_BUDGET:
            print(f"    Stopped: budget reached after {time.perf_counter()-t_gene_tests:.0f}s")
            budget_exceeded = True; break
        ko = dict(bn_dict_pruned); ko[candidate] = "0"
        ko_st, _, _ = simulate(ko, all_ones_sub, candidate, 0)
        if not baseline_perm_on and target_on_any(ko_st): suppressor_releases.append(candidate)
    print(f"  Simulation fallback: {time.perf_counter()-t_sim:.3f}s")
nec_mode = f"BoNesis ({n_pruned} nodes)" if using_bonesis else f"synchronous simulation ({n_pruned} nodes)"

# evaluate a Boolean rule string directly, used to get sink node states for output
def eval_rule_simple(rule, state):
    try: return 1 if eval(_to_py(rule), {"__builtins__": {}}, {g: bool(v) for g, v in state.items()}) else 0
    except: return 0

sink_rules = {g: bn_dict[g] for g in sink_nodes}

# print run metadata
print("\n" + "="*70)
print(f"  RUN METADATA")
print(f"  Network       : filtered_networkL_normalized.csv")
print(f"  Target node   : {target_gene}  [{cat(target_gene)}]")
print(f"  Hops upstream : {MAX_HOPS}")
print(f"  Attractor trace   : {trace_mode}")
print(f"  Necessity / suff  : {nec_mode}")
if not using_bonesis:
    print(f"  NOTE: simulation necessity tests one attractor, results are indicative not guaranteed")
print(f"  Subgraph      : {len(subgraph_nodes)} nodes  |  pruned BN: {n_pruned}  |  sinks: {len(sink_nodes)}")
print(f"  Update scheme : synchronous (all nodes updated simultaneously each step)")
print(f"  Boolean rules : activators combined with OR; suppressors combined with AND NOT")
print(f"  Run time      : {time.strftime('%Y-%m-%d %H:%M:%S')}")
print("\n" + "="*70)
print(f"  UPSTREAM REGULATORY ANALYSIS of '{target_gene}'")
print(f"  {MAX_HOPS} hop(s) upstream  |  {len(all_upstream)} upstream nodes")

# OUTPUT STEP 1: structural direct regulators
print(f"\n{'='*70}")
print(f"  Step 1: Structural direct regulators of '{target_gene}' (1-hop edges in network)")
print(f"    Direct activators  ({len(direct_activators)}):")
for gene in direct_activators:  print(f"      + {gene:40s}  {tags(gene)}")
print(f"    Direct suppressors ({len(direct_suppressors)}):")
for gene in direct_suppressors: print(f"      - {gene:40s}  {tags(gene)}")

# OUTPUT STEP 2: attractor state trace
print(f"\n{'='*70}")
print(f"  Step 2: Target-responsive nodes (attractor state flips when target flips)")
print(f"  Each background is tested independently: only the target lock changes.")
print(f"  Robust = responsive in both dark and permissive backgrounds.")
print(f"    Dark/OFF (all-zeros, target=0): {_att_label(dark_att_off, dark_off_conv)}")
print(f"    Dark/ON  (all-zeros, target=1): {_att_label(dark_att_on,  dark_on_conv)}")
print(f"    Perm/OFF (all-ones,  target=0): {_att_label(perm_att_off, perm_off_conv)}")
print(f"    Perm/ON  (all-ones,  target=1): {_att_label(perm_att_on,  perm_on_conv)}")
print(f"\n  Target-responsive nodes by hop:")
print(f"  {'Hop':>4}  {'Total':>6}  {'Dark-resp':>10}  {'Perm-resp':>10}  {'Robust':>8}")
for h in all_hops:
    genes_h = [g for g in all_upstream if hop_of.get(g) == h]
    nd = len(hop_dark_on[h]); np_ = len(hop_perm_on[h]); nr = len(hop_robust[h])
    print(f"  {h:>4}  {len(genes_h):>6}  {nd:>10}  {np_:>10}  {nr:>8}")
    for gene in sorted(hop_robust[h])[:8]:
        print(f"      ~ {gene:40s}  [robust target-responsive]  {tags(gene)}")
    if nr > 8: print(f"      ... and {nr-8} more robust target-responsive at hop {h}")

# OUTPUT STEP 3: sufficient activators
print(f"\n{'='*70}")
print(f"  Step 3: Sufficient upstream activators")
print(f"  Definition: forcing the node to ON from the dark (all-zeros) initial condition")
print(f"  causes '{target_gene}' to be ON in the resulting attractor, when it would otherwise")
print(f"  be OFF. Tested on {len(sufficiency_candidates)} differential candidates.")
if budget_exceeded: print(f"  NOTE: {GENE_TEST_BUDGET}s budget reached, results may be partial.")
print(f"  Found: {len(sufficient_activators)}")
for gene in sufficient_activators: print(f"    + {gene:40s}  {tags(gene)}")

# OUTPUT STEP 4: necessary activators
print(f"\n{'='*70}")
print(f"  Step 4: Necessary upstream regulators (knockout from permissive background)")
print(f"  Method: {nec_mode}")
if budget_exceeded: print(f"  NOTE: {GENE_TEST_BUDGET}s budget reached, results may be partial.")
print(f"\n  Necessary activators, KO turns '{target_gene}' OFF ({len(necessary_activators)}):")
for gene in necessary_activators: print(f"    ! {gene:40s}  [required activator]  {tags(gene)}")
print(f"\n  Redundant activators, '{target_gene}' stays ON after KO ({len(redundant_activators)}):")
for gene in redundant_activators[:20]: print(f"    ~ {gene:40s}  [dispensable]  {tags(gene)}")
if len(redundant_activators) > 20: print(f"    ... and {len(redundant_activators)-20} more [dispensable]")

# OUTPUT STEP 5: suppressor tests
print(f"\n{'='*70}")
print(f"  Step 5: Suppressor necessity and release")
print(f"  Necessary suppressors: nodes whose ABSENCE (stably OFF in permissive attractor)")
print(f"  is required for '{target_gene}' to be ON. Tested by forcing the node to ON;")
print(f"  target turning OFF confirms the node is a necessary upstream suppressor.")
print(f"\n  Necessary suppressors, forced ON turns '{target_gene}' OFF ({len(necessary_suppressors)}):")
for gene in necessary_suppressors: print(f"    ! {gene:40s}  [necessary suppressor]  {tags(gene)}")
if suppressor_releases:
    print(f"\n  Suppressor release, KO turns '{target_gene}' ON ({len(suppressor_releases)}):")
    for gene in suppressor_releases: print(f"    ~ {gene:40s}  [suppressor release]  {tags(gene)}")

# OUTPUT: sink nodes
print(f"\n{'='*70}")
print(f"  Sink nodes: upstream nodes with no further upstream regulators in this subgraph")
print(f"  ({len(sink_nodes)} nodes). These are candidate master regulators (no in-edges in model).")
if sink_nodes:
    for gene in sorted(sink_nodes)[:30]:
        val = eval_rule_simple(sink_rules[gene], perm_att_on[0])
        print(f"    {gene:40s}  permissive state: {'ON' if val else 'OFF':3s}  {tags(gene)}")
    if len(sink_nodes) > 30: print(f"    ... and {len(sink_nodes)-30} more")

print(f"\n  Total analysis time: {time.perf_counter() - t_total_start:.2f}s")
print(f"\n{'='*70}")

sys.stdout = sys.__stdout__
_logfile.close()
print(f"\nOutput saved to: {_outpath}")
