#!/opt/anaconda3/envs/bachelor_env/bin/python
import csv, os, re, heapq, time, keyword, threading
import graph_tool.all as gt
import bonesis

BASE_DIR             = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BONESIS_TIMEOUT      = 60    # seconds; BoNesis runs first, simulation fallback if stalled
MAX_SIM_STEPS        = 1000  # max synchronous simulation steps before non-convergence warning
SINK_RECOVERY_THRESHOLD = 10000
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

def clean_name(name):
    if not name or not isinstance(name, str): return "unknown"
    name = re.sub(r"_+", "_", re.sub(r"[^a-zA-Z0-9_]", "_", name)).strip("_")
    if not name: return "unknown"
    if name[0].isdigit(): return "g" + name
    if name.upper() in BOOLEAN_RESERVED or keyword.iskeyword(name): return "gene_" + name
    return name

def _to_py(r): return r.replace("!", " not ").replace("|", " or ").replace("&", " and ")

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

def cat(g): c = gene_category.get(g); return CAT_MAP.get(c, c or "?")

cat_counts = {}
for c in gene_category.values(): cat_counts[c] = cat_counts.get(c, 0) + 1
print(f"\nEntity categories ({len(cat_counts)}):")
for c, n in sorted(cat_counts.items(), key=lambda x: x[1], reverse=True): print(f"  {n:>6}  {c}")

t0 = time.perf_counter()
node_list = sorted(all_nodes)
node_idx = {name: i for i, name in enumerate(node_list)}
g_full = gt.Graph(directed=True)
g_full.add_vertex(len(node_list))
gt_etype = g_full.new_edge_property("bool")
for s, t in edges_act: gt_etype[g_full.add_edge(node_idx[s], node_idx[t])] = True
for s, t in edges_sup: gt_etype[g_full.add_edge(node_idx[s], node_idx[t])] = False
g_full.ep["etype"] = gt_etype
print(f"\ngraph-tool: {g_full.num_vertices()} vertices, {g_full.num_edges()} edges  ({time.perf_counter()-t0:.2f}s)")

print(f"\n{'='*70}")
print(f"PHASE 1 — nodes: {len(all_nodes)}  |  edges: {len(edges_act)+len(edges_sup)} ({len(edges_act)} act / {len(edges_sup)} sup)")
t0 = time.perf_counter()
pr = gt.pagerank(g_full)
print(f"  PageRank ({time.perf_counter()-t0:.2f}s) — top 10:")
for i in heapq.nlargest(10, range(g_full.num_vertices()), key=lambda i: pr[i]):
    print(f"    {node_list[i]:40s}  PR: {pr[i]:.6f}  [{cat(node_list[i])}]")
_, hist = gt.label_components(g_full)
print(f"  SCCs: {hist.shape[0]}  |  non-trivial: {int((hist>1).sum())}  |  largest: {int(hist.max())} nodes")

t_total_start = time.perf_counter()

while True:
    source_gene = input("\nSource gene to investigate: ")
    if source_gene not in all_nodes:
        matches = [n for n in all_nodes if source_gene.upper() in n.upper()]
        print(f"WARNING: '{source_gene}' not found.")
        if matches: print(f"  Did you mean: {sorted(matches)[:5]}")
    else:
        print(f"OK: '{source_gene}'  [{cat(source_gene)}]"); break

while True:
    MAX_HOPS = int(input("\nHops to explore downstream: "))
    if MAX_HOPS < 1: print("Must be >= 1")
    else: print(f"OK: {MAX_HOPS} hop(s)"); break

if input('\nType "yolo" to run: ') != "yolo":
    print("Exiting."); exit()

dist = gt.shortest_distance(g_full, source=g_full.vertex(node_idx[source_gene]), directed=True)
subgraph_nodes = {node_list[i] for i in range(g_full.num_vertices()) if dist[i] <= MAX_HOPS}
sub_node_list  = sorted(subgraph_nodes)
sub_node_idx   = {name: i for i, name in enumerate(sub_node_list)}

print(f"\n{'='*70}")
print(f"PHASE 2 — {source_gene}, {MAX_HOPS} hop(s)  |  {len(subgraph_nodes)} nodes, "
      f"{sum(1 for s,t in edges_act+edges_sup if s in subgraph_nodes and t in subgraph_nodes)} edges")
sc = {}
for n in subgraph_nodes: sc[cat(n)] = sc.get(cat(n), 0) + 1
for c, n in sorted(sc.items(), key=lambda x: x[1], reverse=True): print(f"    {n:>4}  {c}")

print(f"\n{'='*70}")
print(f"PHASE 2.5 — Subgraph topology (graph-tool)")
g_sub = gt.Graph(directed=True)
g_sub.add_vertex(len(sub_node_list))
sub_etype = g_sub.new_edge_property("bool")
for s, t in set((s, t) for s, t in edges_act if s in sub_node_idx and t in sub_node_idx):
    sub_etype[g_sub.add_edge(sub_node_idx[s], sub_node_idx[t])] = True
for s, t in set((s, t) for s, t in edges_sup if s in sub_node_idx and t in sub_node_idx):
    sub_etype[g_sub.add_edge(sub_node_idx[s], sub_node_idx[t])] = False
g_sub.ep["etype"] = sub_etype

sub_comp, sub_hist = gt.label_components(g_sub)
large_sub = int((sub_hist > 1).sum())
print(f"  SCCs (Strongly Connected Components): {sub_hist.shape[0]}  |  non-trivial: {large_sub}")
if large_sub:
    scc_members = [sub_node_list[i] for i in range(g_sub.num_vertices()) if int(sub_comp[i]) == int(sub_hist.argmax())]
    print(f"  Largest SCC ({int(sub_hist.max())} nodes):")
    for m in sorted(scc_members)[:10]: print(f"    {m:40s}  [{cat(m)}]")
    if len(scc_members) > 10: print(f"    ... and {len(scc_members)-10} more")

community = {}
if len(subgraph_nodes) <= 5000:
    t0 = time.perf_counter()
    vb, _ = gt.betweenness(g_sub)
    print(f"\n  Betweenness ({time.perf_counter()-t0:.3f}s) — top 10:")
    for name in heapq.nlargest(10, sub_node_list, key=lambda n: vb[sub_node_idx[n]]):
        score = vb[sub_node_idx[name]]
        if score > 0: print(f"    {name:40s}  {score:.6f}  [{cat(name)}]")
else:
    print(f"  Betweenness skipped ({len(subgraph_nodes)} nodes > 5000)")

if len(subgraph_nodes) <= 2000:
    t0 = time.perf_counter()
    print(f"\n  Running SBM (Stochastic Block Model)...")
    state = gt.minimize_blockmodel_dl(g_sub)
    b = state.get_blocks()
    for name in sub_node_list: community[name] = int(b[sub_node_idx[name]])
    comm_sizes = {}
    for c in community.values(): comm_sizes[c] = comm_sizes.get(c, 0) + 1
    print(f"  SBM ({time.perf_counter()-t0:.2f}s) — {len(comm_sizes)} modules:")
    for cid, sz in sorted(comm_sizes.items(), key=lambda x: x[1], reverse=True)[:5]:
        cats = {}
        for m in (n for n, c in community.items() if c == cid): cats[cat(m)] = cats.get(cat(m), 0) + 1
        print(f"    Module {cid}: {sz:>4} nodes  ({', '.join(f'{v} {k}' for k,v in sorted(cats.items(), key=lambda x: x[1], reverse=True))})")
else:
    print(f"\n  SBM (Stochastic Block Model) skipped ({len(subgraph_nodes)} nodes > 2000)")

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
all_off = {g: 0 for g in subgraph_nodes}
all_on  = {g: 1 for g in subgraph_nodes}

def build_simulator(rules, genes, locked, val):
    idx = {g: i for i, g in enumerate(genes)}
    li  = idx[locked]
    cc  = [(idx[g], compile(_to_py(r).strip(), "<string>", "eval"))
           for g, r in rules.items() if g != locked]
    def step(st):
        ns = {g: bool(st[i]) for g, i in idx.items()}
        nw = list(st); nw[li] = val
        for gi, cd in cc: nw[gi] = 1 if eval(cd, {"__builtins__": {}}, ns) else 0
        return nw
    return step, idx

def simulate(rules, start_state, locked, val, max_steps=MAX_SIM_STEPS):
    genes = sorted(start_state.keys())
    step_fn, idx = build_simulator(rules, genes, locked, val)
    state = [start_state.get(g, 0) for g in genes]; state[idx[locked]] = val
    history, seen = [], {}
    for i in range(max_steps):
        key = tuple(state)
        if key in seen:
            return [{g: s[idx[g]] for g in genes} for s in history[seen[key]:]], i, True
        seen[key] = i; history.append(list(state)); state = step_fn(state)
    return [{g: state[idx[g]] for g in genes}], max_steps, False

def _sim_label(states, conv):
    return ("fixed point" if len(states) == 1 else f"cycle/{len(states)}") + (" (conv)" if conv else " (max steps)")

bn_resting   = bonesis.BooleanNetwork(bn_resting_dict)
bn_perturbed = bonesis.BooleanNetwork(bn_perturbed_dict)
dark_start   = {g: 0 for g in bn_dict_pruned}
shared_start = {g: 1 for g in bn_dict_pruned}
print(f"  Attempting BoNesis ({n_pruned} nodes, timeout: {BONESIS_TIMEOUT}s)...")
_bonesis_ok = False; t0 = time.perf_counter(); _box = []
_t = threading.Thread(daemon=True, target=lambda: _box.append((
    list(bn_resting.attractors(reachable_from=dark_start)),
    list(bn_perturbed.attractors(reachable_from=dark_start)),
    list(bn_resting.attractors(reachable_from=shared_start)),
    list(bn_perturbed.attractors(reachable_from=shared_start)),
)))
_t.start(); _t.join(BONESIS_TIMEOUT)
if _box:
    dark_resting_att, dark_perturbed_att, perm_resting_att, perm_perturbed_att = _box[0]
    print(f"  BoNesis: {time.perf_counter()-t0:.2f}s")
    _bonesis_ok = True
else:
    print(f"  BoNesis timed out ({BONESIS_TIMEOUT}s) — falling back to synchronous simulation")
    t0 = time.perf_counter()
    dr_states, _, dr_conv = simulate(bn_dict_pruned, all_off, source_gene, 0)
    dp_states, _, dp_conv = simulate(bn_dict_pruned, all_off, source_gene, 1)
    pr_states, _, pr_conv = simulate(bn_dict_pruned, all_on,  source_gene, 0)
    pp_states, _, pp_conv = simulate(bn_dict_pruned, all_on,  source_gene, 1)
    print(f"  dark  rest: {_sim_label(dr_states, dr_conv)}  pert: {_sim_label(dp_states, dp_conv)}")
    print(f"  perm  rest: {_sim_label(pr_states, pr_conv)}  pert: {_sim_label(pp_states, pp_conv)}")
    print(f"  Simulation: {time.perf_counter()-t0:.3f}s")
    if not all([dr_conv, dp_conv, pr_conv, pp_conv]):
        print(f"  WARNING: not all runs converged within {MAX_SIM_STEPS} steps")
    print(f"  NOTE: simulation finds ONE attractor per condition — use bonesis (fewer hops) for full landscape.")
    dark_resting_att, dark_perturbed_att, perm_resting_att, perm_perturbed_att = dr_states, dp_states, pr_states, pp_states

using_simulation = not _bonesis_ok

def eval_rule_simple(rule, state):
    expr = _to_py(rule)
    try: return 1 if eval(expr, {"__builtins__": {}}, {g: bool(v) for g, v in state.items()}) else 0
    except: return 0

def _topo_sort_sinks(sink_rules):
    """Topological evaluation order for sink genes respecting intra-sink dependencies.
    If sink B's rule references sink A, A is placed before B so recover_sinks
    uses the recovered value of A when evaluating B, not the raw attractor value."""
    sink_set = set(sink_rules)
    deps = {g: set(re.findall(r'\b[a-zA-Z_]\w*\b', str(sink_rules[g]))) & sink_set
            for g in sink_set}
    dependents = {g: [] for g in sink_set}
    in_degree  = {g: 0  for g in sink_set}
    for g in sink_set:
        for d in deps[g]:
            dependents[d].append(g)
            in_degree[g] += 1
    queue = [g for g in sorted(sink_set) if in_degree[g] == 0]
    order = []
    while queue:
        node = queue.pop(0); order.append(node)
        for dep in sorted(dependents[node]):
            in_degree[dep] -= 1
            if in_degree[dep] == 0: queue.append(dep)
    if len(order) < len(sink_set):
        order.extend(g for g in sorted(sink_set) if g not in set(order))
    return order

def recover_sinks(attractors, sink_rules, src, src_val, sink_order):
    """Evaluate sink gene states from attractor in topological order so that
    a sink whose rule references another sink uses the recovered value."""
    result = []
    for att in attractors:
        ext = {**att, src: src_val}
        for sink in sink_order:
            ext[sink] = eval_rule_simple(sink_rules[sink], ext)
        result.append(ext)
    return result

sink_rules = {g: bn_dict[g] for g in sink_nodes}
if len(sink_nodes) <= SINK_RECOVERY_THRESHOLD:
    sink_order = _topo_sort_sinks(sink_rules)
    dark_resting_full   = recover_sinks(dark_resting_att,   sink_rules, source_gene, 0, sink_order)
    dark_perturbed_full = recover_sinks(dark_perturbed_att, sink_rules, source_gene, 1, sink_order)
    perm_resting_full   = recover_sinks(perm_resting_att,   sink_rules, source_gene, 0, sink_order)
    perm_perturbed_full = recover_sinks(perm_perturbed_att, sink_rules, source_gene, 1, sink_order)
else:
    sink_order = []   # empty — KO necessity loop becomes a no-op, no topo-sort overhead
    print(f"  Sink recovery skipped ({len(sink_nodes)} sinks > {SINK_RECOVERY_THRESHOLD}) — sinks excluded from classification")
    dark_resting_full   = [{**a, source_gene: 0} for a in dark_resting_att]
    dark_perturbed_full = [{**a, source_gene: 1} for a in dark_perturbed_att]
    perm_resting_full   = [{**a, source_gene: 0} for a in perm_resting_att]
    perm_perturbed_full = [{**a, source_gene: 1} for a in perm_perturbed_att]

direct_targets = sorted(set(
    g for g, f in bn_resting_dict.items() if re.search(r'\b' + re.escape(source_gene) + r'\b', str(f))
) | {g for g, r in sink_rules.items() if re.search(r'\b' + re.escape(source_gene) + r'\b', str(r))})

def stable_on(g, atts):  return bool(atts) and all(a.get(g, 0) == 1 for a in atts)
def stable_off(g, atts): return bool(atts) and all(a.get(g, 0) == 0 for a in atts)
def tags(g): return "[" + ", ".join([cat(g)] + (["sink"] if g in sink_nodes else [])) + "]"

G = set(bn_dict) | {source_gene}
dark_activated  = sorted({g for g in G if stable_off(g, dark_resting_full)  and stable_on(g,  dark_perturbed_full)} - {source_gene})
dark_suppressed = sorted({g for g in G if stable_on(g,  dark_resting_full)  and stable_off(g, dark_perturbed_full)} - {source_gene})
perm_activated  = sorted({g for g in G if stable_off(g, perm_resting_full)  and stable_on(g,  perm_perturbed_full)} - {source_gene})
perm_suppressed = sorted({g for g in G if stable_on(g,  perm_resting_full)  and stable_off(g, perm_perturbed_full)} - {source_gene})
robust_activated  = sorted(set(dark_activated)  & set(perm_activated))
robust_suppressed = sorted(set(dark_suppressed) & set(perm_suppressed))

conditional_derepressed = sorted(
    (g, sum(1 for a in perm_resting_full if a.get(g, 0) == 0))
    for g in G if g not in perm_activated and g != source_gene
    and any(a.get(g, 0) == 0 for a in perm_resting_full)
    and any(a.get(g, 0) == 1 for a in perm_perturbed_full)
)
conditional_suppressed = sorted(
    (g, sum(1 for a in perm_resting_full if a.get(g, 0) == 1))
    for g in G if g not in perm_suppressed and g != source_gene
    and any(a.get(g, 0) == 1 for a in perm_resting_full)
    and any(a.get(g, 0) == 0 for a in perm_perturbed_full)
)

gstates = {}
for att in perm_perturbed_full:
    for g, v in att.items(): gstates.setdefault(g, set()).add(v)
variable_genes = sorted(g for g, vs in gstates.items() if len(vs) > 1)
# Group genes by their exact ON/OFF pattern across perturbed attractors.
# Genes sharing a pattern co-vary identically; no normalisation so that
# (0,1) and (1,0) remain distinct groups — they move in opposite directions.
decisions = {}
for g in variable_genes:
    p = tuple(a.get(g, 0) for a in perm_perturbed_full)
    if not all(isinstance(v, int) for v in p): continue
    decisions.setdefault(p, []).append(g)

necessary, dispensable = [], []
if perm_activated:
    # Use the same solver as the attractor analysis for methodological consistency.
    # BoNesis finds all reachable attractors; simulation finds one from all-ones.
    nec_method = "BoNesis" if not using_simulation else "synchronous simulation"
    print(f"\n  Necessity test ({len(direct_targets)} direct target(s); method: {nec_method})...")
    t0 = time.perf_counter()
    for candidate in direct_targets:
        # BoNesis path: base on bn_perturbed_dict so source_gene is constant "1",
        # matching the perturbed condition used to define perm_activated.
        # Simulation path: source_gene is fixed externally via locked/val.
        if not using_simulation:
            ko = dict(bn_perturbed_dict); ko[candidate] = "0"
            _bn_ko    = bonesis.BooleanNetwork(ko)
            ko_states = list(_bn_ko.attractors(reachable_from={g: 1 for g in ko}))
        else:
            ko = dict(bn_dict_pruned); ko[candidate] = "0"
            ko_states, _, _ = simulate(ko, all_on, source_gene, 1)
        # Recover sink states for the KO condition. Sinks are absent from the BN
        # and would otherwise remain at their initial value (1), masking any true
        # loss of activation for sink nodes that appear in perm_activated.
        ko_full = []
        for _att in ko_states:
            _ext = dict(_att); _ext[source_gene] = 1; _ext[candidate] = 0
            for _s in sink_order: _ext[_s] = eval_rule_simple(sink_rules[_s], _ext)
            ko_full.append(_ext)
        lost = sorted(g for g in perm_activated if g != candidate and not stable_on(g, ko_full))
        if lost: necessary.append((candidate, lost))
        else:    dispensable.append(candidate)
    print(f"  Necessity: {time.perf_counter()-t0:.3f}s")
else:
    nec_method = "N/A"
    print(f"\n  Necessity test skipped — no genes stably activated in permissive baseline.")

mode = "synchronous simulation" if using_simulation else "BoNesis (complete attractor landscape)"
print(f"\n" + "="*70)
print(f"  RUN METADATA")
print(f"  Network          : filtered_networkL_normalized.csv")
print(f"  Gene             : {source_gene}  [{cat(source_gene)}]")
print(f"  Hops             : {MAX_HOPS}  (downstream)")
print(f"  Attractor solver : {mode}  |  BoNesis timeout: {BONESIS_TIMEOUT}s")
print(f"  Necessity solver : {nec_method}")
print(f"  Update scheme    : synchronous (all genes updated simultaneously per step)")
print(f"  Boolean rules    : activators combined with OR; suppressors combined with AND NOT")
print(f"  Subgraph         : {len(subgraph_nodes)} nodes  |  pruned BN: {n_pruned}  |  sinks: {len(sink_nodes)}")
print(f"  Run time         : {time.strftime('%Y-%m-%d %H:%M:%S')}")
print(f"\n" + "="*70)
print(f"  {source_gene} — {MAX_HOPS} hop(s)  |  mode: {mode}")

print(f"\n  Direct targets: {len(direct_targets)}")
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
    for g, count in conditional_derepressed: print(f"    + ({count}/{len(perm_resting_att)} resting OFF) {g:36s}  {tags(g)}")
    for g, count in conditional_suppressed:  print(f"    - ({count}/{len(perm_resting_att)} resting ON)  {g:36s}  {tags(g)}")

print(f"\n{'='*70}")
print(f"  ROBUST EFFECTS")
print(f"\n  Robustly activated: {len(robust_activated)}")
for g in robust_activated: print(f"    + {g:40s}  {tags(g)}")
print(f"\n  Robustly suppressed: {len(robust_suppressed)}")
for g in robust_suppressed: print(f"    - {g:40s}  {tags(g)}")

if perm_activated:
    print(f"\n{'='*70}")
    print(f"  NECESSITY TEST  (method: {nec_method} from permissive background)")
    print(f"  Definition: a direct target is necessary if knocking it out (fixing to OFF)")
    print(f"  causes at least one permissive-activated gene to lose stable-ON status.")
    print(f"    Necessary  ({len(necessary)}):")
    for g, lost in necessary:
        print(f"    ! {g:38s}  {tags(g)}")
        print(f"        loss of stable activation: {', '.join(lost)}")
    print(f"    Dispensable ({len(dispensable)}):")
    for g in dispensable: print(f"      {g:40s}  {tags(g)}")

print(f"\n{'='*70}")
n_perturbed_att = len(perm_perturbed_full)
print(f"  Context-dependent genes: {len(variable_genes)} gene(s) in {len(decisions)} co-varying pattern(s)")
if variable_genes:
    print(f"  (Genes with different stable values across the {n_perturbed_att} perturbed attractor(s).")
    print(f"   Each group shares an identical ON/OFF pattern across attractors.)")
for i, (pattern, genes) in enumerate(sorted(decisions.items()), 1):
    on_in  = [j + 1 for j, v in enumerate(pattern) if v == 1]
    off_in = [j + 1 for j, v in enumerate(pattern) if v == 0]
    print(f"\n  Pattern {i} — stably ON in attractor(s) {on_in or 'none'}  |  stably OFF in {off_in or 'none'}")
    for g in sorted(genes): print(f"    {g:40s}  {tags(g)}")

print(f"\n{'='*70}")
print(f"  Sink nodes — terminal downstream effectors regulated by the subgraph but not")
print(f"  referenced by any other gene rule ({len(sink_nodes)} node(s)):")
for g in sorted(sink_nodes):
    r = "ON" if stable_on(g, perm_resting_full) else "OFF" if stable_off(g, perm_resting_full) else "var"
    p = "ON" if stable_on(g, perm_perturbed_full) else "OFF" if stable_off(g, perm_perturbed_full) else "var"
    print(f"    {g:40s}  resting: {r:3s}  perturbed: {p:3s}  {tags(g)}")

if community:
    print(f"\n{'='*70}")
    comm_groups = {}
    for g, cid in community.items(): comm_groups.setdefault(cid, []).append(g)
    print(f"  REGULATORY MODULES (SBM) — {len(comm_groups)} modules")
    for cid, members in sorted(comm_groups.items(), key=lambda x: len(x[1]), reverse=True)[:5]:
        print(f"\n  Module {cid}: {len(members)} nodes")
        for g in sorted(members)[:8]: print(f"    {g:40s}  {tags(g)}")
        if len(members) > 8: print(f"    ... and {len(members)-8} more")

print(f"\n  Total: {time.perf_counter()-t_total_start:.2f}s")
print(f"\n{'='*70}")