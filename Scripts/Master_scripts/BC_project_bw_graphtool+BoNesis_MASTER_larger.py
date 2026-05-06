#!/opt/anaconda3/envs/bachelor_env/bin/python
import csv, os, re, heapq, time, keyword
import graph_tool.all as gt
import bonesis

BASE_DIR             = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SIMULATION_THRESHOLD = 400   # n_pruned <= threshold → BoNesis necessity test; else simulation
MAX_SIM_STEPS        = 1000  # max synchronous simulation steps before non-convergence warning
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

node_list = sorted(all_nodes)
node_idx  = {n: i for i, n in enumerate(node_list)}
t_total_start = time.perf_counter()

t0 = time.perf_counter()
g_full = gt.Graph(directed=True)
g_full.add_vertex(len(node_list))
gt_etype = g_full.new_edge_property("bool")
for s, t in edges_act: gt_etype[g_full.add_edge(node_idx[s], node_idx[t])] = True
for s, t in edges_sup: gt_etype[g_full.add_edge(node_idx[s], node_idx[t])] = False
g_full.ep["etype"] = gt_etype
print(f"\ngraph-tool: {g_full.num_vertices()} vertices, {g_full.num_edges()} edges  ({time.perf_counter()-t0:.2f}s)")

print(f"\n{'='*70}")
print(f"PHASE 1 — nodes: {g_full.num_vertices()}  |  edges: {g_full.num_edges()} ({len(edges_act)} act / {len(edges_sup)} sup)")
t0 = time.perf_counter()
pr_map = gt.pagerank(g_full)
print(f"  PageRank ({time.perf_counter()-t0:.2f}s) — top 10:")
for i in heapq.nlargest(10, range(g_full.num_vertices()), key=lambda i: pr_map[i]):
    print(f"    {node_list[i]:40s}  PR: {pr_map[i]:.6f}  [{cat(node_list[i])}]")
_, hist = gt.label_components(g_full)
print(f"  SCCs: {hist.shape[0]}  |  non-trivial: {int((hist>1).sum())}  |  largest: {int(hist.max())} nodes")

while True:
    target_gene = input("\nTarget gene to activate (what gene do you want to turn ON?): ")
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

# Upstream BFS via shortest distance on the reversed graph
g_rev = gt.GraphView(g_full, reversed=True)
dist  = gt.shortest_distance(g_rev, source=g_full.vertex(node_idx[target_gene]), directed=True)
subgraph_nodes = {node_list[i] for i in range(g_full.num_vertices()) if dist[i] <= MAX_HOPS}
sub_node_list  = sorted(subgraph_nodes)
sub_v_idx      = {name: i for i, name in enumerate(sub_node_list)}

print(f"\n{'='*70}")
print(f"PHASE 2 — Upstream subgraph: {target_gene}, {MAX_HOPS} hop(s)  |  {len(subgraph_nodes)} nodes, "
      f"{sum(1 for s,t in edges_act+edges_sup if s in sub_v_idx and t in sub_v_idx)} edges")
sc = {}
for n in subgraph_nodes: sc[cat(n)] = sc.get(cat(n), 0) + 1
for c, n in sorted(sc.items(), key=lambda x: x[1], reverse=True): print(f"    {n:>4}  {c}")

print(f"\n{'='*70}")
print(f"PHASE 2.5 — Subgraph topology (graph-tool)")
g_sub = gt.Graph(directed=True)
g_sub.add_vertex(len(sub_node_list))
sub_etype = g_sub.new_edge_property("bool")
for s, t in set((s, t) for s, t in edges_act if s in sub_v_idx and t in sub_v_idx):
    sub_etype[g_sub.add_edge(sub_v_idx[s], sub_v_idx[t])] = True
for s, t in set((s, t) for s, t in edges_sup if s in sub_v_idx and t in sub_v_idx):
    sub_etype[g_sub.add_edge(sub_v_idx[s], sub_v_idx[t])] = False
g_sub.ep["etype"] = sub_etype

sub_comp, sub_hist = gt.label_components(g_sub)
large_sub = int((sub_hist > 1).sum())
print(f"  SCCs: {sub_hist.shape[0]}  |  non-trivial: {large_sub}")
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
    for name in heapq.nlargest(10, sub_node_list, key=lambda n: vb[sub_v_idx[n]]):
        score = vb[sub_v_idx[name]]
        if score > 0: print(f"    {name:40s}  {score:.6f}  [{cat(name)}]")
else:
    print(f"  Betweenness skipped ({len(subgraph_nodes)} nodes > 5000)")

if len(subgraph_nodes) <= 2000:
    t0 = time.perf_counter()
    print(f"\n  Running SBM...")
    state = gt.minimize_blockmodel_dl(g_sub)
    b = state.get_blocks()
    for name in sub_node_list: community[name] = int(b[sub_v_idx[name]])
    comm_sizes = {}
    for c in community.values(): comm_sizes[c] = comm_sizes.get(c, 0) + 1
    print(f"  SBM ({time.perf_counter()-t0:.2f}s) — {len(comm_sizes)} modules:")
    for cid, sz in sorted(comm_sizes.items(), key=lambda x: x[1], reverse=True)[:5]:
        cats = {}
        for m in (n for n, c in community.items() if c == cid): cats[cat(m)] = cats.get(cat(m), 0) + 1
        print(f"    Module {cid}: {sz:>4} nodes  ({', '.join(f'{v} {k}' for k,v in sorted(cats.items(), key=lambda x: x[1], reverse=True))})")
else:
    print(f"\n  SBM skipped ({len(subgraph_nodes)} nodes > 2000)")

bn_dict = {}
for tgt in set(activators) | set(suppressors):
    if tgt not in subgraph_nodes: continue
    act = " | ".join(s for s in activators.get(tgt, []) if s in subgraph_nodes)
    sup = " & ".join("!" + s for s in suppressors.get(tgt, []) if s in subgraph_nodes)
    if act and sup: bn_dict[tgt] = "(" + act + ") & " + sup
    elif act: bn_dict[tgt] = act
    elif sup: bn_dict[tgt] = sup

referenced = set()
for rule in bn_dict.values(): referenced |= set(re.findall(r'\b[a-zA-Z_]\w*\b', str(rule)))
sink_nodes = {gene for gene in bn_dict if gene not in referenced and gene != target_gene}
bn_dict_pruned = {g: f for g, f in bn_dict.items() if g not in sink_nodes}
n_pruned = len(bn_dict_pruned)
print(f"  Regulated: {len(bn_dict)}  |  pruned: {n_pruned}  |  sinks: {len(sink_nodes)}")

def tags(g):
    parts = [cat(g)]
    if g in sink_nodes: parts.append("sink")
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

def target_on_any(atts): return any(a.get(target_gene, 0) == 1 for a in atts)

all_zeros_sub = {gene: 0 for gene in subgraph_nodes}
all_ones_sub  = {gene: 1 for gene in subgraph_nodes}

mode_str = (f"bonesis ({n_pruned} nodes)" if n_pruned <= SIMULATION_THRESHOLD
            else f"synchronous simulation ({n_pruned} nodes)")
print(f"  Mode: {mode_str}")

direct_activators  = sorted(g for g in activators.get(target_gene, set()) if g in subgraph_nodes)
direct_suppressors = sorted(g for g in suppressors.get(target_gene, set()) if g in subgraph_nodes)

baseline_states, _, _ = simulate(bn_dict_pruned, all_zeros_sub)
baseline_target_on    = target_on_any(baseline_states)
print(f"\n  Running sufficiency test on {len(direct_activators)} direct activator(s)...")
t0 = time.perf_counter()
sufficient_activators = []
for candidate in direct_activators:
    states, _, _ = simulate(bn_dict_pruned, all_zeros_sub, candidate, 1)
    if target_on_any(states) and not baseline_target_on:
        sufficient_activators.append(candidate)
print(f"  Sufficiency test completed in {time.perf_counter()-t0:.3f}s")

if n_pruned <= SIMULATION_THRESHOLD:
    bn_full    = bonesis.BooleanNetwork(bn_dict_pruned)
    all_ones_p = {gene: 1 for gene in bn_dict_pruned}
    baseline_all_on_target_on = target_on_any(list(bn_full.attractors(reachable_from=all_ones_p)))
    using_bonesis = True
else:
    baseline_states_on, _, _ = simulate(bn_dict_pruned, all_ones_sub)
    baseline_all_on_target_on = target_on_any(baseline_states_on)
    using_bonesis = False

necessary_activators, redundant_activators, suppressor_releases = [], [], []
candidates = sorted(set(direct_activators) | set(direct_suppressors))
print(f"\n  Running necessity test on {len(candidates)} direct regulator(s)...")
t0 = time.perf_counter()
for candidate in candidates:
    ko = dict(bn_dict_pruned); ko[candidate] = "0"
    if using_bonesis:
        bn_ko    = bonesis.BooleanNetwork(ko)
        ko_start = {gene: 1 for gene in ko}
        ko_target = target_on_any(list(bn_ko.attractors(reachable_from=ko_start)))
        if candidate not in bn_dict_pruned:
            pd = dict(bn_dict_pruned); pd[candidate] = "1"
            cand_base = target_on_any(list(bonesis.BooleanNetwork(pd).attractors(
                reachable_from={g: 1 for g in pd})))
        else:
            cand_base = baseline_all_on_target_on
    else:
        ko_states, _, _ = simulate(ko, all_ones_sub, candidate, 0)
        ko_target = target_on_any(ko_states)
        cand_base = baseline_all_on_target_on
    if candidate in direct_activators:
        if cand_base and not ko_target: necessary_activators.append(candidate)
        elif cand_base and ko_target:   redundant_activators.append(candidate)
    if candidate in direct_suppressors:
        if not cand_base and ko_target: suppressor_releases.append(candidate)
print(f"  Necessity test completed in {time.perf_counter()-t0:.3f}s")

def eval_rule_simple(rule, state):
    try: return 1 if eval(_to_py(rule), {"__builtins__": {}}, {g: bool(v) for g, v in state.items()}) else 0
    except: return 0

sink_rules = {g: bn_dict[g] for g in sink_nodes}

print("\n" + "="*70)
print(f"  RUN METADATA")
print(f"  Network  : filtered_networkL_normalized.csv")
print(f"  Gene     : {target_gene}  [{cat(target_gene)}]")
print(f"  Hops     : {MAX_HOPS}  (upstream)")
print(f"  Solver   : {mode_str}")
print(f"  Subgraph : {len(subgraph_nodes)} nodes  |  pruned: {n_pruned}  |  sinks: {len(sink_nodes)}")
print(f"  Run time : {time.strftime('%Y-%m-%d %H:%M:%S')}")
print("\n" + "="*70)
print(f"  BACKWARDS ANALYSIS — what regulates '{target_gene}'?")
print(f"  {target_gene} — {MAX_HOPS} hop(s) upstream  |  mode: {mode_str}")

print(f"\n  Step 1 — Direct activators of '{target_gene}': {len(direct_activators)}")
for gene in direct_activators: print(f"    + {gene:40s}  {tags(gene)}")
print(f"\n  Step 1 — Direct suppressors of '{target_gene}': {len(direct_suppressors)}")
for gene in direct_suppressors: print(f"    - {gene:40s}  {tags(gene)}")

print(f"\n{'='*70}")
print(f"  Step 2 — Sufficient activators (sufficient alone to activate {target_gene} from silent network): {len(sufficient_activators)}")
for gene in sufficient_activators: print(f"    + {gene:40s}  {tags(gene)}")

print(f"\n{'='*70}")
print(f"  Step 3 — Necessity test (knockout from permissive background; method: {'BoNesis' if using_bonesis else 'synchronous simulation'})")
print(f"\n  Necessary activators (knockout turns '{target_gene}' OFF): {len(necessary_activators)}")
for gene in necessary_activators: print(f"    ! {gene:40s}  [required]  {tags(gene)}")
print(f"\n  Redundant activators ('{target_gene}' stays ON without them): {len(redundant_activators)}")
for gene in redundant_activators: print(f"    {gene:40s}  [dispensable]  {tags(gene)}")
print(f"\n  Suppressors whose removal de-represses '{target_gene}': {len(suppressor_releases)}")
for gene in suppressor_releases: print(f"    ~ {gene:40s}  [KO turns target ON]  {tags(gene)}")

if community:
    print(f"\n{'='*70}")
    comm_groups = {}
    for gene, cid in community.items(): comm_groups.setdefault(cid, []).append(gene)
    print(f"  REGULATORY MODULES (SBM) — {len(comm_groups)} modules")
    for cid, members in sorted(comm_groups.items(), key=lambda x: len(x[1]), reverse=True)[:5]:
        print(f"\n  Module {cid}: {len(members)} nodes")
        for gene in sorted(members)[:8]: print(f"    {gene:40s}  {tags(gene)}")
        if len(members) > 8: print(f"    ... and {len(members)-8} more")

print(f"\n{'='*70}")
print(f"  Sink nodes — upstream branches that do not reach '{target_gene}': {len(sink_nodes)}")
if sink_nodes:
    sim_base, _, _ = simulate(bn_dict_pruned, all_ones_sub)
    for gene in sorted(sink_nodes):
        val = eval_rule_simple(sink_rules[gene], sim_base[0])
        print(f"    {gene:40s}  permissive state: {'ON' if val else 'OFF':3s}  {tags(gene)}")

print(f"\n  Total analysis time: {time.perf_counter() - t_total_start:.2f}s")
print(f"\n{'='*70}")
