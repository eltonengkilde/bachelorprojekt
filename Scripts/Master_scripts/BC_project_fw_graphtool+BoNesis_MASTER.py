#!/opt/anaconda3/envs/bachelor_env/bin/python
import csv, os, re, heapq, time, keyword
import graph_tool.all as gt
import bonesis

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BOOLEAN_RESERVED = {"TRUE", "FALSE", "NOT", "AND", "OR"}
ACT_REL = "Activation / Induction / Causation / Result"
SUP_REL = "Repression / Inhibition / Negative Regulation"

def clean_name(name):
    if not name or not isinstance(name, str): return "unknown"
    name = re.sub(r"_+", "_", re.sub(r"[^a-zA-Z0-9_]", "_", name)).strip("_")
    if not name: return "unknown"
    if name[0].isdigit(): return "g" + name
    if name.upper() in BOOLEAN_RESERVED or keyword.iskeyword(name): return "gene_" + name
    return name

def _to_py(r): return r.replace("!", " not ").replace("|", " or ").replace("&", " and ")

activators, suppressors, edges_act, edges_sup, all_nodes = {}, {}, [], [], set()
with open(os.path.join(BASE_DIR, "networks_used_by_scripts", "filtered_large_normalized.csv"), newline="", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        rel = row["relationship_category"]
        if rel not in (ACT_REL, SUP_REL): continue
        s, t = clean_name(row["source"]), clean_name(row["target"])
        all_nodes.update((s, t))
        if rel == ACT_REL: activators.setdefault(t, set()).add(s); edges_act.append((s, t))
        else:              suppressors.setdefault(t, set()).add(s); edges_sup.append((s, t))

print(f"Loaded {len(edges_act)+len(edges_sup)} edges  ({len(all_nodes)} nodes)")

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
out_deg = g_full.get_out_degrees(range(g_full.num_vertices()))
in_deg  = g_full.get_in_degrees(range(g_full.num_vertices()))
print(f"  Out-degree — max {int(out_deg.max())}, mean {out_deg.mean():.1f}")
print(f"  In-degree  — max {int(in_deg.max())},  mean {in_deg.mean():.1f}")
print(f"\n  Top 10 hubs by out-degree:")
for deg, name in heapq.nlargest(10, zip(out_deg.tolist(), node_list)):
    print(f"    {name:40s}  out-degree: {int(deg)}")
t0 = time.perf_counter()
pr_map = gt.pagerank(g_full)
print(f"\n  PageRank ({time.perf_counter()-t0:.2f}s) — top 10:")
for i in heapq.nlargest(10, range(g_full.num_vertices()), key=lambda i: pr_map[i]):
    print(f"    {node_list[i]:40s}  PageRank: {pr_map[i]:.4f}")
_, hist = gt.label_components(g_full)
print(f"\n  SCCs: {hist.shape[0]}  |  largest: {int(hist.max())} nodes")
print(f"  Phase 1 completed in {time.perf_counter() - t_total_start:.2f}s")

while True:
    source_gene = input("\nSource gene to investigate: ")
    if source_gene not in node_idx:
        matches = [n for n in node_list if source_gene.upper() in n.upper()]
        print(f"WARNING: '{source_gene}' not found. Try again.")
        if matches: print(f"  Did you mean one of these? {matches[:5]}")
    else:
        print(f"OK: '{source_gene}' found in network"); break

while True:
    MAX_HOPS = int(input("\nNumber of hops to explore downstream: "))
    if MAX_HOPS < 1: print("Must be a whole number larger than 0")
    else: print(f"OK: {MAX_HOPS} hops downstream from '{source_gene}'"); break

if input('\nType and enter "yolo" to run the analysis: ') != "yolo":
    print("Invalid input, exiting."); exit()

dist = gt.shortest_distance(g_full, source=g_full.vertex(node_idx[source_gene]), directed=True)
subgraph_nodes = {node_list[i] for i in range(g_full.num_vertices()) if dist[i] <= MAX_HOPS}
sub_node_list  = sorted(subgraph_nodes)
sub_v_idx      = {name: i for i, name in enumerate(sub_node_list)}

print(f"\n{'='*70}")
print(f"PHASE 2 — {source_gene}, {MAX_HOPS} hop(s)  |  {len(subgraph_nodes)} nodes, "
      f"{sum(1 for s,t in edges_act+edges_sup if s in sub_v_idx and t in sub_v_idx)} edges")

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
    for m in sorted(scc_members)[:10]: print(f"    {m}")
    if len(scc_members) > 10: print(f"    ... and {len(scc_members)-10} more")

community = {}
if len(subgraph_nodes) <= 5000:
    t0 = time.perf_counter()
    vb, _ = gt.betweenness(g_sub)
    print(f"\n  Betweenness ({time.perf_counter()-t0:.3f}s) — top 10:")
    for name in heapq.nlargest(10, sub_node_list, key=lambda n: vb[sub_v_idx[n]]):
        score = vb[sub_v_idx[name]]
        if score > 0: print(f"    {name:40s}  {score:.6f}")
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
        print(f"    Module {cid}: {sz:>4} nodes")
else:
    print(f"\n  SBM skipped ({len(subgraph_nodes)} nodes > 2000)")

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
sink_nodes = {gene for gene in bn_dict if gene not in referenced and gene != source_gene}
bn_dict_pruned = {g: f for g, f in bn_dict.items() if g not in sink_nodes}
n_pruned = len(bn_dict_pruned)
print(f"  Regulated: {len(bn_dict)}  |  pruned: {n_pruned}  |  sinks: {len(sink_nodes)}")

SIMULATION_THRESHOLD = 600
MAX_SIM_STEPS        = 1000
bn_resting_dict   = dict(bn_dict_pruned); bn_resting_dict[source_gene]   = "0"
bn_perturbed_dict = dict(bn_dict_pruned); bn_perturbed_dict[source_gene] = "1"
all_off = {gene: 0 for gene in subgraph_nodes}
all_on  = {gene: 1 for gene in subgraph_nodes}

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

if n_pruned <= SIMULATION_THRESHOLD:
    bn_resting   = bonesis.BooleanNetwork(bn_resting_dict)
    bn_perturbed = bonesis.BooleanNetwork(bn_perturbed_dict)
    dark_start   = {gene: 0 for gene in bn_dict_pruned}
    shared_start = {gene: 1 for gene in bn_dict_pruned}
    print(f"  Mode: bonesis ({n_pruned} nodes)")
    t0 = time.perf_counter()
    dark_resting_att   = list(bn_resting.attractors(reachable_from=dark_start))
    dark_perturbed_att = list(bn_perturbed.attractors(reachable_from=dark_start))
    perm_resting_att   = list(bn_resting.attractors(reachable_from=shared_start))
    perm_perturbed_att = list(bn_perturbed.attractors(reachable_from=shared_start))
    print(f"  Attractors: {time.perf_counter()-t0:.2f}s")
    using_simulation = False
else:
    print(f"  Mode: synchronous simulation ({n_pruned} nodes)")
    t0 = time.perf_counter()
    dr_states, _, dr_conv = simulate(bn_dict_pruned, all_off, source_gene, 0)
    dp_states, _, dp_conv = simulate(bn_dict_pruned, all_off, source_gene, 1)
    pr_states, _, pr_conv = simulate(bn_dict_pruned, all_on,  source_gene, 0)
    pp_states, _, pp_conv = simulate(bn_dict_pruned, all_on,  source_gene, 1)
    print(f"  dark rest: {_sim_label(dr_states, dr_conv)}  pert: {_sim_label(dp_states, dp_conv)}")
    print(f"  perm rest: {_sim_label(pr_states, pr_conv)}  pert: {_sim_label(pp_states, pp_conv)}")
    print(f"  Simulation: {time.perf_counter()-t0:.3f}s")
    if not all([dr_conv, dp_conv, pr_conv, pp_conv]):
        print(f"  WARNING: not all runs converged within {MAX_SIM_STEPS} steps")
    print(f"  NOTE: simulation finds ONE attractor per condition — use bonesis (fewer hops) for full landscape.")
    dark_resting_att, dark_perturbed_att, perm_resting_att, perm_perturbed_att = dr_states, dp_states, pr_states, pp_states
    using_simulation = True

def eval_rule_simple(rule, state):
    try: return 1 if eval(_to_py(rule), {"__builtins__": {}}, {g: bool(v) for g, v in state.items()}) else 0
    except: return 0

def recover_sinks(attractors, sink_rules, src, src_val):
    result = []
    for att in attractors:
        ext = {**att, src: src_val}
        for sink, rule in sink_rules.items(): ext[sink] = eval_rule_simple(rule, ext)
        result.append(ext)
    return result

sink_rules = {g: bn_dict[g] for g in sink_nodes}
dark_resting_full   = recover_sinks(dark_resting_att,   sink_rules, source_gene, 0)
dark_perturbed_full = recover_sinks(dark_perturbed_att, sink_rules, source_gene, 1)
perm_resting_full   = recover_sinks(perm_resting_att,   sink_rules, source_gene, 0)
perm_perturbed_full = recover_sinks(perm_perturbed_att, sink_rules, source_gene, 1)

direct_targets = sorted(set(
    gene for gene, func in bn_resting_dict.items()
    if re.search(r'\b' + re.escape(source_gene) + r'\b', str(func))
) | {gene for gene, rule in sink_rules.items()
     if re.search(r'\b' + re.escape(source_gene) + r'\b', str(rule))})

def stable_on(gene, atts):  return all(a.get(gene, 0) == 1 for a in atts)
def stable_off(gene, atts): return all(a.get(gene, 0) == 0 for a in atts)
def stag(gene): return "  [sink]" if gene in sink_nodes else ""

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
    for gene, val in att.items(): gstates.setdefault(gene, set()).add(val)
variable_genes = sorted(g for g, vs in gstates.items() if len(vs) > 1)
decisions = {}
for g in variable_genes:
    p = tuple(a.get(g, 0) for a in perm_perturbed_full)
    if not all(isinstance(v, int) for v in p): continue
    decisions.setdefault(p if p[0] == 0 else tuple(1-v for v in p), []).append(g)

necessary, dispensable = [], []
if perm_activated:
    print(f"\n  Running necessity test on {len(direct_targets)} direct target(s)...")
    t0 = time.perf_counter()
    for candidate in direct_targets:
        ko = dict(bn_dict_pruned); ko[candidate] = "0"
        ko_states, _, _ = simulate(ko, all_on, source_gene, 1)
        lost = sorted(g for g in perm_activated if g != candidate
                      and not all(s.get(g, 0) == 1 for s in ko_states))
        if lost: necessary.append((candidate, lost))
        else:    dispensable.append(candidate)
    print(f"  Necessity test completed in {time.perf_counter()-t0:.3f}s")
else:
    print(f"\n  Necessity test skipped — no genes stably activated in permissive baseline.")

mode = "synchronous simulation" if using_simulation else "bonesis attractors"
print(f"\n{'='*70}")
print(f"  {source_gene} — {MAX_HOPS} hop(s)  |  mode: {mode}")

print(f"\n  Direct targets: {len(direct_targets)}")
for gene in direct_targets: print(f"    {gene:40s}{stag(gene)}")

print(f"\n{'='*70}")
print(f"  EXP A — dark  |  {len(dark_resting_att)} resting / {len(dark_perturbed_att)} perturbed")
print(f"\n  Activated (OFF→ON): {len(dark_activated)}")
for gene in dark_activated: print(f"    + {gene:40s}{stag(gene)}")
print(f"\n  Suppressed (ON→OFF): {len(dark_suppressed)}")
for gene in dark_suppressed: print(f"    - {gene:40s}{stag(gene)}")

print(f"\n{'='*70}")
print(f"  EXP B — permissive  |  {len(perm_resting_att)} resting / {len(perm_perturbed_att)} perturbed")
print(f"\n  Activated (OFF→ON): {len(perm_activated)}")
for gene in perm_activated: print(f"    + {gene:40s}{stag(gene)}")
print(f"\n  Suppressed (ON→OFF): {len(perm_suppressed)}")
for gene in perm_suppressed: print(f"    - {gene:40s}{stag(gene)}")
if conditional_derepressed or conditional_suppressed:
    print(f"\n  Conditional effects:")
    for gene, count in conditional_derepressed:
        print(f"    + ({count}/{len(perm_resting_att)} resting OFF) {gene:36s}{stag(gene)}")
    for gene, count in conditional_suppressed:
        print(f"    - ({count}/{len(perm_resting_att)} resting ON)  {gene:36s}{stag(gene)}")

print(f"\n{'='*70}")
print(f"  ROBUST EFFECTS")
print(f"\n  Robustly activated: {len(robust_activated)}")
for gene in robust_activated: print(f"    + {gene:40s}{stag(gene)}")
print(f"\n  Robustly suppressed: {len(robust_suppressed)}")
for gene in robust_suppressed: print(f"    - {gene:40s}{stag(gene)}")

if perm_activated:
    print(f"\n{'='*70}")
    print(f"  NECESSITY TEST")
    print(f"    Necessary: {len(necessary)}")
    for gene, lost in necessary:
        print(f"    ! {gene:40s}{stag(gene)}")
        print(f"        causes loss of: {', '.join(lost)}")
    print(f"    Dispensable: {len(dispensable)}")
    for gene in dispensable: print(f"      {gene:40s}{stag(gene)}")

print(f"\n{'='*70}")
print(f"  Context-dependent: {len(variable_genes)} genes, {len(decisions)} decision(s)")
for i, genes in enumerate(sorted(decisions.values(), key=lambda x: x[0]), 1):
    print(f"  Decision {i}:")
    for gene in sorted(genes): print(f"    {gene:40s}{stag(gene)}")

print(f"\n{'='*70}")
print(f"  Sink nodes: {len(sink_nodes)}")
for gene in sorted(sink_nodes):
    r = "ON" if stable_on(gene, perm_resting_full) else "OFF" if stable_off(gene, perm_resting_full) else "variable"
    p = "ON" if stable_on(gene, perm_perturbed_full) else "OFF" if stable_off(gene, perm_perturbed_full) else "variable"
    print(f"    {gene:40s}  resting: {r:8s}  perturbed: {p}")

if community:
    print(f"\n{'='*70}")
    comm_groups = {}
    for gene, cid in community.items(): comm_groups.setdefault(cid, []).append(gene)
    print(f"  REGULATORY MODULES (SBM) — {len(comm_groups)} modules")
    for cid, members in sorted(comm_groups.items(), key=lambda x: len(x[1]), reverse=True)[:5]:
        print(f"\n  Module {cid}: {len(members)} nodes")
        for gene in sorted(members)[:8]: print(f"    {gene:40s}{stag(gene)}")
        if len(members) > 8: print(f"    ... and {len(members)-8} more")

print(f"\n  Total analysis time: {time.perf_counter() - t_total_start:.2f}s")
print(f"\n{'='*70}")