#!/opt/anaconda3/envs/bachelor_env/bin/python
import csv, os, re, heapq, time, keyword, contextlib, threading
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

def eval_rule_simple(rule, state):
    try: return 1 if eval(_to_py(rule), {"__builtins__": {}}, {g: bool(v) for g, v in state.items()}) else 0
    except: return 0

def _topo_sort_sinks(sink_rules):
    """Topological evaluation order for sink genes respecting intra-sink dependencies."""
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
    """Evaluate sink gene states in topological order to respect intra-sink dependencies."""
    result = []
    for att in attractors:
        ext = {**att, src: src_val}
        for sink in sink_order:
            ext[sink] = eval_rule_simple(sink_rules[sink], ext)
        result.append(ext)
    return result

# ── ONE-TIME DATA LOAD ─────────────────────────────────────────────────────────
print(f"Loading: {SCRIPT_NAME}")
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

node_list = sorted(all_nodes)
node_idx  = {n: i for i, n in enumerate(node_list)}

t0 = time.perf_counter()
g_full = gt.Graph(directed=True)
g_full.add_vertex(len(node_list))
gt_etype = g_full.new_edge_property("bool")
for s, t in edges_act: gt_etype[g_full.add_edge(node_idx[s], node_idx[t])] = True
for s, t in edges_sup: gt_etype[g_full.add_edge(node_idx[s], node_idx[t])] = False
g_full.ep["etype"] = gt_etype
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
            print(f"# Gene: {GENE}  |  Hops: {hops}  |  Mode: BoNesis")
            print(f"# Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}")

            dist = gt.shortest_distance(g_full, source=g_full.vertex(node_idx[source_gene]), directed=True)
            subgraph_nodes = {node_list[i] for i in range(g_full.num_vertices()) if dist[i] <= hops}
            sub_node_list  = sorted(subgraph_nodes)
            sub_v_idx      = {name: i for i, name in enumerate(sub_node_list)}
            print(f"  [time] Subgraph extraction: {time.perf_counter()-t_start:.3f}s")
            sc = {}
            for n in subgraph_nodes: sc[cat(n)] = sc.get(cat(n), 0) + 1
            print(f"\n{'='*70}")
            print(f"PHASE 2 — {source_gene}, {hops} hop(s)  |  {len(subgraph_nodes)} nodes, "
                  f"{sum(1 for s,t in edges_act+edges_sup if s in subgraph_nodes and t in subgraph_nodes)} edges")
            for c, n in sorted(sc.items(), key=lambda x: x[1], reverse=True): print(f"    {n:>4}  {c}")

            print(f"\n{'='*70}\nPHASE 2.5 — Subgraph topology")
            g_sub = gt.Graph(directed=True)
            g_sub.add_vertex(len(sub_node_list))
            sub_etype = g_sub.new_edge_property("bool")
            for s, t in set((s,t) for s,t in edges_act if s in sub_v_idx and t in sub_v_idx):
                sub_etype[g_sub.add_edge(sub_v_idx[s], sub_v_idx[t])] = True
            for s, t in set((s,t) for s,t in edges_sup if s in sub_v_idx and t in sub_v_idx):
                sub_etype[g_sub.add_edge(sub_v_idx[s], sub_v_idx[t])] = False
            g_sub.ep["etype"] = sub_etype
            sub_comp, sub_hist = gt.label_components(g_sub)
            large_sub = int((sub_hist > 1).sum())
            print(f"  SCCs: {sub_hist.shape[0]}  |  non-trivial: {large_sub}")
            if large_sub:
                scc_m = [sub_node_list[i] for i in range(g_sub.num_vertices()) if int(sub_comp[i]) == int(sub_hist.argmax())]
                print(f"  Largest SCC ({int(sub_hist.max())} nodes):")
                for m in sorted(scc_m)[:10]: print(f"    {m:40s}  [{cat(m)}]")
                if len(scc_m) > 10: print(f"    ... and {len(scc_m)-10} more")
            community = {}
            if len(subgraph_nodes) <= 5000:
                t0 = time.perf_counter(); vb, _ = gt.betweenness(g_sub)
                print(f"\n  Betweenness ({time.perf_counter()-t0:.3f}s) — top 10:")
                for name in heapq.nlargest(10, sub_node_list, key=lambda n: vb[sub_v_idx[n]]):
                    sc2 = vb[sub_v_idx[name]]
                    if sc2 > 0: print(f"    {name:40s}  {sc2:.6f}  [{cat(name)}]")
            else:
                print(f"  Betweenness skipped ({len(subgraph_nodes)} nodes > 5000)")
            if len(subgraph_nodes) <= 2000:
                t0 = time.perf_counter(); print(f"\n  Running SBM...")
                sbm_s = gt.minimize_blockmodel_dl(g_sub); b = sbm_s.get_blocks()
                for name in sub_node_list: community[name] = int(b[sub_v_idx[name]])
                cs = {}
                for c in community.values(): cs[c] = cs.get(c, 0) + 1
                print(f"  SBM ({time.perf_counter()-t0:.2f}s) — {len(cs)} modules:")
                for cid, sz in sorted(cs.items(), key=lambda x: x[1], reverse=True)[:5]:
                    mems = [n for n, c in community.items() if c == cid]
                    cats = {}
                    for m in mems: cats[cat(m)] = cats.get(cat(m), 0) + 1
                    print(f"    Module {cid}: {sz:>4} nodes  ({', '.join(f'{v} {k}' for k,v in sorted(cats.items(), key=lambda x: x[1], reverse=True))})")
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
            sink_nodes = {g for g in bn_dict if g not in referenced and g != source_gene}
            bn_dict_pruned = {g: f for g, f in bn_dict.items() if g not in sink_nodes}
            n_pruned = len(bn_dict_pruned)
            print(f"  Regulated: {len(bn_dict)}  |  pruned: {n_pruned}  |  sinks: {len(sink_nodes)}")
            def tags(g): return "[" + ", ".join([cat(g)] + (["sink"] if g in sink_nodes else [])) + "]"

            bn_resting_dict   = dict(bn_dict_pruned); bn_resting_dict[source_gene]   = "0"
            bn_perturbed_dict = dict(bn_dict_pruned); bn_perturbed_dict[source_gene] = "1"
            all_on  = {g: 1 for g in subgraph_nodes}

            dk = {g: 0 for g in bn_dict_pruned}
            sh = {g: 1 for g in bn_dict_pruned}
            _remaining = max(0, MAX_TOTAL_SECONDS - (time.perf_counter() - t_benchmark_start))
            print(f"  Running BoNesis ({n_pruned} nodes — budget remaining: {_remaining:.0f}s)...")
            t_bn = time.perf_counter()
            _bn_r = bonesis.BooleanNetwork(bn_resting_dict)
            _bn_p = bonesis.BooleanNetwork(bn_perturbed_dict)
            _box = []
            _t = threading.Thread(daemon=True, target=lambda: _box.append((
                list(_bn_r.attractors(reachable_from=dk)),
                list(_bn_p.attractors(reachable_from=dk)),
                list(_bn_r.attractors(reachable_from=sh)),
                list(_bn_p.attractors(reachable_from=sh)),
            )))
            _t.start(); _t.join(_remaining)
            if _box:
                dark_resting_att, dark_perturbed_att, perm_resting_att, perm_perturbed_att = _box[0]
                print(f"  [time] BoNesis attractors: {time.perf_counter()-t_bn:.3f}s")
            else:
                print(f"  BoNesis timed out — budget exhausted after {time.perf_counter()-t_bn:.0f}s")
                raise _BudgetExhausted()

            t_sink = time.perf_counter()
            sink_rules = {g: bn_dict[g] for g in sink_nodes}
            if len(sink_nodes) <= SINK_RECOVERY_THRESHOLD:
                sink_order = _topo_sort_sinks(sink_rules)
                dark_resting_full   = recover_sinks(dark_resting_att,   sink_rules, source_gene, 0, sink_order)
                dark_perturbed_full = recover_sinks(dark_perturbed_att, sink_rules, source_gene, 1, sink_order)
                perm_resting_full   = recover_sinks(perm_resting_att,   sink_rules, source_gene, 0, sink_order)
                perm_perturbed_full = recover_sinks(perm_perturbed_att, sink_rules, source_gene, 1, sink_order)
            else:
                sink_order = []   # empty — KO necessity loop becomes a no-op, no topo-sort overhead
                print(f"  Sink recovery skipped ({len(sink_nodes)} sinks > {SINK_RECOVERY_THRESHOLD})")
                dark_resting_full   = [{**a, source_gene: 0} for a in dark_resting_att]
                dark_perturbed_full = [{**a, source_gene: 1} for a in dark_perturbed_att]
                perm_resting_full   = [{**a, source_gene: 0} for a in perm_resting_att]
                perm_perturbed_full = [{**a, source_gene: 1} for a in perm_perturbed_att]
            print(f"  [time] Sink recovery: {time.perf_counter()-t_sink:.3f}s")

            direct_targets = sorted(set(
                g for g, f in bn_resting_dict.items() if re.search(r'\b' + re.escape(source_gene) + r'\b', str(f))
            ) | {g for g, r in sink_rules.items() if re.search(r'\b' + re.escape(source_gene) + r'\b', str(r))})

            def stable_on(g, atts):  return bool(atts) and all(a.get(g, 0) == 1 for a in atts)
            def stable_off(g, atts): return bool(atts) and all(a.get(g, 0) == 0 for a in atts)

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
                and any(a.get(g, 0) == 1 for a in perm_perturbed_full))
            conditional_suppressed = sorted(
                (g, sum(1 for a in perm_resting_full if a.get(g, 0) == 1))
                for g in G if g not in perm_suppressed and g != source_gene
                and any(a.get(g, 0) == 1 for a in perm_resting_full)
                and any(a.get(g, 0) == 0 for a in perm_perturbed_full))

            gstates = {}
            for att in perm_perturbed_full:
                for g, v in att.items(): gstates.setdefault(g, set()).add(v)
            variable_genes = sorted(g for g, vs in gstates.items() if len(vs) > 1)
            decisions = {}
            for g in variable_genes:
                p = tuple(a.get(g, 0) for a in perm_perturbed_full)
                if not all(isinstance(v, int) for v in p): continue
                decisions.setdefault(p, []).append(g)

            necessary, dispensable = [], []
            if perm_activated:
                print(f"\n  Running necessity test on {len(direct_targets)} direct target(s) (BoNesis)...")
                t0 = time.perf_counter()
                for candidate in direct_targets:
                    if MAX_TOTAL_SECONDS - (time.perf_counter() - t_benchmark_start) <= 0:
                        print(f"  Necessity test budget exhausted — results partial"); break
                    # Use bn_perturbed_dict so source_gene stays constant "1",
                    # matching the perturbed condition that defined perm_activated.
                    ko = dict(bn_perturbed_dict); ko[candidate] = "0"
                    _bn_ko    = bonesis.BooleanNetwork(ko)
                    ko_states = list(_bn_ko.attractors(reachable_from={g: 1 for g in ko}))
                    # Recover sink states so sink nodes in perm_activated are
                    # evaluated correctly rather than staying at initial value 1.
                    ko_full = []
                    for _att in ko_states:
                        _ext = dict(_att); _ext[source_gene] = 1; _ext[candidate] = 0
                        for _s in sink_order: _ext[_s] = eval_rule_simple(sink_rules[_s], _ext)
                        ko_full.append(_ext)
                    lost = sorted(g for g in perm_activated if g != candidate
                                  and not stable_on(g, ko_full))
                    if lost: necessary.append((candidate, lost))
                    else:    dispensable.append(candidate)
                print(f"  Necessity test completed in {time.perf_counter()-t0:.3f}s")
            else:
                print(f"\n  Necessity test skipped — no genes stably activated in permissive baseline.")

            print(f"\n{'='*70}")
            print(f"  {source_gene} — {hops} hop(s)  |  BoNesis (complete attractor landscape)")
            print(f"  Update: synchronous  |  Rules: activators OR'd, suppressors AND NOT'd")
            print(f"  Necessity solver: BoNesis")
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
                for g, count in conditional_derepressed:
                    print(f"    + ({count}/{len(perm_resting_att)} resting OFF) {g:36s}  {tags(g)}")
                for g, count in conditional_suppressed:
                    print(f"    - ({count}/{len(perm_resting_att)} resting ON)  {g:36s}  {tags(g)}")
            print(f"\n{'='*70}\n  ROBUST EFFECTS")
            print(f"\n  Robustly activated: {len(robust_activated)}")
            for g in robust_activated: print(f"    + {g:40s}  {tags(g)}")
            print(f"\n  Robustly suppressed: {len(robust_suppressed)}")
            for g in robust_suppressed: print(f"    - {g:40s}  {tags(g)}")
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
            n_perturbed_att = len(perm_perturbed_full)
            print(f"  Context-dependent genes: {len(variable_genes)} gene(s) in {len(decisions)} co-varying pattern(s)")
            if variable_genes:
                print(f"  (Different stable values across {n_perturbed_att} perturbed attractor(s); groups share identical ON/OFF patterns.)")
            for i, (pattern, genes) in enumerate(sorted(decisions.items()), 1):
                on_in  = [j + 1 for j, v in enumerate(pattern) if v == 1]
                off_in = [j + 1 for j, v in enumerate(pattern) if v == 0]
                print(f"\n  Pattern {i} — stably ON in attractor(s) {on_in or 'none'}  |  stably OFF in {off_in or 'none'}")
                for g in sorted(genes): print(f"    {g:40s}  {tags(g)}")
            print(f"\n{'='*70}")
            print(f"  Sink nodes — terminal downstream effectors ({len(sink_nodes)} node(s)):")
            for g in sorted(sink_nodes):
                r = "ON" if stable_on(g, perm_resting_full) else "OFF" if stable_off(g, perm_resting_full) else "variable"
                p = "ON" if stable_on(g, perm_perturbed_full) else "OFF" if stable_off(g, perm_perturbed_full) else "variable"
                print(f"    {g:40s}  resting: {r:8s}  perturbed: {p:8s}  {tags(g)}")
            if community:
                print(f"\n{'='*70}")
                cg = {}
                for g, cid in community.items(): cg.setdefault(cid, []).append(g)
                print(f"  REGULATORY MODULES (SBM) — {len(cg)} modules")
                for cid, members in sorted(cg.items(), key=lambda x: len(x[1]), reverse=True)[:5]:
                    print(f"\n  Module {cid}: {len(members)} nodes")
                    for g in sorted(members)[:8]: print(f"    {g:40s}  {tags(g)}")
                    if len(members) > 8: print(f"    ... and {len(members)-8} more")
            print(f"\n======================================================================")
            print(f"  TIMING SUMMARY — {source_gene}, hop {hops}")
            print(f"    Total hop time: {time.perf_counter()-t_start:.3f}s  |  step times logged above")
            print(f"======================================================================")

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
