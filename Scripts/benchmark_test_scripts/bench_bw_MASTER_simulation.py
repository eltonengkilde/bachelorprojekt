#!/opt/anaconda3/envs/bachelor_env/bin/python
import csv, os, re, heapq, time, keyword, threading, contextlib
import graph_tool.all as gt
import bonesis

GENE                    = "MYB46"
MAX_HOPS                = 7
BONESIS_TIMEOUT         = 0   # seconds; set to 0 to force simulation always
MAX_SIM_STEPS           = 1000
PER_HOP_CUTOFF          = 300
SINK_RECOVERY_THRESHOLD = 10000
SCRIPT_NAME = os.path.basename(__file__).replace('.py', '')
OUT_DIR     = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
BASE_DIR    = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.makedirs(OUT_DIR, exist_ok=True)

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

def eval_rule_simple(rule, state):
    try: return 1 if eval(_to_py(rule), {"__builtins__": {}}, {g: bool(v) for g, v in state.items()}) else 0
    except: return 0

# ── ONE-TIME DATA LOAD ─────────────────────────────────────────────────────────
print(f"Loading: {SCRIPT_NAME}")
activators, suppressors, edges_act, edges_sup, all_nodes = {}, {}, [], [], set()
with open(os.path.join(BASE_DIR, "networks_used_by_scripts", "filtered_networkL_normalized.csv"), newline="", encoding="utf-8") as f:
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
print(f"Phase 1: PageRank top hub: {node_list[max(range(g_full.num_vertices()), key=lambda i: pr_map[i])]}, SCCs: {hist.shape[0]}")

if GENE not in node_idx:
    print(f"ERROR: '{GENE}' not found in network. Exiting."); exit(1)
print(f"Gene '{GENE}' confirmed in network.")
print(f"Benchmark (BACKWARDS): hops 1–{MAX_HOPS} | Mode: {'BoNesis (timeout=' + str(BONESIS_TIMEOUT) + 's)' if BONESIS_TIMEOUT > 0 else 'Synchronous simulation'}\n")

# ── BENCHMARK LOOP ─────────────────────────────────────────────────────────────
timings = []
for hops in range(1, MAX_HOPS + 1):
    print(f"  Hops {hops}...", end=" ", flush=True)
    t_start = time.perf_counter()
    target_gene = GENE
    _bonesis_ok = False
    out_file = os.path.join(OUT_DIR, f"{SCRIPT_NAME}_MYB46_hops{hops}.txt")

    try:
        with open(out_file, 'w') as fout, contextlib.redirect_stdout(fout):
            print(f"# Benchmark: {SCRIPT_NAME}.py")
            print(f"# Gene: {GENE}  |  Hops: {hops}  |  Mode: {'BoNesis' if BONESIS_TIMEOUT > 0 else 'Simulation'}")
            print(f"# Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}")

            # Upstream BFS via reversed graph
            g_rev = gt.GraphView(g_full, reversed=True)
            dist  = gt.shortest_distance(g_rev, source=g_full.vertex(node_idx[target_gene]), directed=True)
            subgraph_nodes = {node_list[i] for i in range(g_full.num_vertices()) if dist[i] <= hops}
            sub_node_list  = sorted(subgraph_nodes)
            sub_v_idx      = {name: i for i, name in enumerate(sub_node_list)}
            print(f"\n{'='*70}")
            print(f"PHASE 2 — Upstream: {target_gene}, {hops} hop(s)  |  {len(subgraph_nodes)} nodes, "
                  f"{sum(1 for s,t in edges_act+edges_sup if s in subgraph_nodes and t in subgraph_nodes)} edges")

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
                for m in sorted(scc_m)[:10]: print(f"    {m}")
                if len(scc_m) > 10: print(f"    ... and {len(scc_m)-10} more")
            community = {}
            if len(subgraph_nodes) <= 5000:
                t0 = time.perf_counter(); vb, _ = gt.betweenness(g_sub)
                print(f"\n  Betweenness ({time.perf_counter()-t0:.3f}s) — top 10:")
                for name in heapq.nlargest(10, sub_node_list, key=lambda n: vb[sub_v_idx[n]]):
                    sc = vb[sub_v_idx[name]]
                    if sc > 0: print(f"    {name:40s}  {sc:.6f}")
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
                    print(f"    Module {cid}: {sz:>4} nodes")
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
            sink_nodes = {g for g in bn_dict if g not in referenced and g != target_gene}
            bn_dict_pruned = {g: f for g, f in bn_dict.items() if g not in sink_nodes}
            n_pruned = len(bn_dict_pruned)
            print(f"  Regulated: {len(bn_dict)}  |  pruned: {n_pruned}  |  sinks: {len(sink_nodes)}")
            def stag(gene): return "  [sink]" if gene in sink_nodes else ""

            all_zeros_sub = {g: 0 for g in subgraph_nodes}
            all_ones_sub  = {g: 1 for g in subgraph_nodes}

            def target_on_any(atts): return any(a.get(target_gene, 0) == 1 for a in atts)

            direct_activators  = sorted(g for g in activators.get(target_gene, set()) if g in subgraph_nodes)
            direct_suppressors = sorted(g for g in suppressors.get(target_gene, set()) if g in subgraph_nodes)

            # Step 2: sufficiency test (always simulation)
            baseline_states, _, _ = simulate(bn_dict_pruned, all_zeros_sub)
            baseline_target_on    = target_on_any(baseline_states)
            upstream_genes = sorted(subgraph_nodes - {target_gene})
            print(f"\n  Running sufficiency test on {len(upstream_genes)} upstream gene(s)...")
            t0 = time.perf_counter()
            sufficient_activators = []
            for candidate in upstream_genes:
                states, _, _ = simulate(bn_dict_pruned, all_zeros_sub, candidate, 1)
                if target_on_any(states) and not baseline_target_on:
                    sufficient_activators.append(candidate)
            print(f"  Sufficiency test completed in {time.perf_counter()-t0:.3f}s")

            # Step 3: necessity test (BoNesis or simulation)
            candidates = sorted(set(direct_activators) | set(direct_suppressors))
            necessary_activators, redundant_activators, suppressor_releases = [], [], []
            baseline_all_on_target_on = False

            if BONESIS_TIMEOUT > 0 and candidates:
                print(f"  Attempting BoNesis necessity test ({n_pruned} nodes, timeout: {BONESIS_TIMEOUT}s)...")
                t_bn = time.perf_counter(); _box = []
                def _run_necessity():
                    base_on = target_on_any(list(bonesis.BooleanNetwork(bn_dict_pruned).attractors(
                        reachable_from={g: 1 for g in bn_dict_pruned})))
                    nec, red, rel = [], [], []
                    for cand in candidates:
                        ko = dict(bn_dict_pruned); ko[cand] = "0"
                        ko_t = target_on_any(list(bonesis.BooleanNetwork(ko).attractors(reachable_from={g: 1 for g in ko})))
                        cb = base_on
                        if cand not in bn_dict_pruned:
                            pd = dict(bn_dict_pruned); pd[cand] = "1"
                            cb = target_on_any(list(bonesis.BooleanNetwork(pd).attractors(reachable_from={g: 1 for g in pd})))
                        if cand in direct_activators:
                            if cb and not ko_t: nec.append(cand)
                            elif cb and ko_t:   red.append(cand)
                        if cand in direct_suppressors:
                            if not cb and ko_t: rel.append(cand)
                    _box.append((base_on, nec, red, rel))
                _t = threading.Thread(target=_run_necessity, daemon=True)
                _t.start(); _t.join(BONESIS_TIMEOUT)
                if _box:
                    baseline_all_on_target_on, necessary_activators, redundant_activators, suppressor_releases = _box[0]
                    print(f"  BoNesis necessity: {time.perf_counter()-t_bn:.2f}s")
                    _bonesis_ok = True
                else:
                    print(f"  BoNesis timed out ({BONESIS_TIMEOUT}s) — falling back to simulation for necessity")

            if not _bonesis_ok and candidates:
                t0 = time.perf_counter()
                base_on_sim, _, _ = simulate(bn_dict_pruned, all_ones_sub)
                baseline_all_on_target_on = target_on_any(base_on_sim)
                for cand in candidates:
                    ko = dict(bn_dict_pruned); ko[cand] = "0"
                    ko_states, _, _ = simulate(ko, all_ones_sub, cand, 0)
                    ko_t = target_on_any(ko_states)
                    cb = baseline_all_on_target_on
                    if cand in direct_activators:
                        if cb and not ko_t: necessary_activators.append(cand)
                        elif cb and ko_t:   redundant_activators.append(cand)
                    if cand in direct_suppressors:
                        if not cb and ko_t: suppressor_releases.append(cand)
                print(f"  Necessity test (simulation) completed in {time.perf_counter()-t0:.3f}s")

            mode_str = f"bonesis ({n_pruned} nodes)" if _bonesis_ok else f"synchronous simulation ({n_pruned} nodes)"
            sink_rules = {g: bn_dict[g] for g in sink_nodes}

            print(f"\n{'='*70}")
            print(f"  BACKWARDS ANALYSIS — what regulates '{target_gene}'?")
            print(f"  {target_gene} — {hops} hop(s) upstream  |  mode: {mode_str}")
            print(f"\n  Step 1 — Direct activators of '{target_gene}': {len(direct_activators)}")
            for g in direct_activators: print(f"    + {g:40s}{stag(g)}")
            print(f"\n  Step 1 — Direct suppressors of '{target_gene}': {len(direct_suppressors)}")
            for g in direct_suppressors: print(f"    - {g:40s}{stag(g)}")
            print(f"\n{'='*70}")
            print(f"  Step 2 — Sufficient activators: {len(sufficient_activators)}")
            for g in sufficient_activators: print(f"    + {g:40s}{stag(g)}")
            print(f"\n{'='*70}")
            print(f"  Step 3 — Necessity test")
            print(f"\n  Necessary activators: {len(necessary_activators)}")
            for g in necessary_activators: print(f"    ! {g:40s}  [required]{stag(g)}")
            print(f"\n  Redundant activators: {len(redundant_activators)}")
            for g in redundant_activators: print(f"    {g:40s}  [dispensable]{stag(g)}")
            print(f"\n  Suppressors whose removal de-represses '{target_gene}': {len(suppressor_releases)}")
            for g in suppressor_releases: print(f"    ~ {g:40s}  [KO turns target ON]{stag(g)}")
            if community:
                print(f"\n{'='*70}")
                cg = {}
                for g, cid in community.items(): cg.setdefault(cid, []).append(g)
                print(f"  REGULATORY MODULES (SBM) — {len(cg)} modules")
                for cid, members in sorted(cg.items(), key=lambda x: len(x[1]), reverse=True)[:5]:
                    print(f"\n  Module {cid}: {len(members)} nodes")
                    for g in sorted(members)[:8]: print(f"    {g:40s}{stag(g)}")
                    if len(members) > 8: print(f"    ... and {len(members)-8} more")
            print(f"\n{'='*70}")
            print(f"  Sink nodes: {len(sink_nodes)}")
            if sink_nodes:
                sim_base, _, _ = simulate(bn_dict_pruned, all_ones_sub)
                for g in sorted(sink_nodes):
                    val = eval_rule_simple(sink_rules[g], sim_base[0])
                    print(f"    {g:40s}  permissive state: {'ON' if val else 'OFF'}")
            print(f"\n  Total hop time: {time.perf_counter()-t_start:.2f}s\n{'='*70}")

        elapsed = time.perf_counter() - t_start
        mode_used = "BoNesis" if _bonesis_ok else "Simulation"
        timings.append((hops, elapsed, mode_used))
        print(f"{elapsed:.1f}s [{mode_used}] → {os.path.basename(out_file)}")
        if elapsed > PER_HOP_CUTOFF:
            print(f"  → Exceeded {PER_HOP_CUTOFF}s cutoff. Stopping benchmark."); break

    except Exception as e:
        import traceback
        elapsed = time.perf_counter() - t_start
        print(f"ERROR after {elapsed:.1f}s: {e}")
        timings.append((hops, elapsed, "ERROR"))
        with open(out_file, 'a') as fout: fout.write(f"\n--- ERROR ---\n{traceback.format_exc()}\n")
        break

print(f"\n{'='*70}")
print(f"BENCHMARK SUMMARY: {SCRIPT_NAME}")
print(f"Gene: {GENE}  |  Mode: {'BoNesis (timeout=' + str(BONESIS_TIMEOUT) + 's)' if BONESIS_TIMEOUT > 0 else 'Synchronous simulation'}")
print(f"{'Hops':>6}  {'Time (s)':>10}  {'Mode':>12}  Output file")
print(f"{'-'*6}  {'-'*10}  {'-'*12}  {'-'*45}")
for hops, t, mode in timings:
    print(f"{hops:>6}  {t:>10.2f}  {mode:>12}  {SCRIPT_NAME}_MYB46_hops{hops}.txt")
