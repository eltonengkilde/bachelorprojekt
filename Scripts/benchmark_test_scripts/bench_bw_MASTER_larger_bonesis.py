#!/opt/anaconda3/envs/bachelor_env/bin/python
import csv, os, re, heapq, time, keyword, contextlib, threading
import graph_tool.all as gt
import bonesis

GENE                    = "MYB46"
MAX_TOTAL_SECONDS       = 1800   # 0.5-hour budget; a new hop only starts if time remains
MAX_HOPS                =  30     # safety cap — time limit is the primary stop condition
MAX_SIM_STEPS           = 1000
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
print(f"Benchmark (BACKWARDS): time limit {MAX_TOTAL_SECONDS}s ({MAX_TOTAL_SECONDS//3600}h) | Mode: BoNesis (budget-limited)\n")

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
    target_gene = GENE
    out_file = os.path.join(OUT_DIR, f"{SCRIPT_NAME}_MYB46_hops{hops}.txt")

    try:
        with open(out_file, 'w') as fout, contextlib.redirect_stdout(fout):
            print(f"# Benchmark: {SCRIPT_NAME}.py")
            print(f"# Gene: {GENE}  |  Hops: {hops}  |  Mode: BoNesis")
            print(f"# Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}")

            g_rev = gt.GraphView(g_full, reversed=True)
            dist  = gt.shortest_distance(g_rev, source=g_full.vertex(node_idx[target_gene]), directed=True)
            subgraph_nodes = {node_list[i] for i in range(g_full.num_vertices()) if dist[i] <= hops}
            sub_node_list  = sorted(subgraph_nodes)
            sub_v_idx      = {name: i for i, name in enumerate(sub_node_list)}
            hop_of = {node_list[i]: int(dist[i]) for i in range(g_full.num_vertices()) if int(dist[i]) <= hops}
            print(f"  [time] Subgraph extraction: {time.perf_counter()-t_start:.3f}s")
            sc = {}
            for n in subgraph_nodes: sc[cat(n)] = sc.get(cat(n), 0) + 1
            print(f"\n{'='*70}")
            print(f"PHASE 2 — Upstream: {target_gene}, {hops} hop(s)  |  {len(subgraph_nodes)} nodes, "
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
            def tags(g):
                parts = [cat(g)]
                if g in sink_nodes: parts.append("sink")
                h = hop_of.get(g)
                if h is not None and h > 1: parts.append(f"hop {h}")
                return "[" + ", ".join(parts) + "]"

            all_zeros_sub = {g: 0 for g in subgraph_nodes}
            all_ones_sub  = {g: 1 for g in subgraph_nodes}

            def target_on_any(atts):  return any(a.get(target_gene, 0) == 1 for a in atts)
            def stable_on(g, atts):  return bool(atts) and all(a.get(g, 0) == 1 for a in atts)
            def stable_off(g, atts): return bool(atts) and all(a.get(g, 0) == 0 for a in atts)

            direct_activators  = sorted(g for g in activators.get(target_gene, set()) if g in subgraph_nodes)
            direct_suppressors = sorted(g for g in suppressors.get(target_gene, set()) if g in subgraph_nodes)
            all_upstream       = sorted(subgraph_nodes - {target_gene})

            # Two baseline simulations (dark / permissive)
            t0 = time.perf_counter()
            att_dark, _, dark_conv = simulate(bn_dict_pruned, all_zeros_sub)
            att_perm, _, perm_conv = simulate(bn_dict_pruned, all_ones_sub)
            print(f"  [time] Baseline simulations: {time.perf_counter()-t0:.3f}s")
            if not dark_conv: print(f"  WARNING: dark attractor did not converge within {MAX_SIM_STEPS} steps")
            if not perm_conv: print(f"  WARNING: perm attractor did not converge within {MAX_SIM_STEPS} steps")
            print(f"  NOTE: simulation finds ONE attractor per starting condition.")
            target_in_dark = target_on_any(att_dark)
            target_in_perm = target_on_any(att_perm)
            _lbl = lambda atts, c: ("fixed point" if len(atts)==1 else f"cycle/{len(atts)}") + (" (conv)" if c else " (max)")
            print(f"  Dark: target={'ON' if target_in_dark else 'OFF'}  {_lbl(att_dark, dark_conv)}")
            print(f"  Perm: target={'ON' if target_in_perm else 'OFF'}  {_lbl(att_perm, perm_conv)}")

            # Hop-by-hop stable-state attractor trace
            all_hops = sorted({hop_of[g] for g in all_upstream if hop_of.get(g) is not None})
            hop_pod = {h: [] for h in all_hops}
            hop_opd = {h: [] for h in all_hops}
            for g in all_upstream:
                h = hop_of.get(g)
                if h is None: continue
                if stable_on(g, att_perm) and stable_off(g, att_dark): hop_pod[h].append(g)
                elif stable_off(g, att_perm) and stable_on(g, att_dark): hop_opd[h].append(g)

            perm_on_stable  = sorted(g for g in all_upstream if stable_on(g,  att_perm))
            perm_off_stable = sorted(g for g in all_upstream if stable_off(g, att_perm))
            necessity_act_candidates = perm_on_stable if target_in_perm else []
            sufficiency_candidates   = sorted(g for g in necessity_act_candidates
                                              if stable_off(g, att_dark)) if (target_in_perm and not target_in_dark) else []
            necessity_sup_candidates = perm_off_stable if target_in_perm else []
            suppressor_release_cands = perm_off_stable if not target_in_perm else []
            print(f"  Candidate pools: suff={len(sufficiency_candidates)}  nec-act={len(necessity_act_candidates)}"
                  f"  nec-sup={len(necessity_sup_candidates)}")

            PER_HOP_BUDGET = max(0, (MAX_TOTAL_SECONDS - (time.perf_counter() - t_benchmark_start)) * 0.5)
            t_gene_tests = time.perf_counter()
            budget_exceeded = False

            sufficient_activators = []
            print(f"  Sufficiency test: {len(sufficiency_candidates)} candidates (budget: {PER_HOP_BUDGET:.0f}s)...")
            t0 = time.perf_counter()
            for candidate in sufficiency_candidates:
                if time.perf_counter() - t_gene_tests > PER_HOP_BUDGET:
                    print(f"    Stopped: per-hop budget reached"); budget_exceeded = True; break
                states, _, _ = simulate(bn_dict_pruned, all_zeros_sub, candidate, 1)
                if target_on_any(states) and not target_in_dark:
                    sufficient_activators.append(candidate)
            print(f"  [time] Sufficiency: {time.perf_counter()-t0:.3f}s")

            necessary_activators, redundant_activators, necessary_suppressors, suppressor_releases = [], [], [], []
            _remaining = max(0, MAX_TOTAL_SECONDS - (time.perf_counter() - t_benchmark_start))
            print(f"  BoNesis necessity test: {len(necessity_act_candidates)} act + {len(necessity_sup_candidates)} sup"
                  f" candidates ({n_pruned} nodes — budget: {_remaining:.0f}s)...")
            t_bn = time.perf_counter()
            _box = []
            def _run_bonesis_necessity():
                _base = target_on_any(
                    list(bonesis.BooleanNetwork(bn_dict_pruned).attractors(
                        reachable_from={g: 1 for g in bn_dict_pruned})))
                _nec, _red, _nsup, _sup = [], [], [], []
                # Activator necessity: KO from permissive
                for cand in necessity_act_candidates:
                    ko = dict(bn_dict_pruned); ko[cand] = "0"
                    ko_t = target_on_any(list(bonesis.BooleanNetwork(ko).attractors(
                        reachable_from={g: 1 for g in ko})))
                    cb = _base
                    if cand not in bn_dict_pruned:
                        pd2 = dict(bn_dict_pruned); pd2[cand] = "1"
                        cb = target_on_any(list(bonesis.BooleanNetwork(pd2).attractors(
                            reachable_from={g: 1 for g in pd2})))
                    if cb and not ko_t: _nec.append(cand)
                    elif cb and ko_t:   _red.append(cand)
                # Suppressor necessity: force ON from permissive
                for cand in necessity_sup_candidates:
                    forced = dict(bn_dict_pruned); forced[cand] = "1"
                    f_t = target_on_any(list(bonesis.BooleanNetwork(forced).attractors(
                        reachable_from={g: 1 for g in forced})))
                    if _base and not f_t: _nsup.append(cand)
                # Suppressor release: KO when target is OFF
                for cand in suppressor_release_cands:
                    ko = dict(bn_dict_pruned); ko[cand] = "0"
                    ko_t = target_on_any(list(bonesis.BooleanNetwork(ko).attractors(
                        reachable_from={g: 1 for g in ko})))
                    if not _base and ko_t: _sup.append(cand)
                _box.append((_base, _nec, _red, _nsup, _sup))
            _t = threading.Thread(daemon=True, target=_run_bonesis_necessity)
            _t.start(); _t.join(_remaining)
            if _box:
                baseline_perm_on, necessary_activators, redundant_activators, necessary_suppressors, suppressor_releases = _box[0]
                print(f"  [time] BoNesis necessity: {time.perf_counter()-t_bn:.3f}s")
            else:
                print(f"  BoNesis necessity timed out — budget exhausted after {time.perf_counter()-t_bn:.0f}s")
                raise _BudgetExhausted()

            sink_rules = {g: bn_dict[g] for g in sink_nodes}

            print(f"\n{'='*70}")
            print(f"  UPSTREAM REGULATORY ANALYSIS of '{target_gene}'")
            print(f"  {hops} hop(s)  |  BoNesis + synchronous simulation ({n_pruned} nodes)")
            print(f"  Update: synchronous  |  Rules: activators OR'd, suppressors AND NOT'd")
            print(f"  NOTE: simulation finds ONE attractor per starting condition.")
            print(f"\n  Step 1 — Structural direct regulators (1-hop):")
            print(f"    Activators: {len(direct_activators)}")
            for g in direct_activators:  print(f"      + {g:40s}  {tags(g)}")
            print(f"    Suppressors: {len(direct_suppressors)}")
            for g in direct_suppressors: print(f"      - {g:40s}  {tags(g)}")
            print(f"\n{'='*70}")
            print(f"  Step 2 — Attractor state comparison (dark vs permissive) [CORRELATION ONLY]")
            print(f"    Dark: target={'ON' if target_in_dark else 'OFF'}  |  Perm: target={'ON' if target_in_perm else 'OFF'}")
            print(f"  {'Hop':>4}  {'Total':>6}  {'Perm-ON/Dark-OFF':>17}  {'Perm-OFF/Dark-ON':>17}")
            for h in all_hops:
                genes_h = [g for g in all_upstream if hop_of.get(g) == h]
                print(f"  {h:>4}  {len(genes_h):>6}  {len(hop_pod[h]):>17}  {len(hop_opd[h]):>17}")
                for g in sorted(hop_pod[h])[:8]:
                    print(f"      + {g:40s}  [stably ON perm / OFF dark]  {tags(g)}")
                if len(hop_pod[h]) > 8: print(f"      ... and {len(hop_pod[h])-8} more at hop {h}")
                for g in sorted(hop_opd[h])[:5]:
                    print(f"      - {g:40s}  [stably OFF perm / ON dark]  {tags(g)}")
                if len(hop_opd[h]) > 5: print(f"      ... and {len(hop_opd[h])-5} more at hop {h}")
            print(f"\n{'='*70}")
            print(f"  Step 3 — Sufficient upstream activators:")
            if budget_exceeded: print(f"    NOTE: per-hop budget reached — results may be partial")
            print(f"  {len(sufficient_activators)} found")
            for g in sufficient_activators: print(f"    + {g:40s}  {tags(g)}")
            print(f"\n{'='*70}")
            print(f"  Step 4 — Necessary activators (KO from permissive, BoNesis):")
            if budget_exceeded: print(f"    NOTE: per-hop budget reached — results may be partial")
            print(f"\n  Necessary activators (KO turns '{target_gene}' OFF): {len(necessary_activators)}")
            for g in necessary_activators: print(f"    ! {g:40s}  [required activator]  {tags(g)}")
            print(f"\n  Redundant activators ('{target_gene}' stays ON after KO): {len(redundant_activators)}")
            for g in redundant_activators: print(f"    ~ {g:40s}  [dispensable]  {tags(g)}")
            print(f"\n{'='*70}")
            print(f"  Step 5 — Suppressor necessity and release (BoNesis):")
            print(f"\n  Necessary suppressors (forced ON turns '{target_gene}' OFF): {len(necessary_suppressors)}")
            for g in necessary_suppressors: print(f"    ! {g:40s}  [necessary suppressor]  {tags(g)}")
            print(f"\n  Suppressor release (KO turns '{target_gene}' ON): {len(suppressor_releases)}")
            for g in suppressor_releases: print(f"    ~ {g:40s}  [suppressor release]  {tags(g)}")
            if community:
                print(f"\n{'='*70}")
                cg = {}
                for g, cid in community.items(): cg.setdefault(cid, []).append(g)
                print(f"  REGULATORY MODULES (SBM) — {len(cg)} modules")
                for cid, members in sorted(cg.items(), key=lambda x: len(x[1]), reverse=True)[:5]:
                    print(f"\n  Module {cid}: {len(members)} nodes")
                    for g in sorted(members)[:8]: print(f"    {g:40s}  {tags(g)}")
                    if len(members) > 8: print(f"    ... and {len(members)-8} more")
            print(f"\n{'='*70}")
            print(f"  Sink nodes — upstream genes with no further upstream regulators in subgraph")
            print(f"  ({len(sink_nodes)} nodes). Candidate master regulators (no in-edges in model).")
            if sink_nodes:
                for g in sorted(sink_nodes):
                    val = eval_rule_simple(sink_rules[g], att_perm[0])
                    print(f"    {g:40s}  permissive state: {'ON' if val else 'OFF':3s}  {tags(g)}")
            print(f"\n======================================================================")
            print(f"  TIMING SUMMARY — {target_gene}, hop {hops}")
            print(f"    Total hop time: {time.perf_counter()-t_start:.3f}s  |  step times logged above")
            print(f"======================================================================")

        elapsed = time.perf_counter() - t_start
        mode_used = "BoNesis"
        timings.append((hops, elapsed, mode_used))
        print(f"{elapsed:.1f}s [{mode_used}] | {len(subgraph_nodes)} nodes → {os.path.basename(out_file)}")

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
