#!/opt/anaconda3/envs/bachelor_env/bin/python
import csv, os, re, time, keyword, contextlib, threading, sys
import graph_tool.all as gt
import bonesis

GENE                    = "MYB46"
MAX_TOTAL_SECONDS       = 1800   # 0.5-hour budget; a new hop only starts if time remains
MAX_HOPS                =  30     # safety cap — time limit is the primary stop condition
MAX_SIM_STEPS           = 1000
SCRIPT_NAME  = os.path.basename(__file__).replace('.py', '')
NETWORK_FILE = "filtered_GT_normalized.csv"
OUT_DIR      = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
BASE_DIR    = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.makedirs(OUT_DIR, exist_ok=True)

CATEGORIES_TO_KEEP = {
    "gene",
    "protein",
    "mutant",
    "metabolite",
    "process",
    "phenotype",
}
BOOLEAN_RESERVED = {"TRUE", "FALSE", "NOT", "AND", "OR"}
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

def clean_name(name):
    if not name or not isinstance(name, str): return "unknown"
    name = re.sub(r"_+", "_", re.sub(r"[^a-zA-Z0-9_]", "_", name)).strip("_")
    if not name: return "unknown"
    if name[0].isdigit(): return "n" + name
    if name.upper() in BOOLEAN_RESERVED or keyword.iskeyword(name): return "node_" + name
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
            for s, t in set((s,t) for s,t in edges_act if s in sub_v_idx and t in sub_v_idx):
                g_sub.add_edge(sub_v_idx[s], sub_v_idx[t])
            for s, t in set((s,t) for s,t in edges_sup if s in sub_v_idx and t in sub_v_idx):
                g_sub.add_edge(sub_v_idx[s], sub_v_idx[t])
            sub_comp, sub_hist = gt.label_components(g_sub)
            large_sub = int((sub_hist > 1).sum())
            print(f"  SCCs: {sub_hist.shape[0]}  |  non-trivial: {large_sub}")
            if large_sub:
                scc_m = [sub_node_list[i] for i in range(g_sub.num_vertices()) if int(sub_comp[i]) == int(sub_hist.argmax())]
                print(f"  Largest SCC ({int(sub_hist.max())} nodes):")
                for m in sorted(scc_m)[:10]: print(f"    {m:40s}  [{cat(m)}]")
                if len(scc_m) > 10: print(f"    ... and {len(scc_m)-10} more")

            # Boolean rules: activators OR'd (any activator sufficient); suppressors AND NOT'd (any suppressor dominant).
            # Rule form: (act1 | act2 | ...) & !sup1 & !sup2 & ...  — models functional redundancy + dominant repression.
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
                return "[" + ", ".join(parts) + "]" + mol_flag(g)

            all_zeros_sub = {g: 0 for g in subgraph_nodes}
            all_ones_sub  = {g: 1 for g in subgraph_nodes}

            def target_on_any(atts):  return any(a.get(target_gene, 0) == 1 for a in atts)
            def stable_on(g, atts):  return bool(atts) and all(a.get(g, 0) == 1 for a in atts)
            def stable_off(g, atts): return bool(atts) and all(a.get(g, 0) == 0 for a in atts)

            direct_activators  = sorted(g for g in activators.get(target_gene, set()) if g in subgraph_nodes)
            direct_suppressors = sorted(g for g in suppressors.get(target_gene, set()) if g in subgraph_nodes)
            all_upstream       = sorted(subgraph_nodes - {target_gene})

            # Two baseline simulations — target gene locked, mirroring the forward analysis.
            # Dark: lock target=0 → "what upstream state when target is OFF?"
            # Perm: lock target=1 → "what upstream state when target is ON?"
            t0 = time.perf_counter()
            att_dark, _, dark_conv = simulate(bn_dict_pruned, all_zeros_sub, locked=target_gene, val=0)
            att_perm, _, perm_conv = simulate(bn_dict_pruned, all_ones_sub,  locked=target_gene, val=1)
            print(f"  [time] Baseline simulations: {time.perf_counter()-t0:.3f}s")
            if not dark_conv: print(f"  WARNING: dark attractor did not converge within {MAX_SIM_STEPS} steps")
            if not perm_conv: print(f"  WARNING: perm attractor did not converge within {MAX_SIM_STEPS} steps")
            print(f"  NOTE: simulation finds ONE attractor per starting condition.")
            print(f"  NOTE: target gene locked to 0 (dark) and 1 (perm) — same method as forward analysis.")
            target_in_dark = target_on_any(att_dark)   # always False (locked=0)
            target_in_perm = target_on_any(att_perm)   # always True  (locked=1)
            _lbl = lambda atts, c: ("fixed point" if len(atts)==1 else f"cycle/{len(atts)}") + (" (conv)" if c else " (max)")
            print(f"  Dark (target locked OFF): {_lbl(att_dark, dark_conv)}")
            print(f"  Perm (target locked ON):  {_lbl(att_perm, perm_conv)}")

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

            # Cap per-gene budget at 120 s so one large hop cannot consume the whole run.
            PER_HOP_BUDGET = max(0, (MAX_TOTAL_SECONDS - (time.perf_counter() - t_benchmark_start)) * 0.5)
            t_gene_tests = time.perf_counter()
            budget_exceeded = False

            # _con writes to the real terminal even though stdout is redirected to the file.
            def _con(msg):
                sys.__stdout__.write(msg); sys.__stdout__.flush()

            _con(f"  pruned BN: {n_pruned} nodes  |  target {'ON' if target_in_perm else 'OFF'} in perm attractor\n")
            _con(f"  candidates: suff={len(sufficiency_candidates)}  nec-act={len(necessity_act_candidates)}"
                 f"  nec-sup={len(necessity_sup_candidates)}  budget={PER_HOP_BUDGET:.0f}s\n")

            sufficient_activators = []
            print(f"  Sufficiency test: {len(sufficiency_candidates)} candidates (budget: {PER_HOP_BUDGET:.0f}s)...")
            _con(f"    [sufficiency  {len(sufficiency_candidates)} cands] ")
            t0 = time.perf_counter()
            for candidate in sufficiency_candidates:
                if time.perf_counter() - t_gene_tests > PER_HOP_BUDGET:
                    print(f"    Stopped: per-hop budget reached"); budget_exceeded = True; break
                states, _, _ = simulate(bn_dict_pruned, all_zeros_sub, candidate, 1)
                if target_on_any(states) and not target_in_dark:
                    sufficient_activators.append(candidate)
            _con(f"done {len(sufficient_activators)} found  {time.perf_counter()-t0:.1f}s\n")
            print(f"  [time] Sufficiency: {time.perf_counter()-t0:.3f}s")

            necessary_activators, redundant_activators, necessary_suppressors, suppressor_releases = [], [], [], []
            _remaining = min(PER_HOP_BUDGET, max(0, MAX_TOTAL_SECONDS - (time.perf_counter() - t_benchmark_start)))
            print(f"  BoNesis necessity test: {len(necessity_act_candidates)} act + {len(necessity_sup_candidates)} sup"
                  f" candidates ({n_pruned} nodes — budget: {_remaining:.0f}s)...")
            _con(f"    [BoNesis nec  act={len(necessity_act_candidates)} sup={len(necessity_sup_candidates)}  budget={_remaining:.0f}s] ")
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
                _con(f"done  {time.perf_counter()-t_bn:.1f}s\n")
                print(f"  [time] BoNesis necessity: {time.perf_counter()-t_bn:.3f}s")
            else:
                _con(f"timed out after {time.perf_counter()-t_bn:.0f}s — stopping hop\n")
                print(f"  BoNesis necessity timed out after {time.perf_counter()-t_bn:.0f}s — stopping (BoNesis-only benchmark)")
                raise _BudgetExhausted()
            nec_mode_str = "BoNesis"

            sink_rules = {g: bn_dict[g] for g in sink_nodes}

            print(f"\n{'='*70}")
            print(f"  UPSTREAM REGULATORY ANALYSIS of '{target_gene}'")
            print(f"  {hops} hop(s)  |  {nec_mode_str} ({n_pruned} nodes)")
            print(f"  Update: synchronous  |  Rules: activators OR'd, suppressors AND NOT'd")
            print(f"  ! = non-molecular node (mutant / process / phenotype)")
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
            print(f"  Step 3 — Sufficient upstream activators (single-gene activation drives '{target_gene}' ON):")
            if budget_exceeded: print(f"    NOTE: per-hop budget reached — results may be partial")
            if sufficient_activators:
                for h in all_hops:
                    sa_h = sorted(g for g in sufficient_activators if hop_of.get(g) == h)
                    if sa_h:
                        print(f"    Hop {h}:")
                        for g in sa_h: print(f"      ++ {g:40s}  {tags(g)}")
            else:
                print(f"    None found")
            print(f"\n{'='*70}")
            print(f"  Step 4 — Activator necessity ({nec_mode_str} — KO from permissive):")
            if budget_exceeded: print(f"    NOTE: per-hop budget reached — results may be partial")
            print(f"\n  Necessary activators (KO turns '{target_gene}' OFF): {len(necessary_activators)}")
            if necessary_activators:
                for h in all_hops:
                    na_h = sorted(g for g in necessary_activators if hop_of.get(g) == h)
                    if na_h:
                        print(f"    Hop {h}:")
                        for g in na_h: print(f"      +  {g:40s}  [required]  {tags(g)}")
            print(f"\n  Redundant activators ('{target_gene}' stays ON after KO): {len(redundant_activators)}")
            _red_shown = 0
            for h in all_hops:
                if _red_shown >= 20: break
                red_h = sorted(g for g in redundant_activators if hop_of.get(g) == h)
                if red_h:
                    print(f"    Hop {h}:")
                    for g in red_h:
                        if _red_shown >= 20: break
                        print(f"      ~  {g:40s}  [dispensable]  {tags(g)}"); _red_shown += 1
            if len(redundant_activators) > 20: print(f"    ... and {len(redundant_activators)-20} more [dispensable]")
            print(f"\n{'='*70}")
            print(f"  Step 5 — Suppressor necessity and release ({nec_mode_str}):")
            print(f"\n  Necessary suppressors (forced ON turns '{target_gene}' OFF): {len(necessary_suppressors)}")
            if necessary_suppressors:
                for h in all_hops:
                    ns_h = sorted(g for g in necessary_suppressors if hop_of.get(g) == h)
                    if ns_h:
                        print(f"    Hop {h}:")
                        for g in ns_h: print(f"      -! {g:40s}  [required suppressor]  {tags(g)}")
            print(f"\n  Suppressor release (KO turns '{target_gene}' ON): {len(suppressor_releases)}")
            if suppressor_releases:
                for h in all_hops:
                    sr_h = sorted(g for g in suppressor_releases if hop_of.get(g) == h)
                    if sr_h:
                        print(f"    Hop {h}:")
                        for g in sr_h: print(f"      -~ {g:40s}  [release]  {tags(g)}")
            print(f"\n{'='*70}")
            print(f"  Step 6 — REGULATORY SUMMARY BY HOP ({nec_mode_str})")
            print(f"  Classifies every upstream node by predicted regulatory role.")
            print(f"  ++ = sufficient activator  +  = necessary activator   ~  = redundant")
            print(f"  -! = necessary suppressor  -~ = suppressor release    .+ = correlated activator (attractor only)")
            print(f"  .- = correlated suppressor (attractor only)           no mark = no role detected")
            print(f"  Use this section to compare against known MYB46 regulators hop by hop.")
            _suff_set    = set(sufficient_activators)
            _nec_act_set = set(necessary_activators)
            _red_act_set = set(redundant_activators)
            _nec_sup_set = set(necessary_suppressors)
            _sup_rel_set = set(suppressor_releases)
            _pod_set     = {g for nodes in hop_pod.values() for g in nodes}
            _opd_set     = {g for nodes in hop_opd.values() for g in nodes}
            for h in all_hops:
                nodes_at_h = sorted(g for g in all_upstream if hop_of.get(g) == h)
                _func, _corr, _neut = [], [], []
                for g in nodes_at_h:
                    if   g in _suff_set:    _func.append((g, "++", "sufficient activator"))
                    elif g in _nec_act_set: _func.append((g, "+ ", "necessary activator"))
                    elif g in _red_act_set: _func.append((g, "~ ", "redundant activator"))
                    elif g in _nec_sup_set: _func.append((g, "-!", "necessary suppressor"))
                    elif g in _sup_rel_set: _func.append((g, "-~", "suppressor release"))
                    elif g in _pod_set:     _corr.append((g, ".+", "correlated activator"))
                    elif g in _opd_set:     _corr.append((g, ".-", "correlated suppressor"))
                    else:                   _neut.append(g)
                print(f"\n  Hop {h}  ({len(nodes_at_h)} nodes — {len(_func)} functional  {len(_corr)} correlated  {len(_neut)} neutral):")
                for g, sym, role in sorted(_func, key=lambda x: x[0]):
                    print(f"    {sym}  {g:40s}  [{role}]  {tags(g)}")
                for g, sym, role in sorted(_corr, key=lambda x: x[0]):
                    print(f"    {sym}  {g:40s}  [{role}]  {tags(g)}")
                if _neut:
                    print(f"    ..  {len(_neut)} neutral node(s) at hop {h}")
            print(f"\n{'='*70}")
            print(f"  Sink nodes — upstream genes with no further upstream regulators in subgraph")
            print(f"  ({len(sink_nodes)} nodes). Candidate master regulators (no in-edges in model).")
            if sink_nodes:
                for g in sorted(sink_nodes)[:30]:
                    val = eval_rule_simple(sink_rules[g], att_perm[0])
                    print(f"    {g:40s}  permissive state: {'ON' if val else 'OFF':3s}  {tags(g)}")
                if len(sink_nodes) > 30: print(f"    ... and {len(sink_nodes)-30} more")
            print(f"\n======================================================================")
            print(f"  TIMING SUMMARY — {target_gene}, hop {hops}")
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
