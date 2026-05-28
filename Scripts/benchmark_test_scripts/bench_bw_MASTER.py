#!/opt/anaconda3/envs/bachelor_env/bin/python
import csv, os, re, time, keyword, contextlib
import graph_tool.all as gt

GENE             = "MYB46"
MAX_TOTAL_SECONDS = 1800   # 0.5-hour budget; a new hop only starts if time remains
MAX_HOPS         =  30     # safety cap — time limit is the primary stop condition
MAX_SIGN_ITER    =  20     # max iterations for sign-propagation convergence
SCRIPT_NAME  = os.path.basename(__file__).replace('.py', '')
NETWORK_FILE = "filtered_GT_normalized.csv"
OUT_DIR      = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
BASE_DIR     = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.makedirs(OUT_DIR, exist_ok=True)

CATEGORIES_TO_KEEP = {"gene", "protein", "mutant", "metabolite", "process", "phenotype"}
BOOLEAN_RESERVED   = {"TRUE", "FALSE", "NOT", "AND", "OR"}
ACT_REL = "Activation / Induction / Causation / Result"
SUP_REL = "Repression / Inhibition / Negative Regulation"
CAT_MAP = {
    "gene": "Gene", "protein": "Protein", "mutant": "Mutant",
    "metabolite": "Metabolite", "process": "Process", "phenotype": "Phenotype",
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

# ── ONE-TIME DATA LOAD ─────────────────────────────────────────────────────────
print(f"Loading: {SCRIPT_NAME}")
activators, suppressors, out_act, out_sup = {}, {}, {}, {}
all_nodes, gene_category = set(), {}
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
        if rel == ACT_REL:
            activators.setdefault(t, set()).add(s); out_act.setdefault(s, set()).add(t); _seen_act.add((s, t))
        else:
            suppressors.setdefault(t, set()).add(s); out_sup.setdefault(s, set()).add(t); _seen_sup.add((s, t))
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

if GENE not in node_idx:
    print(f"ERROR: '{GENE}' not found in network. Exiting."); exit(1)
print(f"Gene '{GENE}' confirmed in network [{cat(GENE)}].")
print(f"Benchmark (BACKWARD): time limit {MAX_TOTAL_SECONDS}s | Method: signed path propagation\n")

# Reversed graph for BFS — built once, reused every hop
g_rev  = gt.GraphView(g_full, reversed=True)
dist_all = gt.shortest_distance(g_rev, source=g_full.vertex(node_idx[GENE]), directed=True)

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
        print(f"Time limit ({MAX_TOTAL_SECONDS}s) reached before hop {hops}. Stopping."); break

    print(f"  Hop {hops}... (remaining: {remaining:.0f}s)", end=" ", flush=True)
    t_start = time.perf_counter()
    out_file = os.path.join(OUT_DIR, f"{SCRIPT_NAME}_{GENE}_hops{hops}.txt")

    try:
        with open(out_file, 'w') as fout, contextlib.redirect_stdout(fout):
            print(f"# Benchmark: {SCRIPT_NAME}.py")
            print(f"# Gene: {GENE}  |  Hops: {hops}  |  Method: signed path propagation")
            print(f"# Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}")

            # ── Subgraph extraction ────────────────────────────────────────────
            hop_of = {node_list[i]: int(dist_all[i])
                      for i in range(g_full.num_vertices()) if int(dist_all[i]) <= hops}
            subgraph_nodes = set(hop_of.keys())
            shell_at = {h: sorted(n for n, d in hop_of.items() if d == h and n != GENE)
                        for h in range(1, hops + 1)}
            total_upstream = sum(len(shell_at[h]) for h in range(1, hops + 1))

            n_sub_edges = sum(1 for s, t in edges_act + edges_sup
                              if s in subgraph_nodes and t in subgraph_nodes)
            print(f"\n{'='*70}")
            print(f"  {GENE} — {hops} hop(s)  |  signed path propagation")
            print(f"  Upstream nodes: {total_upstream}  |  Subgraph edges: {n_sub_edges}")

            # ── Iterative sign propagation ─────────────────────────────────────
            net_sign    = {GENE: frozenset({+1})}
            shell_conns = {}
            t_prop = time.perf_counter()
            for _iter in range(MAX_SIGN_ITER):
                _changed = False
                for h in range(1, hops + 1):
                    for node in shell_at[h]:
                        signs = set()
                        conns = []
                        for dn in out_act.get(node, set()):
                            if dn in subgraph_nodes and net_sign.get(dn):
                                signs.update(net_sign[dn])
                                conns.append((dn, '+'))
                        for dn in out_sup.get(node, set()):
                            if dn in subgraph_nodes and net_sign.get(dn):
                                signs.update({-s for s in net_sign[dn]})
                                conns.append((dn, '-'))
                        new_sign = frozenset(signs)
                        if new_sign != net_sign.get(node, frozenset()):
                            _changed = True
                        net_sign[node]    = new_sign
                        shell_conns[node] = sorted(conns, key=lambda x: (hop_of.get(x[0], 999), x[0]))
                if not _changed:
                    print(f"  Sign propagation: converged in {_iter + 1} iteration(s)  "
                          f"({time.perf_counter()-t_prop:.3f}s)")
                    break
            else:
                print(f"  WARNING: sign propagation did not converge after {MAX_SIGN_ITER} iterations")

            def classify(node):
                s = net_sign.get(node, frozenset())
                if   s == frozenset({+1}): return "activator"
                elif s == frozenset({-1}): return "suppressor"
                elif {+1, -1} <= s:        return "mixed"
                else:                      return "no_path"

            def tags(g):
                return "[" + ", ".join([cat(g)]) + "]" + mol_flag(g)

            # ── Direct regulators ─────────────────────────────────────────────
            direct_act = sorted(g for g in activators.get(GENE, set()) if g in hop_of)
            direct_sup = sorted(g for g in suppressors.get(GENE, set()) if g in hop_of)
            print(f"\n  Direct regulators — activators: {len(direct_act)}  suppressors: {len(direct_sup)}")

            # ── Per-hop classification ─────────────────────────────────────────
            summary_act, summary_sup, summary_mix, summary_none = {}, {}, {}, {}
            for h in range(1, hops + 1):
                acts   = sorted(n for n in shell_at[h] if classify(n) == "activator")
                sups   = sorted(n for n in shell_at[h] if classify(n) == "suppressor")
                mixed  = sorted(n for n in shell_at[h] if classify(n) == "mixed")
                nopath = sorted(n for n in shell_at[h] if classify(n) == "no_path")
                for g in acts:   summary_act[g] = h
                for g in sups:   summary_sup[g] = h
                for g in mixed:  summary_mix[g] = h
                for g in nopath: summary_none[g] = h

                print(f"\n{'='*70}")
                print(f"  HOP {h}  |  {len(shell_at[h])} nodes  —  "
                      f"activators: {len(acts)}  suppressors: {len(sups)}  "
                      f"mixed: {len(mixed)}  no_path: {len(nopath)}")
                for g in acts:   print(f"    + {g:50s}  {tags(g)}")
                for g in sups:   print(f"    - {g:50s}  {tags(g)}")
                for g in mixed:  print(f"    ~ {g:50s}  {tags(g)}")
                if nopath:
                    print(f"    ? {len(nopath)} node(s) with no resolved signed path")

            # ── Summary ───────────────────────────────────────────────────────
            print(f"\n{'='*70}")
            print(f"  SUMMARY — {GENE}  {hops} hop(s)  |  signed path propagation")
            print(f"  Net activators  ({len(summary_act)}):  "
                  + ", ".join(f"{g} (h{h})" for g, h in sorted(summary_act.items(), key=lambda x: x[1])))
            print(f"  Net suppressors ({len(summary_sup)}):  "
                  + ", ".join(f"{g} (h{h})" for g, h in sorted(summary_sup.items(), key=lambda x: x[1])))
            if summary_mix:
                print(f"  Mixed effect    ({len(summary_mix)}):  "
                      + ", ".join(f"{g} (h{h})" for g, h in sorted(summary_mix.items(), key=lambda x: x[1])))
            print(f"  Total hop time: {time.perf_counter()-t_start:.3f}s")
            print(f"{'='*70}")

        elapsed = time.perf_counter() - t_start
        timings.append((hops, elapsed, len(subgraph_nodes), len(summary_act), len(summary_sup), len(summary_mix)))
        print(f"{elapsed:.3f}s  |  {total_upstream} upstream nodes  →  {os.path.basename(out_file)}")

    except Exception as e:
        import traceback
        elapsed = time.perf_counter() - t_start
        print(f"ERROR after {elapsed:.1f}s: {e}")
        timings.append((hops, elapsed, 0, 0, 0, 0))
        with open(out_file, 'a') as fout: fout.write(f"\n--- ERROR ---\n{traceback.format_exc()}\n")
        break

print(f"\n{'='*70}")
print(f"BENCHMARK SUMMARY: {SCRIPT_NAME}")
print(f"Gene: {GENE}  |  Method: signed path propagation")
print(f"{'Hops':>6}  {'Time (s)':>10}  {'Nodes':>7}  {'Act':>5}  {'Sup':>5}  {'Mix':>5}  Output file")
print(f"{'-'*6}  {'-'*10}  {'-'*7}  {'-'*5}  {'-'*5}  {'-'*5}  {'-'*45}")
for h, t, n, a, s, m in timings:
    print(f"{h:>6}  {t:>10.3f}  {n:>7}  {a:>5}  {s:>5}  {m:>5}  {SCRIPT_NAME}_{GENE}_hops{h}.txt")
