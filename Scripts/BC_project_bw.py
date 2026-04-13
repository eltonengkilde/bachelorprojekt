# Cell A
# Import of standard CSV module for reading CSV files in python
# import os and BASE_DIR at the top — this always points to the project folder regardless of where VS Code runs Python from.
# Both file paths (kg_all_2005-2025_1000papers_combined.csv and network.txt) now use BASE_DIR so they always resolve correctly.
import csv
import os
import re
import mpbn

# Base directory — always the folder where this script lives, regardless of where VS Code launches from
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Load pre-filtered data produced by the separate filtering script
with open(os.path.join(BASE_DIR, "networks_used_by_scripts", "filtered_large.csv"), newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    last_rows = list(reader)

print(f"Loaded {len(last_rows)} rows from filtered network")

# boolean.py reserved words that cannot be used as gene names
BOOLEAN_RESERVED = {"TRUE", "FALSE", "NOT", "AND", "OR"}
# Helper: clean a name so mpbn can parse it
def clean_name(name):
    if not name or not isinstance(name, str):  # guard against None, empty string, or non-string input
        return "unknown"

    name = re.sub(r"[^a-zA-Z0-9_]", "_", name)  # replace any illegal character with underscore
    name = re.sub(r"_+", "_", name)              # collapse multiple consecutive underscores into one
    name = name.strip("_")                        # remove leading/trailing underscores

    if not name:  # catch the edge case where cleaning stripped the entire name to an empty string
        return "unknown"
    elif name[0].isdigit():   # mpbn variable names cannot start with a digit — prefix with "g"
        return "g" + name
    elif name.upper() in BOOLEAN_RESERVED:  # avoid clashing with mpbn parser keywords
        return "gene_" + name
    return name

# Cell H
# Group rows by target, collecting activators and suppressors
activators = {}
suppressors = {}

# Maps the full relationship category strings from the KG to short internal labels.
# "act" = activation, "sup" = suppression. Any relationship not in this map is skipped.
RELATIONSHIP_MAPPING = {
    "Activation / Induction / Causation / Result": "act",
    "Repression / Inhibition / Negative Regulation": "sup",
}

# Loop through every row and sort each source into either the activators or suppressors
# dictionary for that target. Each target ends up with a set of activating sources
# and a set of suppressing sources, which are then used to build its Boolean formula.
for row in last_rows:
    rel = RELATIONSHIP_MAPPING.get(row["relationship_category"])
    if rel is None:
        continue  # skip anything that isn't act or sup (is this nessary given the previous filtering? maybe not but it adds safety)

    target = clean_name(row["target"])
    source = clean_name(row["source"])

    if rel == "act":
        activators.setdefault(target, set()).add(source)
    else:
        suppressors.setdefault(target, set()).add(source)

# Build Boolean expression for each target by looping through each target and gluing or logic to all activators and and logic to all supressors of the target and also adda a ! (not) to all these supressors as we want them to not be active for target active to be true.
# Logic: (A1 | A2) & !S1 & !S2 (this means that either A1 or A2 will activate a target as long as any supressor of the target is not active thus overriding the signal with a supression)
bn_dict = {}

# Union of all nodes that appear as a target in any activation or suppression relationship.
# Every node in this set needs a Boolean formula built for it.
all_targets = set(activators) | set(suppressors)

for target in all_targets:
    act = " | ".join(activators.get(target, []))
    sup = " & ".join("!" + s for s in suppressors.get(target, []))

# This is where the boolean logic can later be made more complex.
# This loop sorts the logic. if theres both activators and supressors for a gene write it like the first part, else if it only has activators ignore the supressors and write only the or logic for the activators, else only write the supressor logic. this logic is given by the code above.
    if act and sup:
        bn_dict[target] = "(" + act + ") & " + sup
    elif act:
        bn_dict[target] = act # if only activators, then the rule is just the or logic of the activators
    else:
        bn_dict[target] = sup

# Source-only nodes are self-sustaining.
# The code above only creates rules for the targets on the background of the act/sup information from the relationship column that is modified. but mpbn needs rules for all nodes to work. this chunk creates rules for the sources that does not appear as targets and thus doesnt have a relationship rule
# these will therefore not appear in the bn_dict, it then adds these sources to the bn_dict as well but with the rule "Genex":"Genex" meaning that gene x activates itself.
# these genes are the control genes that can be set for 0 or 1 to check what the consequence is downstream.
for row in last_rows:
    cleaned_source_name = clean_name(row["source"])
    if cleaned_source_name not in bn_dict:
        bn_dict[cleaned_source_name] = cleaned_source_name #!!!is this self activation a problem resulting in the large output, this should be investigated.

# Remove source-only nodes (rule = self) — no incoming regulation, no network context
# Save them first: they are the external inputs and matter for the backwards analysis.
source_only_nodes = {g for g, f in bn_dict.items() if f == g}
source_only = len(source_only_nodes)
bn_dict = {g: f for g, f in bn_dict.items() if f != g}

print(f"Source-only nodes removed: {source_only}")
print(f"Regulated nodes remaining: {len(bn_dict)}")
print("\nSample (first 10):")
for i, (gene, formula) in enumerate(list(bn_dict.items())[:10]):
    print(f"  {gene} = {formula}")
# Can be different each time because of the way python rendomises the internal order of sets every time a program is run. No relationships are changed though. it is a safety feature of python that it is scrambled which is called hash randomization and prevents attacks. The set wont change but can be scrambled.
# set is used to avoid duplicates but the order of the elements in the set is not guaranteed. This prevents broken boolean rules so that the same gene does not appear many times and confuse the program.

# Count total activation and suppression edges across all targets
total_act_edges = sum(len(srcs) for srcs in activators.values())
total_sup_edges = sum(len(srcs) for srcs in suppressors.values())
total_edges = total_act_edges + total_sup_edges
print(f"\nActivation edges:  {total_act_edges} ({100*total_act_edges//total_edges}%)")
print(f"Suppression edges: {total_sup_edges} ({100*total_sup_edges//total_edges}%)")

# Print a few example edges as a sanity check so the user can verify
# that source → target relationships were parsed in the correct direction.
print("\nExample activations (first 3 targets):")
for target, srcs in list(activators.items())[:3]:
    print(f"  {list(srcs)[:3]} --> {target}")

print("\nExample suppressions (first 3 targets):")
for target, srcs in list(suppressors.items())[:3]:
    print(f"  {list(srcs)[:3]} --| {target}")


# User Input step ______________________________________________________
# BACKWARDS MODE: ask which gene the user wants to activate, then find what upstream genes can do it.
# Keeps asking until the user enters a gene name that exists in the network
while True:
    target_gene = input("\nTarget gene to activate (what gene do you want to turn ON?): ")
    if target_gene not in bn_dict:
        matches = [g for g in bn_dict if target_gene.upper() in g.upper()]
        print(f"WARNING: '{target_gene}' not found in network! Try again.")
        if matches:
            print(f"  Did you mean one of these? {matches[:5]}")
    else:
        print(f"OK: '{target_gene}' found in network")
        break

# Ask how many regulatory steps upstream to search from the target gene.
# 1 hop = direct regulators only. More hops = larger subnetwork but slower analysis.
while True:
    MAX_HOPS = int(input("\nNumber of hops to explore upstream from target gene: "))
    if MAX_HOPS < 1:
        print(f"Must be a whole number larger than 0")
    else:
        print(f"OK: Searching {MAX_HOPS} hops upstream of '{target_gene}'")
        break

# Safety gate — requires the user to type "yolo" before the slow mpbn computation begins.
executor = input('type and enter "yolo" to run the mpbn attractor analysis')
if executor == "yolo":
    pass
else:
    print("Invalid input, exiting.")
    exit()

# Build upstream map: gene -> set of genes that appear in its Boolean formula (its regulators).
# Include source_only_nodes in the lookup — they are external inputs that appear in formulas
# but were stripped from bn_dict. Without them the BFS misses the first upstream layer entirely.
all_known_genes = set(bn_dict) | source_only_nodes
upstream_map = {}
for gene, formula in bn_dict.items():
    regulators = set(re.findall(r'\b[a-zA-Z_]\w*\b', str(formula))) & all_known_genes
    upstream_map[gene] = regulators

# BFS backwards from target_gene up to MAX_HOPS steps upstream
nodes_before_subnetwork = len(bn_dict)
upstream_subnetwork = {target_gene}
frontier = {target_gene}
for hop in range(MAX_HOPS):
    next_frontier = set()
    for node in frontier:
        for regulator in upstream_map.get(node, set()):
            if regulator not in upstream_subnetwork:
                next_frontier.add(regulator)
    upstream_subnetwork |= next_frontier
    frontier = next_frontier

bn_dict = {g: f for g, f in bn_dict.items() if g in upstream_subnetwork}
nodes_after_subnetwork = len(bn_dict)
nodes_removed_subnetwork = nodes_before_subnetwork - nodes_after_subnetwork
print(f"Nodes before upstream filter: {nodes_before_subnetwork}")
print(f"Nodes removed (not upstream of '{target_gene}'): {nodes_removed_subnetwork}")
print(f"Nodes remaining in upstream subnetwork ({MAX_HOPS} hops): {nodes_after_subnetwork}")

# Build the unperturbed Boolean network from the upstream subnetwork — used for attractor searches in Step 2.
bn = mpbn.MPBooleanNetwork(bn_dict)

# Step 1: direct regulators of target_gene from the raw KG data
direct_activators_of_target = sorted(activators.get(target_gene, set()))
direct_suppressors_of_target = sorted(suppressors.get(target_gene, set()))

# Print the analysis header and begin reporting results for each step.
print(f"\n{'='*80}")
print(f"  Backwards analysis: what activates '{target_gene}'?")

print(f"\n Step 1: Direct activators of '{target_gene}' (from KG): {len(direct_activators_of_target)} total")
for gene in direct_activators_of_target:
    print(f"    {gene}")

print(f"\n Step 1: Direct suppressors of '{target_gene}' (from KG): {len(direct_suppressors_of_target)} total")
for gene in direct_suppressors_of_target:
    print(f"    {gene}")

# Step 2: baseline — all genes OFF, is target_gene naturally ON?
# Then for each upstream gene, set it to 1 (all others 0) and check if target turns ON.
baseline_state = {gene: 0 for gene in bn_dict}
baseline_target_on = any(
    a.get(target_gene, 0) == 1
    for a in bn.attractors(reachable_from=baseline_state)
)

# All upstream nodes except the target itself — these are the candidates to test in Step 2.
upstream_genes = sorted(upstream_subnetwork - {target_gene})

single_activators = []    # setting this gene ON (alone) turns target ON

print(f"\n Step 2: Testing {len(upstream_genes)} upstream genes one at a time (all others OFF)...")
for candidate in upstream_genes:
    # Source-only nodes have no rule in bn_dict — add them back locked to "1" for this test.
    if candidate in source_only_nodes:
        test_bn_dict = dict(bn_dict)
        test_bn_dict[candidate] = "1"
        bn_test = mpbn.MPBooleanNetwork(test_bn_dict)
        test_state = {gene: 0 for gene in test_bn_dict}
        test_state[candidate] = 1
    else:
        bn_test = bn
        test_state = {gene: 0 for gene in bn_dict}
        test_state[candidate] = 1

    # Run attractor search from this test state and check if target_gene is ON in any attractor.
    # If the target is ON here but was OFF in the baseline, this candidate is a sufficient activator.
    candidate_target_on = any(
        a.get(target_gene, 0) == 1
        for a in bn_test.attractors(reachable_from=test_state)
    )
    if candidate_target_on and not baseline_target_on:
        single_activators.append(candidate)

# Print Step 2 results: genes that are individually sufficient to activate the target from silence.
print(f"\n  Genes that, when activated alone, turn ON '{target_gene}': {len(single_activators)} total")
for gene in single_activators:
    print(f"    {gene}")

# Step 3: Necessity test — knockout each direct activator and suppressor one at a time.
# Start from all-ones (everything active = permissive background), clamp the candidate to 0,
# and check whether the target goes OFF in the resulting attractors.
# This mirrors standard in silico knock-out experiments used in Boolean network analysis:
# "A prominent structural disruption, applied in different studies of BNs, is fixing a compound
#  of the system to either false or true. This kind of perturbation corresponds to in silico
#  knock-out or overexpression experiments."
# — Schwab et al. (2020), Computational and Structural Biotechnology Journal
#   https://pmc.ncbi.nlm.nih.gov/articles/PMC7096748/
#
# Necessity (knockout, clamp to 0) directly answers "what do I NEED?" because it reveals
# which activators the target cannot compensate for. Sufficiency (Step 2, clamp to 1 from
# all-zeros) answers a different question: "what alone is enough?" — which is overly strict
# since it ignores all co-factors present in real cellular contexts.

# Starting state for Step 3: all nodes ON (maximally permissive background).
# Knockout tests knock one candidate out from this state to test necessity.
all_on_start = {gene: 1 for gene in bn_dict}

# Check baseline: is target ON when everything is active (no knockouts)?
baseline_all_on_target_on = any(
    a.get(target_gene, 0) == 1
    for a in bn.attractors(reachable_from=all_on_start)
)

# Test every direct activator and suppressor for necessity
candidates_to_test = set(direct_activators_of_target) | set(direct_suppressors_of_target)
necessary_activators = []    # knocking this out turns target OFF
redundant_activators = []    # target stays ON even without this gene
suppressor_releases = []     # knocking this suppressor out turns target ON (de-repression)

print(f"\n Step 3: Necessity test: knocking out each direct regulator (all-ones background)...")
for candidate in sorted(candidates_to_test):
    # For source-only nodes (not in bn_dict), the global all-ones baseline doesn't include
    # them at all, so it can't serve as the "present" reference. Build a per-candidate
    # baseline that locks the node to "1" so the comparison is meaningful.
    if candidate in source_only_nodes:
        present_dict = dict(bn_dict)
        present_dict[candidate] = "1"
        bn_present = mpbn.MPBooleanNetwork(present_dict)
        present_start = {gene: 1 for gene in present_dict}
        candidate_baseline_on = any(
            a.get(target_gene, 0) == 1
            for a in bn_present.attractors(reachable_from=present_start)
        )
    else:
        candidate_baseline_on = baseline_all_on_target_on

    # Knockout: clamp candidate to 0.
    ko_dict = dict(bn_dict)
    ko_dict[candidate] = "0"
    bn_ko = mpbn.MPBooleanNetwork(ko_dict)
    ko_start = {gene: 1 for gene in ko_dict}
    ko_start[candidate] = 0

    ko_target_on = any(
        a.get(target_gene, 0) == 1
        for a in bn_ko.attractors(reachable_from=ko_start)
    )

    if candidate in direct_activators_of_target:
        if candidate_baseline_on and not ko_target_on:
            necessary_activators.append(candidate)   # target goes OFF without it
        elif candidate_baseline_on and ko_target_on:
            redundant_activators.append(candidate)   # target stays ON — another activator covers

    if candidate in direct_suppressors_of_target:
        if not candidate_baseline_on and ko_target_on:
            suppressor_releases.append(candidate)    # removing this suppressor releases the target

# Print Step 3 results — three categories based on what happened when each regulator was knocked out.

# Necessary: knocking this activator out turned the target OFF — the target depends on it.
print(f"\n  Necessary activators (knockout turns '{target_gene}' OFF): {len(necessary_activators)} total")
for gene in necessary_activators:
    print(f"    ! {gene}  [required]")

# A redundant activator is a gene that the KG says directly activates my target, but when it is knowked out in the boolean network the target stays ON anyway. this is because of the "target's" boolean rule because of the OR logic between activators.
print(f"\n  Redundant activators ('{target_gene}' stays ON without them): {len(redundant_activators)} total")
for gene in redundant_activators:
    print(f"    {gene}  [dispensable — another activator compensates]")

# Suppressor release: this suppressor was keeping the target OFF — removing it de-represses the target.
print(f"\n  Suppressors whose removal de-represses '{target_gene}': {len(suppressor_releases)} total")
for gene in suppressor_releases:
    print(f"    {gene}  [knocking this out turns target ON]")

print(f"\n{'='*80}")