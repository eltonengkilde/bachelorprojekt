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
with open(os.path.join(BASE_DIR, "networks_used_by_scripts", "filtered_largerexperiment_normalized.csv"), newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    last_rows = list(reader)

print(f"Loaded {len(last_rows)} rows from filtered network")

# Category overview — shows all unique source and target categories in the loaded data
# Use this to identify which categories to map to genes, phenotypes, and metabolites
source_counts = {}
target_counts = {}
for row in last_rows:
    source_counts[row["source_category"]] = source_counts.get(row["source_category"], 0) + 1
    target_counts[row["target_category"]] = target_counts.get(row["target_category"], 0) + 1

print(f"\nUnique source categories ({len(source_counts)} total):")
for category, count in sorted(source_counts.items(), key=lambda x: x[1], reverse=True):
    print(f"  {count:>6}  {category}")

print(f"\nUnique target categories ({len(target_counts)} total):")
for category, count in sorted(target_counts.items(), key=lambda x: x[1], reverse=True):
    print(f"  {count:>6}  {category}")

# DATA fitting
# boolean.py reserved words that cannot be used as gene names
BOOLEAN_RESERVED = {"TRUE", "FALSE", "NOT", "AND", "OR"}
# Helper: clean a name so mpbn can parse it
def clean_name(name):
    if not name or not isinstance(name, str):
        return "unknown"

    name = re.sub(r"[^a-zA-Z0-9_]", "_", name)  # replace any illegal character with underscore
    name = re.sub(r"_+", "_", name)              # collapse multiple consecutive underscores into one
    name = name.strip("_")                        # remove leading/trailing underscores

    if not name:  # catch the edge case where cleaning stripped the entire name to an empty string
        return "unknown"
    elif name[0].isdigit():
        return "g" + name
    elif name.upper() in BOOLEAN_RESERVED:
        return "gene_" + name
    return name

#Data preprocessing
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

# Keep only rows where both source and target belong to the relevant biological categories
CATEGORIES_TO_KEEP = {
    "Gene / Protein",
    "Phenotype / Trait / Disease",
    "Chemical / Metabolite / Cofactor / Ligand",
    "Biological Process / Pathway / Function / Regulatory / Signaling Mechanism",
}
rows_before = len(last_rows)
last_rows = [r for r in last_rows if r["source_category"] in CATEGORIES_TO_KEEP
                                  and r["target_category"] in CATEGORIES_TO_KEEP]
print(f"Rows after category filter: {len(last_rows)} (removed {rows_before - len(last_rows)}) \n Keeping categories: genes/proteins, Phenotypes, metabolites and pathways")

# Build a lookup: cleaned gene name -> category (source or target category from the KG)
gene_category = {}
for row in last_rows:
    gene_category[clean_name(row["source"])] = row["source_category"]
    gene_category[clean_name(row["target"])] = row["target_category"]


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

all_targets = set(activators) | set(suppressors)

for target in all_targets:
    act = " | ".join(activators.get(target, []))
    sup = " & ".join("!" + s for s in suppressors.get(target, []))

# This is where the boolean logic is made.
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
source_only = sum(1 for g, f in bn_dict.items() if f == g)
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

# Analysis
# User Input step ______________________________________________________
# SET THE GENE TO ACTIVATE AND THE TARGET TO CHECK
# Keeps asking until the user enters a gene name that exists in the network
while True:
    source_gene = input("\nSource gene to investigate: ") # Gene to activate # Use "VND_and_NST_proteins" or "AtMYB46"
    if source_gene not in bn_dict:
        matches = [g for g in bn_dict if source_gene.upper() in g.upper()]
        print(f"WARNING: '{source_gene}' not found in network! Try again.")
        if matches:
            print(f"  Did you mean one of these? {matches[:5]}")
    else:
        print(f"OK: '{source_gene}' found in network")
        break

# Ask how many regulatory steps downstream to explore from the source gene.
# 1 hop = direct targets only. More hops = larger subnetwork but slower analysis.
while True:
    MAX_HOPS = int(input("\nNumber of hops to explore downstream from source gene: "))
    if MAX_HOPS < 1:
        print(f"Must be a whole number larger than 0")
    else:
        print(f"OK: Continuing with subnetwork with the depth of '{MAX_HOPS}' downstream from source gene")
        break


# Cell J
# On the basis of mareks code this code aims to predict the immedeate downstream consequences of activating a specific gene from the KG
executor = input('type and enter "yolo" to run the mpbn attractor analysis')
if executor == "yolo":
    pass
else:
    print("Invalid input, exiting.")
    exit()

# Restrict bn_dict to the subnetwork reachable downstream from source_gene.
# First build a forward map: regulator -> set of genes whose formula mentions it.
downstream_map = {g: set() for g in bn_dict}
for gene, formula in bn_dict.items():
    for token in set(re.findall(r'\b[a-zA-Z_]\w*\b', str(formula))):
        if token in downstream_map:
            downstream_map[token].add(gene)

# BFS from source_gene up to MAX_HOPS steps downstream
nodes_before_subnetwork = len(bn_dict)
subnetwork = {source_gene}
frontier = {source_gene}
for hop in range(MAX_HOPS):
    next_frontier = set()
    for node in frontier:
        for target in downstream_map.get(node, set()):
            if target not in subnetwork:
                next_frontier.add(target)
    subnetwork |= next_frontier
    frontier = next_frontier

bn_dict = {g: f for g, f in bn_dict.items() if g in subnetwork}
nodes_after_subnetwork = len(bn_dict)
nodes_removed_subnetwork = nodes_before_subnetwork - nodes_after_subnetwork
print(f"Nodes before subnetwork filter: {nodes_before_subnetwork}")
print(f"Nodes removed (not reachable from '{source_gene}'): {nodes_removed_subnetwork}")
print(f"Nodes remaining in subnetwork ({MAX_HOPS} hops): {nodes_after_subnetwork}")

print(f"\n{'='*80}")
# Build two versions of the network with source_gene locked OFF or ON via its rule.
# This is necessary because mpbn follows Boolean rules during attractor search,
# so setting the starting state alone is not enough to hold a gene fixed.
bn_resting_dict = dict(bn_dict)
bn_resting_dict[source_gene] = "0"
bn_resting = mpbn.MPBooleanNetwork(bn_resting_dict)

bn_perturbed_dict = dict(bn_dict)
bn_perturbed_dict[source_gene] = "1"
bn_perturbed = mpbn.MPBooleanNetwork(bn_perturbed_dict)

# Step 1: find direct targets (genes whose Boolean rule mentions source_gene)
bn = mpbn.MPBooleanNetwork(bn_dict)
direct_targets = [gene for gene, func in bn.items() if re.search(r'\b' + re.escape(source_gene) + r'\b', str(func))]
 
# Step 2: find downstream activations (OFF->ON) from all-zero baseline
# Baseline: all genes OFF (no perturbation)
# Step 2: find downstream activations (OFF->ON) from all-zero baseline
# Baseline: all genes OFF (no perturbation)
baseline_state = {gene: 0 for gene in bn_dict}
baseline_active = set()
for attractor in bn.attractors(reachable_from=baseline_state):
    baseline_active.update(g for g, v in attractor.items() if v == 1)

# Perturbed: source gene = 1
initial_state = {gene: 0 for gene in bn_dict}
initial_state[source_gene] = 1

downstream_activated = set()
for attractor in bn.attractors(reachable_from=initial_state):
    downstream_activated.update(g for g, v in attractor.items() if v == 1 and g not in baseline_active)

# Compute attractors for conditional effects (Step 3) and context-dependent genes (Step 4).
# Both start from all-ones with source_gene locked OFF (resting) or ON (perturbed).
shared_start = {gene: 1 for gene in bn_dict}
resting_attractors   = list(bn_resting.attractors(reachable_from=shared_start))
perturbed_attractors = list(bn_perturbed.attractors(reachable_from=shared_start))
n_resting = len(resting_attractors)

# Print a header and summary line for the analysis output
print(f"\n  {source_gene} effect summary ({n_resting} resting states tested)")

#Step 1 Direct targets 
print(f"\n Direct targets of '{source_gene}': {len(direct_targets)} total")
for gene in sorted(direct_targets):
    print(f" {gene} [{gene_category.get(gene, '?')}]")
#Step 2 Downstream activations
print(f"\n Downstream activations (OFF->ON): {len(downstream_activated)} total")
for gene in sorted(downstream_activated):
    print(f" {gene} [{gene_category.get(gene, '?')}]")

# Step 3: downstream suppressions (ON->OFF) — mirror of Step 2
# Genes naturally ON in the baseline that turn OFF when source_gene is locked ON
perturbed_active = set()
for attractor in bn_perturbed.attractors(reachable_from=initial_state):
    perturbed_active.update(g for g, v in attractor.items() if v == 1)

downstream_suppressed = sorted(baseline_active - perturbed_active - {source_gene})

print(f"\n Step 3 - Downstream suppressions (ON->OFF): {len(downstream_suppressed)} total")
for gene in downstream_suppressed:
    print(f" {gene} [{gene_category.get(gene, '?')}]")

# Step 3 conditional: genes that shift in only some attractors, not all
conditional_suppressed = sorted(
    (g, sum(1 for a in resting_attractors if a.get(g, 0) == 1))
    for g in bn_dict
    if g not in downstream_suppressed and g != source_gene
    and any(a.get(g, 0) == 1 for a in resting_attractors)
    and any(a.get(g, 0) == 0 for a in perturbed_attractors)
    and not all(a.get(g, 0) == 1 for a in perturbed_attractors)
)
conditional_derepressed = sorted(
    (g, sum(1 for a in resting_attractors if a.get(g, 0) == 0))
    for g in bn_dict
    if g not in downstream_activated and g != source_gene
    and any(a.get(g, 0) == 0 for a in resting_attractors)
    and any(a.get(g, 0) == 1 for a in perturbed_attractors)
    and not all(a.get(g, 0) == 0 for a in perturbed_attractors)
)

# Print conditional effects only if any exist — these are nodes that shift state
# in some attractors but not all, meaning the effect of source_gene depends on
# which resting state the network was already in before the perturbation.
if conditional_suppressed or conditional_derepressed:
    print(f"\n  Conditional effects (only in some resting states):")
    for gene, count in conditional_derepressed:
        print(f"    ? ({count}/{n_resting}) + {gene}  [de-repressed]  [{gene_category.get(gene, '?')}]")
    for gene, count in conditional_suppressed:
        print(f"    ? ({count}/{n_resting}) - {gene}  [{gene_category.get(gene, '?')}]")

# Step 4: context-dependent genes — differ across perturbed attractors
# Genes with the same ON/OFF pattern always switch together (one independent decision)
gene_states = {}
for attractor in perturbed_attractors:
    for gene, val in attractor.items():
        gene_states.setdefault(gene, set()).add(val)
variable_genes = sorted(g for g, vals in gene_states.items() if len(vals) > 1)

patterns  = {g: tuple(a.get(g, 0) for a in perturbed_attractors) for g in variable_genes}
canonical = {g: p if p[0] == 0 else tuple(1 - v for v in p) for g, p in patterns.items()}  # normalise each pattern so it starts with 0, allowing genes that always flip together to share the same key
decisions = {}
for gene, key in canonical.items():
    decisions.setdefault(key, []).append(gene)

# Print Step 4: list the context-dependent genes grouped by decision.
# Each decision is a group of genes that always switch ON/OFF together across attractors,
# meaning they represent one independent regulatory choice rather than separate events.
print(f"\n Step 4 - Context-dependent genes: {len(variable_genes)} total, {len(decisions)} independent decision(s)")
for i, genes in enumerate(sorted(decisions.values(), key=lambda g: g[0]), 1):
    print(f"  Decision {i}:")
    for gene in sorted(genes):
        print(f" {gene} [{gene_category.get(gene, '?')}]")

print(f"\n{'='*80}")