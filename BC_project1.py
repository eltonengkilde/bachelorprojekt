# Cell A
# Import of standard CSV module for reading CSV files in python
# import os and BASE_DIR at the top — this always points to the project folder regardless of where VS Code runs Python from.
# Both file paths (kg_all_2005-2025_1000papers_combined.csv and network.txt) now use BASE_DIR so they always resolve correctly.
import csv
import os
import re
import mpbn

# Base directory — always the folder where this script lives, regardless of where VS Code launches from
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# opens the CSV file, handles special signs and ignores empty lines, it identifies and saves header lines and reads every line data as cells corresponding to the header title and does this for all the rows in a loop.
def new_func(BASE_DIR):
    with open(os.path.join(BASE_DIR, "kg_all_2005-2025_1000papers_combined.csv"), newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames
        rows = list(reader)
    return headers,rows

headers, rows = new_func(BASE_DIR)

# Verification
# print headers
print("Headers:")
print(headers)

# print first 3 rows
print("\nFirst 3 rows:")
for i, row in enumerate(rows[:3]):
    print(f"\nRow {i+1}:")
    print(row)

# Cell B
COLUMNS_TO_KEEP = [
    "source",
    "source_category",
    "source_identifier",
    "relationship",
    "relationship_category",
    "target",
    "target_category",
    "target_identifier"
]

filtered_rows = [{col: row[col] for col in COLUMNS_TO_KEEP} for row in rows]

#verification
print("\nFiltered rows (first 3):")
for row in filtered_rows[:3]:
    print(row)

# Cell C
# This block is for filtering rows for source_category and Target_category to only include "Gene / Protein"
# This chunk allows for filtering of specific values in the relationship_category column
# EDIT THIS LIST with the relationship values to keep
# (use the exact names printed in the cell above)

RELATIONSHIPS_CATEGORIES_TO_KEEP = [
    "Activation / Induction / Causation / Result",
    "Repression / Inhibition / Negative Regulation",
]

# Filter rows to only keep selected relationships
final_rows = [row for row in filtered_rows if row["relationship_category"] in RELATIONSHIPS_CATEGORIES_TO_KEEP]

# Verification
print(f"Rows before filtering: {len(filtered_rows)}")
print(f"Rows after filtering:  {len(final_rows)}")
print(f"Relationships kept: {RELATIONSHIPS_CATEGORIES_TO_KEEP}\n")

print("Filtered rows (first 5):")
for row in final_rows[:5]:
    print(row)

# Cell D
# This block is for filtering rows for source_category and Target_category to only include "Gene / Protein"

# This block is for filtering rows for source_category and target_category to only include "Gene / Protein"

SOURCE_CATEGORIES_TO_KEEP = [
    "Gene / Protein"
]

TARGET_CATEGORIES_TO_KEEP = [
    "Gene / Protein"
]

# Filter rows to only keep selected source_category AND target_category
last_rows = [
    row for row in final_rows
    if row["source_category"] in SOURCE_CATEGORIES_TO_KEEP
    and row["target_category"] in TARGET_CATEGORIES_TO_KEEP
]

# Verification
print(f"Rows before filtering: {len(final_rows)}")
print(f"Rows after filtering:  {len(last_rows)}")
print(f"Sources kept: {SOURCE_CATEGORIES_TO_KEEP}")
print(f"Targets kept: {TARGET_CATEGORIES_TO_KEEP}\n")

print("Filtered rows (first 5):")
for row in last_rows[:5]:
    print(row)

# Create a network file with the rows and colums picked out "last-rows" becomes a text file which can be used in cytoscape to visualize the network and also to check if the filtering is correct. this file is made by looping through the rows and joining the values of the columns with a tab and then writing these lines to a text file.
#_______________________________________________
save = []
for row in last_rows:
    save.append('\t'.join(list(row.values()))+'\n')
v = open(os.path.join(BASE_DIR, 'network.txt'),'w', encoding="utf-8")
v.writelines(save)
v.close()
# lets count the amount of activators and repressors

# boolean.py reserved words that cannot be used as gene names
BOOLEAN_RESERVED = {"TRUE", "FALSE", "NOT", "AND", "OR"}
# Helper: clean a name so mpbn can parse it
def clean_name(name):
    if not name or not isinstance(name, str):
        return "unknown"

    name = re.sub(r"[^a-zA-Z0-9_]", "_", name)
    name = re.sub(r"_+", "_", name)
    name = name.strip("_")

    if not name:
        return "unknown"
    elif name[0].isdigit():
        return "g" + name
    elif name.upper() in BOOLEAN_RESERVED:
        return "gene_" + name
    return name

# Cell H
# Group rows by target, collecting activators and suppressors
activators = {}
suppressors = {}

RELATIONSHIP_MAPPING = {
    "Activation / Induction / Causation / Result": "act",
    "Repression / Inhibition / Negative Regulation": "sup",
}

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

# Verification
print(f"Total nodes: {len(bn_dict)}")
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

print("\nExample activations (first 3 targets):")
for target, srcs in list(activators.items())[:3]:
    print(f"  {list(srcs)[:3]} --> {target}")

print("\nExample suppressions (first 3 targets):")
for target, srcs in list(suppressors.items())[:3]:
    print(f"  {list(srcs)[:3]} --| {target}")


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

# Cell J
# On the basis of mareks code this code aims to predict the immedeate downstream consequences of activating a specific gene from the KG
executor = input('type and enter "yolo" to run the mpbn attractor analysis')
if executor == "yolo":
    pass
else:
    print("Invalid input, exiting.")
    exit()

bn = mpbn.MPBooleanNetwork(bn_dict)

# Cell J
# Step 1: find direct targets (genes whose formula mentions source_gene)
direct_targets = [gene for gene, func in bn.items() if re.search(r'\b' + re.escape(source_gene) + r'\b', str(func))]

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

# Step 3: find suppressions (ON->OFF) via resting-state perturbation
resting_start     = {gene: (1 if gene in baseline_active else 0) for gene in bn_dict}
resting_attractors = list(bn.attractors(reachable_from=resting_start))

suppressed_in = {}  # gene -> count of resting attractors where it turns OFF (ON->OFF)
activated_in  = {}  # gene -> count of resting attractors where it turns ON  (OFF->ON, de-repression)
for r_att in resting_attractors:
    perturbed = {gene: (1 if val == 1 else 0) for gene, val in r_att.items()}
    perturbed[source_gene] = 1
    r_suppressed, r_activated = set(), set()
    for p_att in bn.attractors(reachable_from=perturbed):
        for gene in bn_dict:
            if r_att.get(gene, 0) == 1 and p_att.get(gene, 0) == 0:
                r_suppressed.add(gene)
            elif r_att.get(gene, 0) == 0 and p_att.get(gene, 0) == 1:
                r_activated.add(gene)
    for gene in r_suppressed:
        suppressed_in[gene] = suppressed_in.get(gene, 0) + 1
    for gene in r_activated:
        activated_in[gene]  = activated_in.get(gene, 0)  + 1

n_resting = len(resting_attractors)
consistent_suppressed = sorted(g for g, c in suppressed_in.items() if c == n_resting)
consistent_derepressed = sorted(g for g, c in activated_in.items() if c == n_resting and g not in downstream_activated)
conditional_suppressed = sorted((g, c) for g, c in suppressed_in.items() if c < n_resting)
conditional_derepressed = sorted((g, c) for g, c in activated_in.items() if c < n_resting and g not in downstream_activated)

# Step 4: variable genes (differ between resting attractors, unaffected by source gene)
all_perturbed_attractors = []
for r_att in resting_attractors:
    perturbed = {gene: (1 if val == 1 else 0) for gene, val in r_att.items()}
    perturbed[source_gene] = 1
    all_perturbed_attractors.extend(bn.attractors(reachable_from=perturbed))

gene_states = {}
for attractor in all_perturbed_attractors:
    for gene, val in attractor.items():
        if gene not in gene_states:
            gene_states[gene] = set()
        gene_states[gene].add(val)
variable_genes = sorted(g for g, vals in gene_states.items() if len(vals) > 1)

# Cluster variable genes by correlation — genes in the same decision always move together
# Build a pattern (tuple of 0/1) for each variable gene across all perturbed attractors
patterns = {g: tuple(att.get(g, 0) for att in all_perturbed_attractors) for g in variable_genes}
# Normalise: flip pattern if it starts with 1 so identical and opposite patterns map to the same key
canonical = {g: p if p[0] == 0 else tuple(1 - v for v in p) for g, p in patterns.items()}
# Group genes by canonical pattern
decisions = {}
for gene, key in canonical.items():
    if key not in decisions:
        decisions[key] = []
    decisions[key].append(gene)

# Print unified output
print(f"\n{'='*60}")
print(f"  {source_gene} effect summary ({n_resting} resting states tested)")

#Step 1
print(f"\n Direct targets of '{source_gene}': {len(direct_targets)} total")
for gene in sorted(direct_targets):
    print(f"    {gene}")
#Step 2
print(f"\n Downstream activations (OFF->ON): {len(downstream_activated)} total")
for gene in sorted(downstream_activated):
    print(f"  {gene}")
if consistent_derepressed:
    print(f"  (additional via de-repression: {len(consistent_derepressed)})")
    for gene in consistent_derepressed:
        print(f"  {gene}  [de-repressed]")

#Step 3
print(f"\n Downstream suppressions (ON->OFF): {len(consistent_suppressed)} total")
for gene in consistent_suppressed:
    print(f"  {gene}")

#Step 4
#Recording each variable gene's ON/OFF pattern across all perturbed attractors as a tuple e.g. (0,0,1,1,0,0,1,1)
#Normalising — flipping the pattern if it starts with 1, so genes that are always opposite still map to the same group key
#Grouping genes with identical normalised patterns into the same decision
print(f"\n Consistently variable genes (context-dependent): {len(variable_genes)} total, {len(decisions)} independent decision(s)")
for i, genes in enumerate(sorted(decisions.values(), key=lambda g: g[0]), 1):
    print(f"  Decision {i}:")
    for gene in sorted(genes):
        print(f"    {gene}")

if conditional_suppressed or conditional_derepressed:
    print(f"\n  Conditional effects (only in some resting states):")
    for gene, count in conditional_derepressed:
        print(f"    ? ({count}/{n_resting}) + {gene}  [de-repressed]")
    for gene, count in conditional_suppressed:
        print(f"    ? ({count}/{n_resting}) - {gene}")

indirect_targets = [gene for gene, func in bn.items() if source_gene in str(func) and gene not in direct_targets]
if len(indirect_targets) > 0:
    print(f"\n '{source_gene}' also appears in the formulas of {len(indirect_targets)} other genes (potential indirect effects):")
    for gene in sorted(indirect_targets):
        print(f"    {gene}")

print(f"\n{'='*60}")