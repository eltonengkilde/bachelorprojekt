import csv
import os
import re
# Base directory — always the folder where this script lives, regardless of where VS Code launches from
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Cell B – filter criteria (defined before get_data so they can be passed in)
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

# Cell C – relationship categories to keep
RELATIONSHIPS_CATEGORIES_TO_KEEP = [
    "Activation / Induction / Causation / Result",
    "Repression / Inhibition / Negative Regulation",
]

# Cell D – source/target categories to keep
SOURCE_CATEGORIES_TO_KEEP = [
    "Gene / Protein"
]

TARGET_CATEGORIES_TO_KEEP = [
    "Gene / Protein"
]

OUTPUT_FILE = os.path.join(BASE_DIR, "filtered_dummy.csv")

def get_data(BASE_DIR, columns, rel_categories, src_categories, tgt_categories):
    print("Reading full source file")
    rel_set = set(rel_categories)
    src_set = set(src_categories)
    tgt_set = set(tgt_categories)
    rows = []
    with open(os.path.join(BASE_DIR, "oup_wly_elv_spr_sci_sciadv_all_kg.csv"), newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames
        for row in reader:
            if (row["relationship_category"] in rel_set
                    and row["source_category"] in src_set
                    and row["target_category"] in tgt_set):
                rows.append({col: row[col] for col in columns})

    with open(OUTPUT_FILE, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved to {OUTPUT_FILE}")
    return headers, rows

headers, last_rows = get_data(
    BASE_DIR,
    COLUMNS_TO_KEEP,
    RELATIONSHIPS_CATEGORIES_TO_KEEP,
    SOURCE_CATEGORIES_TO_KEEP,
    TARGET_CATEGORIES_TO_KEEP,
)

# Verification
print("Headers:")
print(headers)

print(f"\nRows after all filters: {len(last_rows)}")
print(f"Relationships kept: {RELATIONSHIPS_CATEGORIES_TO_KEEP}")
print(f"Sources kept: {SOURCE_CATEGORIES_TO_KEEP}")
print(f"Targets kept: {TARGET_CATEGORIES_TO_KEEP}\n")

print("First 3 rows:")
for i, row in enumerate(last_rows[:3]):
    print(f"\nRow {i+1}:")
    print(row)

#Create a Compacted .CSV file for the analysis code
save = []
for row in last_rows:
    save.append('\t'.join(list(row.values()))+'\n')
v = open(os.path.join(BASE_DIR, 'filtere_dummy_network.txt'),'w', encoding="utf-8")
v.writelines(save)
v.close()