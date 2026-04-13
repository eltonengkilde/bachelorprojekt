# Import csv for reading and writing CSV files, os for building file paths
import csv
import os

# BASE_DIR points two levels up from this script (Scripts/ -> bachelorprojekt/)
# so all file paths resolve correctly regardless of where VS Code runs Python from
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# The 8 columns to keep from the raw KG — everything else in each row is discarded.
# This covers the source node, the relationship, and the target node, each with category and identifier.
COLUMNS_TO_KEEP = [
    "source", "source_category", "relationship_category",
    "target", "target_category"
]

# Read the full raw knowledge graph CSV row by row.
# Unlike the simple filter (CSVfiltering_simple.py), this version does NOT restrict
# source/target category — it keeps all entity types (genes, phenotypes, metabolites, etc.)
# as long as the relationship is an activation or repression.
# total_rows counts every row seen so we can report how many were filtered out.
rows = []
total_rows = 0
with open(os.path.join(BASE_DIR, "raw_data_used_by_scripts", "oup_wly_elv_spr_sci_sciadv_all_kg.csv"), newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    for row in reader:
        total_rows += 1
        # Keep only rows where the relationship is an activation or a repression.
        # Rows with any other relationship type (binding, transport, correlation, etc.) are skipped.
        # Each surviving row is trimmed to only the 8 columns defined above.
        if (row["relationship_category"] in ("Activation / Induction / Causation / Result",
                                              "Repression / Inhibition / Negative Regulation")):
            rows.append({col: row[col] for col in COLUMNS_TO_KEEP})

print(f"Rows before filtering: {total_rows}")

# Write the filtered rows to a CSV file in the networks folder.
# This is the file loaded by BC_project_fw_largerexperiment.py and BC_project_bw_largerexperiment.py.
OUTPUT_FILE = os.path.join(BASE_DIR, "networks_used_by_scripts", "filtered_largerexperiment.csv")
with open(OUTPUT_FILE, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=COLUMNS_TO_KEEP)
    writer.writeheader()
    writer.writerows(rows)

print(f"Rows after filtering: {len(rows)}")
print(f"Saved to {OUTPUT_FILE}")

# Print the first 3 rows as a sanity check so the user can verify the output looks correct.
print("\nFirst 3 rows:")
for i, row in enumerate(rows[:3]):
    print(f"\nRow {i+1}:")
    print(row)

# Also save a tab-separated version of the same data as a .txt network file.
# Some graph analysis tools expect tab-separated input rather than comma-separated.
with open(os.path.join(BASE_DIR, "networks_used_by_scripts", "filtered_largerexperiment_network.txt"), "w", encoding="utf-8") as f:
    for row in rows:
        f.write("\t".join(row.values()) + "\n")
