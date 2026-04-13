import csv
import os

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

COLUMNS_TO_KEEP = [
    "source", "source_category", "source_identifier",
    "relationship", "relationship_category",
    "target", "target_category", "target_identifier"
]

rows = []
total_rows = 0
with open(os.path.join(BASE_DIR, "raw_data_used_by_scripts", "kg_all_2005-2025_1000papers_combined.csv"), newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    for row in reader:
        total_rows += 1
        if (row["relationship_category"] in ("Activation / Induction / Causation / Result",
                                              "Repression / Inhibition / Negative Regulation")
                and row["source_category"] == "Gene / Protein"
                and row["target_category"] == "Gene / Protein"):
            rows.append({col: row[col] for col in COLUMNS_TO_KEEP})

print(f"Rows before filtering: {total_rows}")

OUTPUT_FILE = os.path.join(BASE_DIR, "networks_used_by_scripts", "filtered_small.csv")
with open(OUTPUT_FILE, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=COLUMNS_TO_KEEP)
    writer.writeheader()
    writer.writerows(rows)

print(f"Rows after filtering: {len(rows)}")
print(f"Saved to {OUTPUT_FILE}")
print("\nFirst 3 rows:")
for i, row in enumerate(rows[:3]):
    print(f"\nRow {i+1}:")
    print(row)

# Save as tab-separated network file for analysis
with open(os.path.join(BASE_DIR, "networks_used_by_scripts", "filtered_small_network.txt"), "w", encoding="utf-8") as f:
    for row in rows:
        f.write("\t".join(row.values()) + "\n")
