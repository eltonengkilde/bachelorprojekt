# This script removes redundant edges from the filtered network CSV files.
# Redundant edges are rows where the same (source, relationship_category, target) triple
# appears more than once — this happens because the KG was built from 1000 papers,
# and the same biological relationship can be reported in multiple papers.
# Keeping duplicates does not add biological knowledge — it just inflates the network
# and can bias the Boolean model by making some edges appear artificially stronger.
# The unique edge (source -> relationship -> target) is kept; all extra copies are removed.
#
# NOTE: "At" prefix stripping (e.g. AtMYB46 -> MYB46) was investigated and rejected.
# AtMYB46 and MYB46 are treated as separate nodes in the KG — merging them caused
# nodes to disappear. The deduplication here only removes exact duplicate triples.

import csv
import os

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Files to deduplicate — add more if needed
FILES_TO_DEDUPLICATE = [
    ("filtered_large.csv",            "filtered_large_dedup.csv"),
    ("filtered_largerexperiment.csv", "filtered_largerexperiment_dedup.csv"),
]

for input_file, output_file in FILES_TO_DEDUPLICATE:
    input_path  = os.path.join(BASE_DIR, "networks_used_by_scripts", input_file)
    output_path = os.path.join(BASE_DIR, "networks_used_by_scripts", output_file)

    if not os.path.exists(input_path):
        print(f"Skipping {input_file} — file not found")
        continue

    with open(input_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = list(reader)

    # Deduplicate on the three columns that define a unique biological relationship.
    # All other columns (identifiers, full relationship string) are kept from the first
    # occurrence of each triple — they will be identical across duplicates anyway.
    seen = set()
    unique_rows = []
    for row in rows:
        key = (row["source"], row["relationship_category"], row["target"])
        if key not in seen:
            seen.add(key)
            unique_rows.append(row)

    duplicates_removed = len(rows) - len(unique_rows)

    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(unique_rows)

    print(f"\n{'='*60}")
    print(f"File:               {input_file}")
    print(f"Rows before:        {len(rows)}")
    print(f"Duplicates removed: {duplicates_removed}  ({100*duplicates_removed//len(rows)}%)")
    print(f"Rows after:         {len(unique_rows)}")
    print(f"Saved to:           {output_file}")

print("\nDone.")
