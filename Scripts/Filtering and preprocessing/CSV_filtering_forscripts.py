import pandas as pd
import os

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

INPUT_FILE  = os.path.join(BASE_DIR, "raw_data_used_by_scripts", "oup_wly_elv_spr_sci_sciadv_all_kg.csv")
OUTPUT_CSV  = os.path.join(BASE_DIR, "networks_used_by_scripts", "filtered_networkL.csv")
OUTPUT_TXT  = os.path.join(BASE_DIR, "networks_used_by_scripts", "filtered_networkL.txt")

KEEP = ["source", "source_category", "relationship_category", "target", "target_category"]

ACT = "Activation / Induction / Causation / Result"
SUP = "Repression / Inhibition / Negative Regulation"
KEEP_VALUES = {ACT, SUP}

# Large file (11M rows) — read and filter in chunks to avoid loading it all into memory.
chunks = []
total_rows = 0
for chunk in pd.read_csv(INPUT_FILE, encoding="utf-8", chunksize=100_000, low_memory=False):
    total_rows += len(chunk)
    chunk = chunk.apply(lambda c: c.str.strip() if c.dtype == object else c)
    filtered = chunk[chunk["relationship_category"].isin(KEEP_VALUES)][KEEP].copy()
    if len(filtered):
        chunks.append(filtered)

df = pd.concat(chunks, ignore_index=True)

print(f"Rows before filtering: {total_rows}")
print(f"Rows after filtering:  {len(df)}")
print(df["relationship_category"].value_counts().to_string())

df.to_csv(OUTPUT_CSV, index=False)
df.to_csv(OUTPUT_TXT, sep="\t", index=False)

print(f"\nSaved CSV: {OUTPUT_CSV}")
print(f"Saved TXT: {OUTPUT_TXT}")
print("\nFirst 3 rows:")
print(df.head(3).to_string())
