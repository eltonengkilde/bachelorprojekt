import pandas as pd
import os

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

INPUT_FILE  = os.path.join(BASE_DIR, "raw_data_used_by_scripts", "oup_wly_elv_spr_sci_sciadv_all_kg.csv")
OUTPUT_CSV  = os.path.join(BASE_DIR, "networks_used_by_scripts", "filtered_networkL.csv")
OUTPUT_TXT  = os.path.join(BASE_DIR, "networks_used_by_scripts", "filtered_networkL.txt")

KEEP = ["source", "source_category", "relationship_category", "target", "target_category"]
ACT  = "Activation / Induction / Causation / Result"
SUP  = "Repression / Inhibition / Negative Regulation"

df = pd.read_csv(INPUT_FILE, encoding="utf-8")
total_rows = len(df)

df = df[df["relationship_category"].isin([ACT, SUP])].copy()
df = df[KEEP]

print(f"Rows before filtering: {total_rows}")
print(f"Rows after filtering:  {len(df)}")

df.to_csv(OUTPUT_CSV, index=False)
df.to_csv(OUTPUT_TXT, sep="\t", index=False)

print(f"Saved to {OUTPUT_CSV}")
print("\nFirst 3 rows:")
print(df.head(3).to_string())
