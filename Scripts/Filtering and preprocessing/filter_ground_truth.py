import pandas as pd
import os

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

INPUT_FILE  = os.path.join(BASE_DIR, "raw_data_used_by_scripts", "Ground_truth_experiment.csv")
OUTPUT_CSV  = os.path.join(BASE_DIR, "networks_used_by_scripts", "filtered_GT.csv")
OUTPUT_TXT  = os.path.join(BASE_DIR, "networks_used_by_scripts", "filtered_GT.txt")

KEEP = ["source", "source_category", "relationship_category", "target", "target_category"]

ACT = "Activation / Induction / Causation / Result"
SUP = "Repression / Inhibition / Negative Regulation"

# Ground truth uses free-text relationship descriptions rather than pre-labeled
# categories, so each synonym is mapped explicitly to the canonical label.
activations = {
    "activates", "activate", "activated", "activates expression of",
    "promote", "promotes expression of", "induce", "induces", "induced in",
    "enhance", "stimulates", "triggers", "cause", "drives", "facilitates",
    "upregulates", "up-regulates", "positively regulates",
    "increased expression of", "upregulated", "upregulated in",
    "up-regulated", "up-regulated in", "was up-regulated in", "AR-up-regulated",
}
suppressions = {
    "represses", "repress", "repressed", "inhibit", "inhibits",
    "suppresses", "suppress", "suppressed", "negatively regulates",
    "downregulates", "down-regulates", "downregulated", "downregulated in",
    "down-regulated", "down-regulated in", "antagonizes", "prevents",
    "AR-down-regulation",
}
KEEP_VALUES = {**{k: ACT for k in activations}, **{k: SUP for k in suppressions}}

df = pd.read_csv(INPUT_FILE, sep=";", encoding="utf-8-sig")
total_rows = len(df)

df = df.apply(lambda c: c.str.strip() if c.dtype == object else c)
df = df[df["relationship"].isin(KEEP_VALUES)].copy()
df["relationship_category"] = df["relationship"].map(KEEP_VALUES)
df = df[KEEP]

print(f"Rows before filtering: {total_rows}")
print(f"Rows after filtering:  {len(df)}")
print(df["relationship_category"].value_counts().to_string())

df.to_csv(OUTPUT_CSV, index=False)
df.to_csv(OUTPUT_TXT, sep="\t", index=False)

print(f"\nSaved CSV: {OUTPUT_CSV}")
print(f"Saved TXT: {OUTPUT_TXT}")
print("\nFirst 3 rows:")
print(df.head(3).to_string())
