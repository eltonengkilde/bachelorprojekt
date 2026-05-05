import pandas as pd
import os
import re

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

INPUT_FILE  = os.path.join(BASE_DIR, "networks_used_by_scripts", "filtered_networkL.csv")
OUTPUT_CSV  = os.path.join(BASE_DIR, "networks_used_by_scripts", "filtered_networkL_normalized.csv")
OUTPUT_TXT  = os.path.join(BASE_DIR, "networks_used_by_scripts", "filtered_networkL_normalized.txt")

df = pd.read_csv(INPUT_FILE).dropna(subset=["source", "target"])
df["source"] = df["source"].astype(str)
df["target"] = df["target"].astype(str)
input_rows = len(df)
print(f"Loaded {input_rows} rows from {os.path.basename(INPUT_FILE)}")

# ── Build alias map ───────────────────────────────────────────────────────────
# Regex alias detection cannot be vectorised — it requires examining each node
# name individually and checking membership in a set.  The loop stays in pure
# Python; only the application of the resulting map moves to pandas.
all_names = set(df["source"]) | set(df["target"])
name_lookup = {n.upper(): n for n in all_names}

SKIP = {
    "AND", "OR", "THE", "FOR", "OF", "IN", "TO", "AT",
    "GENE", "GENES", "PROTEIN", "PROTEINS", "FACTOR", "FACTORS",
    "FAMILY", "COMPLEX", "DOMAIN", "LIKE", "TYPE", "CLASS",
    "PATHWAY", "PROCESS", "RESPONSE", "SIGNALING", "REGULATION",
    "SCW", "SWN", "SWNS", "TF", "TFS",
}

alias_map = {}
for name in sorted(all_names):
    m = re.search(r'\(([^)]+)\)\s*$', name)
    if not m:
        continue
    abbrev = m.group(1).strip()
    if not re.match(r'^[A-Za-z][A-Za-z0-9]{1,11}$', abbrev):
        continue
    if abbrev.upper() in SKIP:
        continue
    if abbrev.upper() in name_lookup and name_lookup[abbrev.upper()] != name:
        alias_map[name] = name_lookup[abbrev.upper()]

print(f"Alias pairs detected: {len(alias_map)}")
print(f"\nSample aliases (first 20):")
for long, short in list(alias_map.items())[:20]:
    print(f"  {long!r:70s} -> {short!r}")

# ── Apply alias map, remove self-loops, re-deduplicate ────────────────────────
replacements = (df["source"].isin(alias_map)).sum() + (df["target"].isin(alias_map)).sum()
df["source"] = df["source"].replace(alias_map)
df["target"] = df["target"].replace(alias_map)

before_selfloop = len(df)
df = df[df["source"] != df["target"]]
self_loops_removed = before_selfloop - len(df)

before_dedup = len(df)
df = df.drop_duplicates(subset=["source", "relationship_category", "target"])
new_duplicates = before_dedup - len(df)

print(f"\nCell replacements made: {replacements}")
print(f"Self-loops removed: {self_loops_removed}")
print(f"New duplicates created by alias merging (removed): {new_duplicates}")
print(f"Final row count: {len(df)}")

# ── Save ──────────────────────────────────────────────────────────────────────
df.to_csv(OUTPUT_CSV, index=False)
df.to_csv(OUTPUT_TXT, sep="\t", index=False)

print(f"\nSaved CSV: {os.path.basename(OUTPUT_CSV)}")
print(f"Saved TXT: {os.path.basename(OUTPUT_TXT)}")
print(f"\n{'='*60}")
print(f"Input rows:              {input_rows}")
print(f"Alias pairs merged:      {len(alias_map)}")
print(f"Cell replacements:       {replacements}")
print(f"Self-loops removed:      {self_loops_removed}")
print(f"New duplicates removed:  {new_duplicates}")
print(f"Output rows:             {len(df)}")
print(f"Net reduction:           {input_rows - len(df)} rows")
