import pandas as pd
import os
import re

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

INPUT_FILE = os.path.join(BASE_DIR, "networks_used_by_scripts", "filtered_GT.csv")
OUTPUT_CSV = os.path.join(BASE_DIR, "networks_used_by_scripts", "filtered_GT_normalized.csv")
OUTPUT_TXT = os.path.join(BASE_DIR, "networks_used_by_scripts", "filtered_GT_normalized.txt")

DEDUP = ["source", "relationship_category", "target"]

# ── Detection rules ───────────────────────────────────────────────────────────

SKIP_WORDS = {
    "AND", "OR", "THE", "FOR", "OF", "IN", "TO", "AT",
    "GENE", "GENES", "PROTEIN", "PROTEINS", "FACTOR", "FACTORS",
    "FAMILY", "COMPLEX", "DOMAIN", "LIKE", "TYPE", "CLASS",
    "PATHWAY", "PROCESS", "RESPONSE", "SIGNALING", "REGULATION",
    "SCW", "SWN", "SWNS", "TF", "TFS",
}

AT_PREFIX   = re.compile(r'^At([A-Za-z].*)$')
PAREN_RE    = re.compile(r'\(([^)]+)\)\s*$')
GENE_SYM_RE = re.compile(r'^[A-Za-z][A-Za-z0-9]{1,11}$')

SUFFIX_RE = re.compile(r"""(
      \s+protein\(s\)\s+activity | \s+protein\s+activity
    | \s+expression\s+levels?    | \s+transcript\s+levels?
    | \s+mrna\s+levels?          | \s+promoter\s+activity
    | \s+overexpression          | \s+expression
    | \s+promoter                | \s+function
    | \s+activity                | \s+pathway
    | \s+proteins?               | \s+transcripts?
    | \s+mrna                    | \s+levels?
)\s*$""", re.IGNORECASE | re.VERBOSE)

# Last words that can appear at the end of a sub-entity name — used as a fast
# pre-filter before running the full SUFFIX_RE on each name.
SUFFIX_LAST_WORDS = frozenset({
    "activity", "expression", "promoter", "function", "pathway",
    "protein", "proteins", "overexpression", "transcript", "transcripts",
    "mrna", "level", "levels",
})


def build_alias_map(df):
    """Scan every node name and return a dict of long→canonical for all three
    detection strategies: parenthesis abbreviation, At-prefix, and sub-entity suffix."""
    # pd.unique on a concatenated Series is 5-10× faster than Python set() on
    # object-dtype columns because it uses a Cython hash table internally.
    names  = pd.unique(pd.concat([df["source"], df["target"]]))
    lookup = {n.upper(): n for n in names}
    alias_map = {}

    for name in names:
        # Strategy 1: "FULL NAME (ABBREV)" → ABBREV
        if "(" in name:
            m = PAREN_RE.search(name)
            if m:
                abbrev = m.group(1).strip()
                if (GENE_SYM_RE.match(abbrev)
                        and abbrev.upper() not in SKIP_WORDS
                        and abbrev.upper() in lookup
                        and lookup[abbrev.upper()] != name):
                    alias_map[name] = lookup[abbrev.upper()]
                    continue

        # Strategy 2: AtGENE → GENE when bare name exists
        if name.startswith("At") and len(name) > 2 and name[2].isalpha():
            m = AT_PREFIX.match(name)
            if m:
                bare = m.group(1)
                if bare.upper() in lookup and lookup[bare.upper()] != name:
                    alias_map[name] = lookup[bare.upper()]
                    continue

        # Strategy 3: "GENE suffix" → GENE when bare name exists
        # Two-level pre-filter: space present AND last word is a known suffix.
        # This avoids running SUFFIX_RE on the vast majority of node names.
        if " " in name:
            last_word = name.rsplit(None, 1)[-1].lower().rstrip("s")
            if last_word in SUFFIX_LAST_WORDS or (last_word + "s") in SUFFIX_LAST_WORDS:
                stripped = SUFFIX_RE.sub("", name).strip()
                if stripped and stripped != name and stripped.upper() in lookup and lookup[stripped.upper()] != name:
                    alias_map[name] = lookup[stripped.upper()]

    return alias_map


def apply_and_clean(df, alias_map):
    df = df.copy()
    # Mask-based replacement is faster than .replace(dict) on object-dtype
    # columns: isin() finds the rows to change, map() applies only to those.
    src_mask = df["source"].isin(alias_map)
    tgt_mask = df["target"].isin(alias_map)
    df["source"] = df["source"].where(~src_mask, df["source"].map(alias_map))
    df["target"] = df["target"].where(~tgt_mask, df["target"].map(alias_map))
    df = df[df["source"] != df["target"]]
    return df.drop_duplicates(subset=DEDUP)


# ── Load + initial dedup ──────────────────────────────────────────────────────
df = pd.read_csv(INPUT_FILE)
rows_raw = len(df)
df = df.drop_duplicates(subset=DEDUP)
print(f"Loaded {rows_raw} rows from {os.path.basename(INPUT_FILE)}", flush=True)
print(f"Initial dedup: {len(df)} rows kept\n", flush=True)

# ── Alias normalisation — iterate until convergence ───────────────────────────
# Each iteration re-scans the current node set, so chained aliases resolve
# naturally (e.g. "AtABI3 expression" → "ABI3 expression" → "ABI3").
iteration = 0
while True:
    alias_map = build_alias_map(df)
    if not alias_map:
        break
    replacements = (df["source"].isin(alias_map)).sum() + (df["target"].isin(alias_map)).sum()
    df = apply_and_clean(df, alias_map)
    iteration += 1
    print(f"Iteration {iteration}: {len(alias_map):4d} aliases detected, "
          f"{replacements:5d} replacements → {len(df)} rows", flush=True)

print(f"\nConverged after {iteration} iteration(s)")

# ── Case-variant collapsing ───────────────────────────────────────────────────
# Merge names that are identical when uppercased (e.g. "MYB46" / "Myb46").
# Pick the most-frequent variant; ties prefer all-caps.
freq = pd.concat([df["source"], df["target"]]).value_counts()
names_df = pd.DataFrame({"name": list(set(df["source"]) | set(df["target"]))})
names_df["upper"]    = names_df["name"].str.upper()
names_df["freq"]     = names_df["name"].map(freq).fillna(0)
names_df["is_upper"] = names_df["name"] == names_df["name"].str.upper()
canonical = (names_df.sort_values(["freq", "is_upper"], ascending=[False, False])
             .groupby("upper")["name"].first())
case_map = {v: canonical[v.upper()] for v in names_df["name"] if canonical.get(v.upper(), v) != v}

if case_map:
    df = apply_and_clean(df, case_map)
    print(f"Case variants collapsed: {len(case_map)}")

# ── Save ──────────────────────────────────────────────────────────────────────
df.to_csv(OUTPUT_CSV, index=False)
df.to_csv(OUTPUT_TXT, sep="\t", index=False)

print(f"\nSaved CSV: {os.path.basename(OUTPUT_CSV)}")
print(f"Saved TXT: {os.path.basename(OUTPUT_TXT)}")
print(f"\n{'='*60}")
print(f"Input rows:   {rows_raw}")
print(f"Output rows:  {len(df)}")