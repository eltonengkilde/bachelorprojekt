#!/opt/anaconda3/envs/bachelor_env/bin/python
"""
extract_predictions_bw.py
=========================
Backward-analysis version of extract_predictions.py.
Reads one or more backward benchmark output files and writes results
to the SAME CSV format as extract_predictions.py so both can be combined
in one Excel sheet.

HOW TO USE
----------
1. Add your backward benchmark output file paths to FILES below.
2. Run:  python extract_predictions_bw.py
3. The CSV appends to (or creates) Predictions_all_genes.csv — same file
   used by extract_predictions.py — so you get one combined sheet.

Prediction types written to CSV:
  Direct_Activator       — structural 1-hop activator (Step 1)
  Direct_Suppressor      — structural 1-hop suppressor (Step 1)
  Perm_ON_Dark_OFF       — stably ON when target ON, OFF when target OFF (Step 2 attractor trace)
  Perm_OFF_Dark_ON       — stably OFF when target ON, ON when target OFF (Step 2)
  Sufficient_Activator   — forcing gene ON causes target to turn ON (Step 3)
  Necessary_Activator    — KO turns target OFF (Step 4)
  Redundant_Activator    — KO keeps target ON (Step 4)
  Necessary_Suppressor   — forcing gene ON turns target OFF (Step 5)
  Suppressor_Release     — KO of gene turns target ON (Step 5)
"""

import re, csv, os

# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURE HERE
# ══════════════════════════════════════════════════════════════════════════════
BASE = "/Users/elton_1stboot/Documents/Bachelor Git/bachelorprojekt"

FILES = [
    BASE + "/Benchmark BoNesis simulation BW/bench_bw_MASTER_larger_bonesis_MYB46_hops1.txt",
    BASE + "/Benchmark BoNesis simulation BW/bench_bw_MASTER_larger_bonesis_MYB46_hops2.txt",

    
    
    
    # Add more files here, e.g.:
    # BASE + "/Benchmark synchronous simulation BW/bench_bw_MASTER_larger_simulation_MYB46_hops3.txt",
]

# Same output file as extract_predictions.py — append so both scripts feed one sheet.
OUT_CSV    = BASE + "/Predictions_all_genes.csv"
GENES_ONLY = False  # False = include all entity types (Gene, Pathway, Phenotype, Metabolite)
                    # True  = keep only [Gene] rows (filter in script rather than Excel)

# Set to True to overwrite the CSV; False to append to an existing one.
OVERWRITE  = False
# ══════════════════════════════════════════════════════════════════════════════


def is_gene(cat_str):
    return 'gene' in {c.strip().lower() for c in cat_str.split(',')}


# ── Section header patterns ───────────────────────────────────────────────────
_STEP1     = re.compile(r'Step 1.*Structural direct regulator', re.I)
_STEP1_ACT = re.compile(r'Activators:\s*\d+', re.I)
_STEP1_SUP = re.compile(r'Suppressors:\s*\d+', re.I)
_STEP2     = re.compile(r'Step 2.*Attractor state', re.I)
_STEP3     = re.compile(r'Step 3.*Sufficient upstream', re.I)
_STEP4     = re.compile(r'Step 4.*Necessity.*activator', re.I)
_STEP5     = re.compile(r'Step 5.*Suppressor', re.I)

# Step 1 & 2: 6-space indent, optional second bracket for Step 2
# "      + gene   [Category]"
# "      + gene   [stably ON perm / OFF dark]  [Category]"
_GENE6 = re.compile(
    r'^\s{6}([+\-])\s+(\S+)\s+\[([^\]]+)\](?:\s+\[([^\]]+)\])?'
)

# Steps 3-5: 4-space indent, signs can be +, -, !, ~
# "    ! gene   [required activator]  [Category]"
_GENE4 = re.compile(
    r'^\s{4}([+\-!~])\s+(\S+)\s+\[([^\]]+)\](?:\s+\[([^\]]+)\])?'
)


def parse_file(filepath, genes_only=True):
    """Return list of (pred_type, gene, category) tuples."""
    rows = []
    zone = subzone = None

    with open(filepath, encoding='utf-8', errors='ignore') as f:
        for line in f:
            # Section detection
            if _STEP1.search(line):     zone, subzone = 'step1', None; continue
            if _STEP2.search(line):     zone, subzone = 'step2', None; continue
            if _STEP3.search(line):     zone, subzone = 'step3', None; continue
            if _STEP4.search(line):     zone, subzone = 'step4', None; continue
            if _STEP5.search(line):     zone, subzone = 'step5', None; continue
            if _STEP1_ACT.search(line) and zone == 'step1': subzone = 'act'; continue
            if _STEP1_SUP.search(line) and zone == 'step1': subzone = 'sup'; continue

            # 6-space gene lines (Steps 1 and 2)
            m = _GENE6.match(line)
            if m and zone in ('step1', 'step2'):
                sign, gene, bkt1, bkt2 = m.group(1), m.group(2), m.group(3), m.group(4)
                cat = bkt2 if bkt2 else bkt1   # category is the last bracket
                if genes_only and not is_gene(cat): continue
                pred = None
                if zone == 'step1':
                    if   subzone == 'act' and sign == '+': pred = 'Direct_Activator'
                    elif subzone == 'sup' and sign == '-': pred = 'Direct_Suppressor'
                elif zone == 'step2':
                    if   sign == '+': pred = 'Perm_ON_Dark_OFF'
                    elif sign == '-': pred = 'Perm_OFF_Dark_ON'
                if pred:
                    rows.append((pred, gene, cat))
                continue

            # 4-space gene lines (Steps 3–5)
            m2 = _GENE4.match(line)
            if m2 and zone in ('step3', 'step4', 'step5'):
                sign, gene, bkt1, bkt2 = m2.group(1), m2.group(2), m2.group(3), m2.group(4)
                cat = bkt2 if bkt2 else bkt1
                if genes_only and not is_gene(cat): continue
                pred = None
                if   zone == 'step3' and sign == '+':  pred = 'Sufficient_Activator'
                elif zone == 'step4' and sign == '!':  pred = 'Necessary_Activator'
                elif zone == 'step4' and sign == '~':  pred = 'Redundant_Activator'
                elif zone == 'step5' and sign == '!':  pred = 'Necessary_Suppressor'
                elif zone == 'step5' and sign == '~':  pred = 'Suppressor_Release'
                if pred:
                    rows.append((pred, gene, cat))

    return rows


def metadata_from_filename(filepath):
    name = os.path.basename(filepath)
    m = re.search(r'_(bonesis|simulation)_([A-Za-z0-9]+)_hops(\d+)', name, re.I)
    if m:
        solver = m.group(1).capitalize()
        source = m.group(2)
        hop    = int(m.group(3))
    else:
        solver, source, hop = 'Unknown', 'Unknown', 0
    return source, hop, solver


def main():
    all_rows = []

    for filepath in FILES:
        if not os.path.exists(filepath):
            print(f"  SKIP (not found): {os.path.basename(filepath)}")
            continue

        source, hop, solver = metadata_from_filename(filepath)
        gene_rows = parse_file(filepath, genes_only=GENES_ONLY)

        for pred, gene, cat in gene_rows:
            all_rows.append({
                'Direction':       'Backward',
                'Source_Gene':     source,
                'Hop':             hop,
                'Solver':          solver,
                'Prediction_Type': pred,
                'Gene':            gene,
                'Category':        cat,
                'Is_Gene':         'Yes' if is_gene(cat) else 'No',
            })

        counts = {}
        for pred, _, _ in gene_rows:
            counts[pred] = counts.get(pred, 0) + 1

        if counts:
            print(f"  [BW] {source}  hop {hop}  [{solver}]  →  "
                  + "  ".join(f"{k.replace('_',' ')}: {v}" for k, v in sorted(counts.items())))
        else:
            print(f"  [BW] {source}  hop {hop}  [{solver}]  →  "
                  f"(no prediction sections — file may have timed out before output was written)")

    if not all_rows:
        print("\nNo data extracted. Check file paths are correct.")
        return

    fieldnames = ['Direction', 'Source_Gene', 'Hop', 'Solver',
                  'Prediction_Type', 'Gene', 'Category', 'Is_Gene']

    # Decide whether to write a fresh CSV or append to an existing one
    file_exists = os.path.exists(OUT_CSV)
    mode = 'w' if OVERWRITE or not file_exists else 'a'
    write_header = (mode == 'w')

    with open(OUT_CSV, mode, newline='', encoding='utf-8-sig') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
        writer.writerows(all_rows)

    action = "Written" if write_header else "Appended"
    print(f"\n  {action} {len(all_rows)} rows to: {OUT_CSV}")
    print(f"  Open '{os.path.basename(OUT_CSV)}' in Excel — Forward and Backward results in one sheet.")
    print(f"  Filter by 'Direction' column to separate Forward from Backward.")


if __name__ == '__main__':
    print(f"Extracting BACKWARD predictions from {len(FILES)} file(s)...\n")
    main()
