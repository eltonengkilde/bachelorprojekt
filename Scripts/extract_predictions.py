#!/opt/anaconda3/envs/bachelor_env/bin/python
"""
extract_predictions.py
======================
Reads one or more benchmark output text files, extracts all prediction
gene lists, and writes everything to a single CSV file that opens
directly in Excel.

HOW TO USE
----------
1. Add your benchmark output file paths to the FILES list below.
2. Set OUT_CSV to where you want the Excel-ready file saved.
3. Run:  python extract_predictions.py
4. Open the CSV in Excel.

Each row in the CSV is one predicted gene, with columns:
  Source_Gene | Hop | Solver | Prediction_Type | Gene | Category | Is_Gene
"""

import re, csv, os

# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURE HERE — add / remove files as needed
# ══════════════════════════════════════════════════════════════════════════════
BASE = "/Users/elton_1stboot/Documents/Bachelor Git/bachelorprojekt"

FILES = [
    # ── FORWARD files only (bench_fw_...) ─────────────────────────────────
    # For backward files use extract_predictions_bw.py instead.
    BASE + "/Benchmark synchronous simulation FW/bench_fw_MASTER_simulation_MYB46_hops1.txt",
    BASE + "/Benchmark synchronous simulation FW/bench_fw_MASTER_simulation_MYB46_hops2.txt",
    BASE + "/Benchmark synchronous simulation FW/bench_fw_MASTER_simulation_MYB46_hops3.txt",
    BASE + "/Benchmark synchronous simulation FW/bench_fw_MASTER_simulation_MYB46_hops4.txt",
    BASE + "/Benchmark synchronous simulation FW/bench_fw_MASTER_simulation_MYB46_hops5.txt",
    BASE + "/Benchmark synchronous simulation FW/bench_fw_MASTER_simulation_MYB46_hops6.txt",
    BASE + "/Benchmark synchronous simulation FW/bench_fw_MASTER_simulation_MYB46_hops7.txt",
    BASE + "/Benchmark synchronous simulation FW/bench_fw_MASTER_simulation_MYB46_hops8.txt",
    BASE + "/Benchmark synchronous simulation FW/bench_fw_MASTER_simulation_MYB46_hops9.txt",
    BASE + "/Benchmark synchronous simulation FW/bench_fw_MASTER_simulation_MYB46_hops10.txt",
    BASE + "/Benchmark synchronous simulation FW/bench_fw_MASTER_simulation_MYB46_hops11.txt",
    BASE + "/Benchmark synchronous simulation FW/bench_fw_MASTER_simulation_MYB46_hops12.txt",
    BASE + "/Benchmark synchronous simulation FW/bench_fw_MASTER_simulation_MYB46_hops13.txt",
    BASE + "/Benchmark synchronous simulation FW/bench_fw_MASTER_simulation_MYB46_hops14.txt",
    BASE + "/Benchmark synchronous simulation FW/bench_fw_MASTER_simulation_MYB46_hops15.txt",
    BASE + "/Benchmark synchronous simulation FW/bench_fw_MASTER_simulation_MYB46_hops16.txt",
    BASE + "/Benchmark synchronous simulation FW/bench_fw_MASTER_simulation_MYB46_hops17.txt",
    BASE + "/Benchmark synchronous simulation FW/bench_fw_MASTER_simulation_MYB46_hops18.txt",
    BASE + "/Benchmark synchronous simulation FW/bench_fw_MASTER_simulation_MYB46_hops19.txt",
    BASE + "/Benchmark synchronous simulation FW/bench_fw_MASTER_simulation_MYB46_hops20.txt",
    BASE + "/Benchmark synchronous simulation FW/bench_fw_MASTER_simulation_MYB46_hops21.txt",
    BASE + "/Benchmark synchronous simulation FW/bench_fw_MASTER_simulation_MYB46_hops22.txt",
    BASE + "/Benchmark synchronous simulation FW/bench_fw_MASTER_simulation_MYB46_hops23.txt",
    BASE + "/Benchmark synchronous simulation FW/bench_fw_MASTER_simulation_MYB46_hops24.txt",
    BASE + "/Benchmark synchronous simulation FW/bench_fw_MASTER_simulation_MYB46_hops25.txt",
    BASE + "/Benchmark synchronous simulation FW/bench_fw_MASTER_simulation_MYB46_hops26.txt",
    BASE + "/Benchmark synchronous simulation FW/bench_fw_MASTER_simulation_MYB46_hops27.txt",
    BASE + "/Benchmark synchronous simulation FW/bench_fw_MASTER_simulation_MYB46_hops28.txt",
    BASE + "/Benchmark synchronous simulation FW/bench_fw_MASTER_simulation_MYB46_hops29.txt",
    BASE + "/Benchmark synchronous simulation FW/bench_fw_MASTER_simulation_MYB46_hops30.txt",




    # Add more genes here when you run the benchmark for them, e.g.:
    # BASE + "/Benchmark BoNesis simulation FW/bench_fw_MASTER_larger_bonesis_HY5_hops2.txt",
    # BASE + "/Benchmark BoNesis simulation FW/bench_fw_MASTER_larger_bonesis_WRKY33_hops2.txt",
    # BASE + "/Benchmark BoNesis simulation FW/bench_fw_MASTER_larger_bonesis_EIN3_hops2.txt",
    # BASE + "/Benchmark BoNesis simulation FW/bench_fw_MASTER_larger_bonesis_pif4_hops2.txt",
]

# Where to save the Excel-ready CSV
OUT_CSV = BASE + "/Predictions_all_genes.csv"

# True = only keep [Gene] rows; False = include Pathways, Phenotypes etc.
GENES_ONLY = False  # False = include all entity types (Gene, Pathway, Phenotype, Metabolite)
                    # True  = keep only [Gene] rows (filter in script rather than Excel)
# ══════════════════════════════════════════════════════════════════════════════


# ── Regex patterns ────────────────────────────────────────────────────────────
RE_EXP_A      = re.compile(r'EXP A.*dark background', re.I)
RE_EXP_B      = re.compile(r'EXP B.*permissive background', re.I)
RE_ROBUST     = re.compile(r'ROBUST EFFECTS', re.I)
RE_NECESSARY  = re.compile(r'NECESSITY TEST', re.I)
RE_ACTIVATED  = re.compile(r'Activated \(OFF.*ON\)|Robustly activated', re.I)
RE_SUPPRESSED = re.compile(r'Suppressed \(ON.*OFF\)|Robustly suppressed', re.I)
RE_COND       = re.compile(r'Conditional effects', re.I)
RE_NEC_LIST   = re.compile(r'Necessary\s+\(|Necessary  \(', re.I)
RE_GENE_LINE  = re.compile(
    r'^\s{4}([+\-])\s+(?:\(\d+/\d+[^)]+\)\s+)?(\S+)\s+\[([^\]]+)\]'
)
RE_NEC_LINE   = re.compile(r'^\s{4}!\s+(\S+)\s+\[([^\]]+)\]')


def is_gene(cat_str):
    return 'gene' in {c.strip().lower() for c in cat_str.split(',')}


def parse_file(filepath, genes_only=True):
    """Return a list of dicts, one per predicted gene."""
    rows = []
    zone, subzone = None, None

    with open(filepath, encoding='utf-8', errors='ignore') as f:
        for line in f:
            if RE_EXP_A.search(line):     zone, subzone = 'exp_a', None; continue
            if RE_EXP_B.search(line):     zone, subzone = 'exp_b', None; continue
            if RE_ROBUST.search(line):    zone, subzone = 'robust', None; continue
            if RE_NECESSARY.search(line): zone, subzone = 'necessary', None; continue
            if RE_ACTIVATED.search(line):  subzone = 'activated'; continue
            if RE_SUPPRESSED.search(line): subzone = 'suppressed'; continue
            if RE_COND.search(line):       subzone = 'conditional'; continue
            if RE_NEC_LIST.search(line):   subzone = 'nec_list'; continue

            m = RE_GENE_LINE.match(line)
            if m:
                sign, gene, cat = m.group(1), m.group(2), m.group(3)
                if genes_only and not is_gene(cat): continue

                pred_type = None
                if zone == 'exp_a':
                    if subzone == 'activated'   and sign == '+': pred_type = 'Dark_Activated'
                    elif subzone == 'suppressed' and sign == '-': pred_type = 'Dark_Suppressed'
                elif zone == 'exp_b':
                    if subzone == 'activated'   and sign == '+': pred_type = 'Perm_Activated'
                    elif subzone == 'suppressed' and sign == '-': pred_type = 'Perm_Suppressed'
                    elif subzone == 'conditional':
                        pred_type = 'Conditional_Activated' if sign == '+' else 'Conditional_Suppressed'
                elif zone == 'robust':
                    if subzone == 'activated'   and sign == '+': pred_type = 'Robust_Activated'
                    elif subzone == 'suppressed' and sign == '-': pred_type = 'Robust_Suppressed'

                if pred_type:
                    rows.append({'pred_type': pred_type, 'gene': gene, 'category': cat})
                continue

            m2 = RE_NEC_LINE.match(line)
            if m2 and zone == 'necessary':
                gene, cat = m2.group(1), m2.group(2)
                if not genes_only or is_gene(cat):
                    rows.append({'pred_type': 'Necessary_Direct_Target', 'gene': gene, 'category': cat})

    return rows


def metadata_from_filename(filepath):
    """Extract source gene, hop, and solver from the filename."""
    name = os.path.basename(filepath)
    # Pattern: bench_fw_MASTER_larger_<solver>_<GENE>_hops<N>.txt
    m = re.search(r'_(bonesis|simulation)_([A-Za-z0-9]+)_hops(\d+)', name, re.I)
    if m:
        solver     = m.group(1).capitalize()   # Bonesis or Simulation
        source     = m.group(2)                # e.g. MYB46
        hop        = int(m.group(3))
    else:
        solver, source, hop = "Unknown", "Unknown", 0
    return source, hop, solver


def main():
    all_rows = []

    for filepath in FILES:
        if not os.path.exists(filepath):
            print(f"  SKIP (not found): {os.path.basename(filepath)}")
            continue

        source, hop, solver = metadata_from_filename(filepath)
        gene_rows = parse_file(filepath, genes_only=GENES_ONLY)

        for row in gene_rows:
            all_rows.append({
                'Source_Gene':     source,
                'Hop':             hop,
                'Solver':          solver,
                'Prediction_Type': row['pred_type'],
                'Gene':            row['gene'],
                'Category':        row['category'],
                'Is_Gene':         'Yes' if is_gene(row['category']) else 'No',
            })

        # Summary for this file
        counts = {}
        for row in gene_rows:
            counts[row['pred_type']] = counts.get(row['pred_type'], 0) + 1
        print(f"  {source}  hop {hop}  [{solver}]  →  "
              + "  ".join(f"{k.replace('_',' ')}: {v}" for k, v in sorted(counts.items())))

    if not all_rows:
        print("No data extracted — check that your file paths are correct.")
        return

    # Write CSV
    fieldnames = ['Source_Gene', 'Hop', 'Solver', 'Prediction_Type', 'Gene', 'Category', 'Is_Gene']
    with open(OUT_CSV, 'w', newline='', encoding='utf-8-sig') as f:
        # utf-8-sig adds the BOM that Excel needs to open UTF-8 correctly
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_rows)

    print(f"\n  Total rows written: {len(all_rows)}")
    print(f"  Saved to: {OUT_CSV}")
    print(f"\n  Open '{os.path.basename(OUT_CSV)}' in Excel.")
    print(f"  Use Data → Filter or a PivotTable to compare genes across hops and genes.")


if __name__ == '__main__':
    print(f"Extracting predictions from {len(FILES)} file(s)...\n")
    main()
