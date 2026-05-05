import pandas as pd
import os

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

FILES_TO_DEDUPLICATE = [
    ("filtered_networkL.csv", "filtered_networkL_dedup.csv"),
]

for input_file, output_file in FILES_TO_DEDUPLICATE:
    input_path  = os.path.join(BASE_DIR, "networks_used_by_scripts", input_file)
    output_path = os.path.join(BASE_DIR, "networks_used_by_scripts", output_file)
    txt_path    = output_path.replace(".csv", ".txt")

    if not os.path.exists(input_path):
        print(f"Skipping {input_file} — file not found")
        continue

    df = pd.read_csv(input_path)
    before = len(df)
    df = df.drop_duplicates(subset=["source", "relationship_category", "target"])
    removed = before - len(df)

    df.to_csv(output_path, index=False)
    df.to_csv(txt_path, sep="\t", index=False)

    print(f"\n{'='*60}")
    print(f"File:               {input_file}")
    print(f"Rows before:        {before}")
    print(f"Duplicates removed: {removed}  ({100*removed//before}%)")
    print(f"Rows after:         {len(df)}")
    print(f"Saved to:           {output_file}")
    print(f"Saved to:           {os.path.basename(txt_path)}")

print("\nDone.")
