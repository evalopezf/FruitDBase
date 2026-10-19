"""
Filter the curated metadata by Kallisto pseudoalignment rate.

Reads `p_pseudoaligned` from each BioSample's Kallisto `run_info.json` and keeps
only the BioSamples at or above `min_pseudoalign_pct`.

Outputs:
    metadata_filtered.tsv  rows of the metadata whose BioSample passed
    removed_samples.tsv    BioSamples removed, with their rate and the reason
"""

import json
import os

import pandas as pd

metadata_path = snakemake.input.metadata
kallisto_dir = snakemake.params.kallisto_dir
min_pct = float(snakemake.params.min_pct)

df = pd.read_csv(metadata_path, sep="\t", encoding="latin1", dtype=str)

folder = df["BioSample"].str.strip() + "__" + df["tissue_name"].str.strip()


def pseudoaligned_pct(name):
    path = os.path.join(kallisto_dir, name, "run_info.json")
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return json.load(fh).get("p_pseudoaligned")


rates = {name: pseudoaligned_pct(name) for name in folder.unique()}
pct = folder.map(rates)

missing = pct.isna()
low = (~missing) & (pct < min_pct)
keep = ~(missing | low)

removed = pd.DataFrame(
    {
        "BioSample": df.loc[~keep, "BioSample"],
        "tissue_name": df.loc[~keep, "tissue_name"],
        "p_pseudoaligned": pct[~keep],
        "reason": ["no_run_info" if m else f"below_{min_pct:g}" for m in missing[~keep]],
    }
).drop_duplicates(subset="BioSample")

df[keep].to_csv(snakemake.output.metadata_filtered, sep="\t", index=False, encoding="latin1")
removed.to_csv(snakemake.output.removed, sep="\t", index=False)

n_total = df["BioSample"].nunique()
n_kept = df.loc[keep, "BioSample"].nunique()
print(f"BioSamples kept: {n_kept}/{n_total} (min pseudoalignment {min_pct:g}%)")
print(f"BioSamples removed: {len(removed)}")
