import numpy as np, pandas as pd, cellxgene_census
from config import *
from common import load_surface_genes, META_COLS

surface = load_surface_genes()
print("surface genes in list:", len(surface))

with cellxgene_census.open_soma() as census:
    human = census["census_data"]["homo_sapiens"]
    cols = ["soma_joinid", "dataset_id", "donor_id", "tissue_general", "cell_type", "disease"]
    read = lambda f: human.obs.read(value_filter=f, column_names=cols).concat().to_pandas()

    tumor = read(TUMOR_FILTER)
    print("tumor cells found:", len(tumor))
    tumor = tumor.sample(min(len(tumor), MAX_TUMOR_CELLS), random_state=SEED)
    tumor["label"] = "tumor"

    parts = []
    for t in HEALTHY_TISSUES:
        df = read(f"{HEALTHY_FILTER} and tissue_general == '{t}'")
        if len(df) == 0:
            print("  no healthy data for:", t); continue
        print(f"  {t}: {len(df)} cells")
        parts.append(df.sample(min(len(df), HEALTHY_PER_TISSUE), random_state=SEED))
    healthy = pd.concat(parts); healthy["label"] = "healthy"
    meta = pd.concat([tumor, healthy], ignore_index=True)

    var = human.ms["RNA"].var.read(column_names=["soma_joinid", "feature_name"]).concat().to_pandas()
    var = var[var.feature_name.isin(surface)].drop_duplicates("feature_name")
    print("surface genes found in Census:", len(var))

    adata = cellxgene_census.get_anndata(
        census, organism="Homo sapiens",
        obs_coords=np.sort(meta.soma_joinid.unique()),
        var_coords=var.soma_joinid.values,
        column_names={"obs": ["soma_joinid"]},
    )

# align expression rows with our metadata
order = pd.Series(range(adata.n_obs), index=adata.obs["soma_joinid"].values)
X = adata.X[order[meta.soma_joinid].values]
pos = pd.DataFrame((X >= MIN_COUNT).toarray(), columns=adata.var["feature_name"].values)
pos = pos.loc[:, ~pos.columns.duplicated()]

# split by DONOR (patient), never by cell
donor_key = meta.dataset_id + "|" + meta.donor_id
donors = donor_key.unique()
test_donors = set(np.random.default_rng(SEED).permutation(donors)[: int(len(donors) * TEST_FRACTION)])
meta["split"] = np.where(donor_key.isin(test_donors), "test", "train")

# shortlist genes using TRAIN tumor cells only
is_tr_tumor = ((meta.split == "train") & (meta.label == "tumor")).values
keep = pos.columns[pos[is_tr_tumor].mean() >= MIN_TUMOR_PREVALENCE]
print("genes kept (prevalent in train tumor):", len(keep))

out = pd.concat([meta[META_COLS], pos[keep]], axis=1)
out.to_parquet(DATA / "cells.parquet")
print(out.groupby(["split", "label"]).size())
print("train/test donors:", (len(donors) - len(test_donors)), "/", len(test_donors))