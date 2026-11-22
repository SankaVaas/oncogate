import cellxgene_census

with cellxgene_census.open_soma() as census:
    obs = census["census_data"]["homo_sapiens"].obs.read(
        value_filter="tissue_general == 'breast' and is_primary_data == True",
        column_names=["dataset_id", "donor_id", "cell_type", "disease"],
    ).concat().to_pandas()

print("Cells per disease label:")
print(obs["disease"].value_counts(), "\n")

print("Donors per disease label:")
print(obs.groupby("disease")["donor_id"].nunique(), "\n")

print("Top cell types:")
print(obs["cell_type"].value_counts().head(15))