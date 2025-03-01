import pandas as pd
from anndata import AnnData
import os
import muon as mu
import numpy as np


def load_cll(data_path="data/"):
    """Load and return CLL dataset."""
    data_dict = {}

    obs = pd.read_csv(os.path.join(data_path, "cll_metadata.csv"), index_col="Sample")
    obs["Gender"] = obs["Gender"].replace({"m" : 0, "f" : 1})
    obs["Gender"].fillna(-1, inplace=True)
    obs["died"] = obs["died"].replace({False : 0, True : 1})
    obs["died"].fillna(-1, inplace=True)
    obs["IGHV"] = obs["IGHV"].replace({"M" : 1, "U" : 0})
    obs["IGHV"].fillna(-1, inplace=True)
    obs.rename(columns={"Gender" : "Sex"}, inplace=True)


    for ome in ["drugs", "methylation", "mrna", "mutations"]:
        modality = pd.read_csv(os.path.join(data_path, f"cll_{ome}.csv"), sep=",", index_col=0, encoding="utf_8").T
        data_dict[ome] = AnnData(X=modality)

    # Replace with gene ID with gene name
    gene_ids = pd.read_csv(os.path.join(data_path, "cll_gene_ids.csv"), index_col=0)
    cols = list(data_dict["mrna"].var_names)

    # Replace each value in cols with the corrsponding value in the gene_ids dataframe
    cols = [gene_ids.loc[gene_ids["GENEID"] == gene, "SYMBOL"].item() for gene in cols]
    data_dict["mrna"].var_names = cols

    # avoid duplicated names with the Mutations view
    data_dict["mutations"].var_names = [f"m_{x}" for x in data_dict["mutations"].var_names]

    # Replace drug names
    # Create mapping from drug_id to name
    drug_names = pd.read_csv(os.path.join(data_path, "drugs.txt"), sep=",", index_col=0)
    mapping = drug_names["name"].to_dict()

    # Replace all substrings in drugs.columns as keys with the corresponding values in the mapping
    cols = []
    for k in data_dict["drugs"].var_names:
        for v in mapping.keys():
            if v in k:
                cols.append(k.replace(v, mapping[v]))
                break

    data_dict["drugs"].var_names = cols

    for view in data_dict.keys():
        data_dict[view] = data_dict[view][data_dict[view].obs_names.argsort()]
        data_dict[view] = data_dict[view][:, data_dict[view].var_names.argsort()]

    mdata = mu.MuData(data_dict)
    mdata.obs = obs.loc[mdata.obs_names]

    return mdata