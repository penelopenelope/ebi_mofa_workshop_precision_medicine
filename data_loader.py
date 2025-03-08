import pandas as pd
from anndata import AnnData
import os
import muon as mu
import numpy as np
import anndata as ad
import mudata as md


def load_cll(data_path="data/"):
    """Load and return CLL dataset."""
    metadata_df = pd.read_csv("data/cll_metadata.csv", index_col=0)

    # assuming you called the DataFrame `metadata_df`
    metadata_df.rename(columns={"Gender" : "Sex"}, inplace=True)
    metadata_df.replace({"Sex" : {"m" : 0, "f" : 1}}, inplace=True)
    metadata_df.replace({"died" : {False : 0, True : 1}}, inplace=True)
    metadata_df.replace({"IGHV" : {"U" : 0, "M" : 1}}, inplace=True)
    metadata_df.fillna(-1, inplace=True)

    mrna_df = pd.read_csv("data/cll_mrna.csv", index_col=0).T
    mrna_adata = ad.AnnData(mrna_df)

    # assuming you called the AnnData object `adata_mrna`
    gene_ids = pd.read_csv(os.path.join("data/cll_gene_ids.csv"), index_col=0)
    cols = list(mrna_adata.var_names)
    cols = [gene_ids.loc[gene_ids["GENEID"] == gene, "SYMBOL"].item() for gene in cols]
    mrna_adata.var_names = cols

    mutations_df = pd.read_csv("data/cll_mutations.csv", index_col=0).T
    mutations_adata = ad.AnnData(mutations_df)

    # assuming you called the AnnData object `adata_mutations`
    mutations_adata.var_names = [f"m_{x}" for x in mutations_adata.var_names]

    methylation_df = pd.read_csv("data/cll_methylation.csv", index_col=0).T
    methylation_adata = ad.AnnData(methylation_df)

    drugs_df = pd.read_csv("data/cll_drugs.csv", index_col=0).T
    drugs_adata = ad.AnnData(drugs_df)

    # assuming you called the AnnData object `adata_drugs`
    drug_names = pd.read_csv("data/drugs.txt", sep=",", index_col=0)
    mapping = drug_names["name"].to_dict()
    cols = []
    for k in drugs_adata.var_names:
        for v in mapping.keys():
            if v in k:
                cols.append(k.replace(v, mapping[v]))
                break

    drugs_adata.var_names = cols

    modality_dict = {
        "mrna" : mrna_adata,
        "mutations" : mutations_adata,
        "methylation" : methylation_adata,
        "drugs" : drugs_adata
    }

    mdata = md.MuData(modality_dict)
    mdata.obs = metadata_df.loc[mdata.obs_names]

    return mdata