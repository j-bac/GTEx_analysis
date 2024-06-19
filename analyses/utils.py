import pandas as pd
import pathlib
import scipy
import scanpy as sc
import numpy as np


def get_protein_coding_genes(
    file_path="../data/gencode.v26.annotation.gtf.gz",
    save_path="../results/protein_coding_genes.txt",
    overwrite=False,
):
    if pathlib.Path(save_path).is_file() and not overwrite:
        print("File already exists and overwrite is False")
    else:
        print("Loading gencode file with gene metadata...")
        gtf = pd.read_csv(file_path, sep="\t", comment="#", header=None)

        protein_coding_genes = gtf[
            (gtf[2] == "gene") & (gtf[8].str.contains('gene_type "protein_coding"'))
        ].copy()
        protein_coding_genes[["gene_id", "gene_name"]] = protein_coding_genes[
            8
        ].str.extract('gene_id "([^"]+)"; .* gene_name "([^"]+)"')

        print(
            len(protein_coding_genes),
            f"protein coding genes found. Saving to {save_path}.",
        )
        protein_coding_genes.to_csv(save_path)


def load_rnaseq_sample_selected_tissues_anndata(
    X_path="../results/rnaseq_sample_selected_tissues.parquet",
    genes_rows_path="../results/rnaseq_sample_selected_tissues_genes_rows.csv",
    samples_columns_path="../results/rnaseq_sample_selected_tissues_samples_columns.csv",
    samples_metadata_path="../results/rnaseq_sample_selected_tissues_metadata.csv",
):
    X = pd.read_parquet(X_path)
    rows = pd.read_csv(genes_rows_path)
    cols = pd.read_csv(samples_columns_path, index_col=0)["0"].values
    cols_metadata = pd.read_csv(samples_metadata_path, index_col=0)

    assert all(cols == cols_metadata.index)

    adata = sc.AnnData(X)
    adata.var_names = cols
    adata.var = cols_metadata
    adata.obs_names = rows["Name"].values
    if "Description" in rows.columns:
        adata.obs["Description"] = rows["Description"].values
    adata = adata.T.copy()
    return adata


def is_outlier(adata, metric: str, nmads: int, positive_only=False):
    M = adata.obs[metric]
    threshold = nmads * scipy.stats.median_abs_deviation(M)
    if positive_only:
        outlier = np.median(M) + threshold < M
    else:
        outlier = (M < np.median(M) - threshold) | (np.median(M) + threshold < M)
    return outlier, (np.median(M), threshold)
