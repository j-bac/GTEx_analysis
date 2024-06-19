library(Matrix)
library(data.table)

load_rnaseq_sample_selected_tissues <- function(
  X_path = "results/rnaseq_sample_selected_tissues.mm",
  genes_rows_path = "results/rnaseq_sample_selected_tissues_genes_rows.csv",
  samples_columns_path = "results/rnaseq_sample_selected_tissues_samples_columns.csv",
  samples_metadata_path = "results/rnaseq_sample_selected_tissues_metadata.csv"
) {
  # Load the sparse matrix
  X <- read_parquet(X_path)
  
  # Load the rows (genes)
  rows <- fread(genes_rows_path)
  rownames(X) <- rows[['Name']]  # Assuming the first column contains gene names

  # Load the columns (samples)
  cols <- fread(samples_columns_path)[['V2']]
  colnames(X) <- cols[2:length(cols)]

  # Load the sample metadata
  cols_metadata <- fread(samples_metadata_path)

  # Ensure the column order matches the metadata index
  stopifnot(all(colnames(X) == cols_metadata[[1]]))

  # Create a sparse data frame-like object
  sparse_df <- list(X = X,  metadata = cols_metadata)
  
  return(sparse_df)
}