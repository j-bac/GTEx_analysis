library(Matrix)
library(data.table)

load_rnaseq_sample_selected_tissues <- function(
  X_path = "../results/rnaseq_sample_selected_tissues.mm",
  genes_rows_path = "../results/rnaseq_sample_selected_tissues_genes_rows.csv",
  samples_columns_path = "../results/rnaseq_sample_selected_tissues_samples_columns.csv",
  samples_metadata_path = "../results/rnaseq_sample_selected_tissues_metadata.csv"
) {
  # Load the sparse matrix
  X <- readMM(X_path)
  
  # Load the rows (genes)
  rows <- fread(genes_rows_path)
  rownames(X) <- rows[[2]]  # Assuming the first column contains gene names
  
  # Load the columns (samples)
  cols <- fread(samples_columns_path, select = 1)[[1]]
  colnames(X) <- cols
  
  # Load the sample metadata
  cols_metadata <- fread(samples_metadata_path)
  
  # Ensure the column order matches the metadata index
  stopifnot(all(cols == cols_metadata[[1]]))
  
  # Create a sparse data frame-like object
  sparse_df <- list(matrix = X, rows = rows, cols = cols, metadata = cols_metadata)
  
  return(sparse_df)
}

# Example usage
sparse_df <- load_rnaseq_sample_selected_tissues()
str(sparse_df)
