library(Matrix)
library(data.table)
library(arrow)

load_rnaseq_sample_selected_tissues <- function(
  X_path = "results/0_rnaseq_sample_selected_tissues_QC.parquet",
  samples_metadata_path = "results/0_rnaseq_sample_selected_tissues_QC_metadata.csv"
) {
  # Load the sparse matrix
  X <- read_parquet(X_path)

  # assign rownames
  rownames <- X$Name
  X$Name <- NULL
  X <- as.matrix(X)
  rownames(X) <- rownames 

  # Load the sample metadata
  cols_metadata <- fread(samples_metadata_path)

  # Ensure the column order matches the metadata index
  stopifnot(all(colnames(X) == cols_metadata[[1]]))

  # Create a sparse data frame-like object
  out <- list(X = X,  metadata = cols_metadata)
  
  return(out)
}

