library(Matrix)
library(data.table)

load_rnaseq_sample_selected_tissues <- function(
  X_path = "results/rnaseq_sample_selected_tissues.mm",
  genes_rows_path = "results/rnaseq_sample_selected_tissues_genes_rows.csv",
  samples_columns_path = "results/rnaseq_sample_selected_tissues_samples_columns.csv",
  samples_metadata_path = "results/rnaseq_sample_selected_tissues_metadata.csv"
) {
  # Load the sparse matrix
  X <- readMM(X_path)
  
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

# rnaseq_sample_selected_tissues <- load_rnaseq_sample_selected_tissues(
#     X_path = "results/rnaseq_sample_selected_tissues.mm",
#   genes_rows_path = "results/rnaseq_sample_selected_tissues_genes_rows.csv",
#   samples_columns_path = "results/rnaseq_sample_selected_tissues_samples_columns.csv",
#   samples_metadata_path = "results/rnaseq_sample_selected_tissues_metadata.csv")

# metadata <- data.frame(rnaseq_sample_selected_tissues$metadata)
# rownames(metadata) <- metadata$SAMPID
# expr_set <- Biobase::ExpressionSet(as.matrix(rnaseq_sample_selected_tissues$X),phenoData=Biobase::AnnotatedDataFrame(metadata))
# # yarn::checkMisAnnotation(expr_set,"GENDER",controlGenes="Y",legendPosition="topleft")
# yarn::checkTissuesToMerge(expr_set,"SMTS","SMTSD")
# expr_set_filtered <- yarn::filterLowGenes(expr_set,"SMTSD")
# expr_set_filtered_norm <- yarn::normalizeTissueAware(expr_set_filtered,"SMTSD")

# # Convert to sparse matrix format
# rnaseq_sample_selected_tissues_yarn_normalized_sparse <- as( expr_set_filtered_norm@assayData$normalizedMatrix, "sparseMatrix")
# writeMM(rnaseq_sample_selected_tissues_yarn_normalized_sparse, "results/rnaseq_sample_selected_tissues_yarn_normalized.mm")
# write.csv(data.frame(genes = rownames(rnaseq_sample_selected_tissues_yarn_normalized_sparse)), "rnaseq_sample_selected_tissues_yarn_normalized_genes_rows.csv", quote=FALSE)
# write.csv(data.frame(samples =  colnames(rnaseq_sample_selected_tissues_yarn_normalized_sparse)), "rnaseq_sample_selected_tissues_yarn_normalized_samples_columns.csv", quote=FALSE)
