here::set_here()
library(yarn)
library(quantro)
library(qsmooth)
library(dplyr)
source("analyses/preprocessing.R")

rnaseq_sample_selected_tissues <- load_rnaseq_sample_selected_tissues(
  X_path = "results/rnaseq_sample_selected_tissues.mm",
  genes_rows_path = "results/rnaseq_sample_selected_tissues_genes_rows.csv",
  samples_columns_path = "results/rnaseq_sample_selected_tissues_samples_columns.csv",
  samples_metadata_path = "results/rnaseq_sample_selected_tissues_metadata.csv"
)

# convert to ExpressionSet
metadata <- data.frame(rnaseq_sample_selected_tissues$metadata)
rownames(metadata) <- metadata$SAMPID
expr_set <- Biobase::ExpressionSet(
    as.matrix(rnaseq_sample_selected_tissues$X),
    phenoData=Biobase::AnnotatedDataFrame(metadata)
)

# process with yarn qsmooth
# yarn::checkMisAnnotation(expr_set,"GENDER",controlGenes="Y",legendPosition="topleft")
# yarn::checkTissuesToMerge(expr_set,"SMTS","SMTSD")
expr_set_filtered <- yarn::filterLowGenes(expr_set,"SMTSD")
expr_set_filtered_norm <- yarn::normalizeTissueAware(expr_set_filtered,"SMTSD")
# qs_norm_e1 <- qsmooth(object = expr_set_filtered@assayData$exprs, group_factor = expr_set$SMTSD)


# save to MM (just to keep single format - matrix is now quite dense)
rnaseq_sample_selected_tissues_yarn_normalized_sparse <- as( expr_set_filtered_norm@assayData$normalizedMatrix, "sparseMatrix")
writeMM(rnaseq_sample_selected_tissues_yarn_normalized_sparse, "results/rnaseq_sample_selected_tissues_yarn_normalized.mm")
write.csv(data.frame(Name = rownames(expr_set_filtered)), "results/rnaseq_sample_selected_tissues_yarn_normalized_genes_rows.csv", quote=FALSE)
write.csv(data.frame('0' =  colnames(expr_set_filtered)), "results/rnaseq_sample_selected_tissues_yarn_normalized_samples_columns.csv", quote=FALSE)
