here::set_here()
library(yarn)
library(quantro)
library(qsmooth)
library(dplyr)
library(arrow)
source("analyses/utils.R")

# load data
rnaseq_sample_selected_tissues <- load_rnaseq_sample_selected_tissues(
  X_path = "results/0_rnaseq_sample_selected_tissues_QC.parquet",
  samples_metadata_path = "results/0_rnaseq_sample_selected_tissues_QC_metadata.csv"
)

# convert to ExpressionSet
metadata <- data.frame(rnaseq_sample_selected_tissues$metadata)
rownames(metadata) <- metadata$SAMPID
expr_set <- Biobase::ExpressionSet(
    as.matrix(rnaseq_sample_selected_tissues$X),
    phenoData=Biobase::AnnotatedDataFrame(metadata)
)
rm(rnaseq_sample_selected_tissues)

# filer genes, process with qsmooth and quantile normalization
expr_set_filtered <- yarn::filterLowGenes(expr_set,"SMTSD_grouped")
rm(expr_set)
expr_set_filtered_qsmooth <- yarn::normalizeTissueAware(expr_set_filtered,"SMTSD_grouped", normalizationMethod='qsmooth')
# expr_set_filtered_quantile <- yarn::normalizeTissueAware(expr_set_filtered,"SMTSD_grouped", normalizationMethod='quantile')

# get qsmooth weights
# qs <- qsmooth(expr_set_filtered@assayData$exprs,"SMTSD_grouped")
# expr_set_filtered_qsmooth_data <- qsmoothData(qs) # extract smoothed quantile normalized data
# expr_set_filtered_qsmooth_weights qsmoothWeights(qs) # extract smoothed quantile normalized weights

# check assumptions with quantro
# qtest <- quantro(object = expr_set_filtered, groupFactor = expr_set_filtered$SMTSD_grouped)
# qtestPerm <- quantro(object = expr_set_filtered, groupFactor = expr_set_filtered$SMTSD_grouped, B = 100)
# quantroPlot(qtestPerm)

# summary(qtest)
# anova(qtest)
# quantroStat(qtest)

# summary(qtestPerm)
# anova(qtestPerm)
# quantroStat(qtestPerm)

# Write data
write_parquet(as.data.frame(expr_set_filtered_qsmooth@assayData$normalizedMatrix), "results/1_rnaseq_sample_selected_tissues_qsmooth.parquet")
# Write gene and sample names
write.csv(data.frame(Name = rownames(expr_set_filtered)), "results/1_rnaseq_sample_selected_tissues_qsmooth_genes.csv")
