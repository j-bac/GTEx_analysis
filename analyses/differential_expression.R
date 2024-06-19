metadata <- data.frame(
  SampleID = colnames(expr_matrix),
  Condition = rep(c("Control", "Treatment"), each = 5),
  PatientID = rep(c("Patient1", "Patient2", "Patient3", "Patient4", "Patient5"), 2)
)
metadata$SMTSD_grouped <- factor(metadata$SMTSD_grouped)
metadata$donor_ids <- factor(metadata$donor_ids)

# Create the design matrix
design <- model.matrix(~ SMTSD_grouped + donor_ids, data = metadata)
colnames(design) <- make.names(colnames(design))
print(design)

fit <- lmFit(expr_matrix, design)
fit <- eBayes(fit)
results <- topTable(fit, coef = "ConditionTreatment", adjust.method = "BH")
head(results)