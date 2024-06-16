metadata <- data.frame(
  SampleID = colnames(expr_matrix),
  Condition = rep(c("Control", "Treatment"), each = 5),
  PatientID = rep(c("Patient1", "Patient2", "Patient3", "Patient4", "Patient5"), 2)
)
metadata$Condition <- factor(metadata$Condition)
metadata$PatientID <- factor(metadata$PatientID)

# Create the design matrix
design <- model.matrix(~ Condition + PatientID, data = metadata)
colnames(design) <- make.names(colnames(design))
print(design)

# Fit the linear model
fit <- lmFit(expr_matrix, design)

# Apply empirical Bayes moderation
fit <- eBayes(fit)

# Extract the top table with differential expression results for the condition
results <- topTable(fit, coef = "ConditionTreatment", adjust.method = "BH")
head(results)