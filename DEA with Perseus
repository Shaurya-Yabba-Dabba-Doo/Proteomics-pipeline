library(tidyverse)
library(limma)

setwd("path/to/your/directory")

data <- read.table("maxquant_protein.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

head(data)

protein_data <- data %>%
  select(Protein.IDs, starts_with("Intensity"))

protein_matrix <- as.matrix(protein_data[,-1])
rownames(protein_matrix) <- protein_data$Protein.IDs

protein_matrix_log2 <- log2(protein_matrix + 1)

normalized_data <- normalizeBetweenArrays(protein_matrix_log2, method = "quantile")

normalized_df <- as.data.frame(normalized_data)
normalized_df$Protein.IDs <- rownames(normalized_df)

write.table(normalized_df, "normalized_protein_data.txt", sep = "\t", row.names = FALSE, quote = FALSE)

condition <- factor(c(rep("Control", n_control_samples), rep("Treatment", n_treatment_samples)))
design <- model.matrix(~ condition)

fit <- lmFit(normalized_data, design)

fit <- eBayes(fit)

results <- topTable(fit, adjust = "fdr", sort.by = "P", number = Inf)

write.table(results, "differential_expression_results.txt", sep = "\t", row.names = TRUE, quote = FALSE)

print(head(results))
