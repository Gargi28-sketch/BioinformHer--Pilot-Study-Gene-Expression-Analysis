# Pilot Study: 8 genes × 10 patients

set.seed(123)  # for reproducibility

#Generate Dataset 
patient_id <- paste0("P", 1:10)
status <- sample(c("Responder", "Non-Responder"), size = 10, replace = TRUE)

# simulate expression (normal distribution, mean=10, sd=3)
gene_data <- matrix(rnorm(10 * 8, mean = 10, sd = 3), nrow = 10, ncol = 8)
colnames(gene_data) <- paste0("G", 1:8)

expr <- data.frame(patient_id, status, gene_data)

cat("Dataset generated:\n")
print(expr)
cat("\n-------\n")

# 1. Structure & summary 
cat("Structure of expr:\n")
str(expr)
cat("\nSummary of gene columns:\n")
print(summary(expr[, 3:10]))
cat("\n-------\n")

# 2. Per-gene overall means 
gene_means <- colMeans(expr[, 3:10])
cat("Per-gene overall means:\n")
print(gene_means)
cat("Top 2 highest-mean genes:\n")
print(sort(gene_means, decreasing = TRUE)[1:2])
cat("\n---------\n")

# 3. Responder vs Non-Responder means (tapply) 
group_means <- sapply(expr[, 3:10], function(g) tapply(g, expr$status, mean))
cat("Responder vs Non-Responder means:\n")
print(group_means)
cat("\n------------\n")

# 4. Genes higher in Responders 
higher_in_resp <- names(which(group_means["Responder", ] > group_means["Non-Responder", ]))
cat("Genes higher in Responders:\n")
print(higher_in_resp)
cat("\n-----------\n")

# 5. Total expression per patient
expr$total <- rowSums(expr[, 3:10])
cat("Total expression per patient:\n")
print(expr[, c("patient_id", "total")])
cat("\n-----------\n")

# 6. Label HighExpr if total > median 
median_total <- median(expr$total)
expr$HighExpr <- ifelse(expr$total > median_total, "HighExpr", "LowExpr")
cat("HighExpr / LowExpr count per group:\n")
print(table(expr$status, expr$HighExpr))
cat("\n---------\n")

# 7. Patient with highest total 
cat("Patient with highest total:\n")
print(expr[which.max(expr$total), ])
cat("\n--------\n")

# 8. Per-gene variance 
gene_vars <- apply(expr[, 3:10], 2, var)
cat("Per-gene variance:\n")
print(gene_vars)
cat("Top 2 most variable genes:\n")
top_var_genes <- names(sort(gene_vars, decreasing = TRUE)[1:2])
print(top_var_genes)
cat("\n-----------\n")

# 9. Highest single patient × gene 
which_max <- which(expr[, 3:10] == max(expr[, 3:10]), arr.ind = TRUE)
cat("Patient × gene combination with highest single value:\n")
print(expr[which_max[1], 1:2])   # patient & status
cat("Gene:", colnames(expr)[2 + which_max[2]], "\n")
cat("Value:", max(expr[, 3:10]), "\n")
cat("\n----------\n")

# 10. Responder-only means vs overall 
resp_means <- colMeans(expr[expr$status == "Responder", 3:10])
cat("Responder-only per-gene means:\n")
print(resp_means)
cat("Comparison to overall means:\n")
print(data.frame(overall = gene_means, responder_only = resp_means))
cat("\n-----------\n")

# 11. Rename G4 to G4A 
colnames(expr)[colnames(expr) == "G4"] <- "G4A"
cat("G4 renamed to G4A, columns now:\n")
print(colnames(expr))
cat("\n---------\n")

# 12. Create gene_mat (10×8) 
gene_mat <- as.matrix(expr[, 3:10])
rownames(gene_mat) <- expr$patient_id
cat("gene_mat created (first 5 rows):\n")
print(gene_mat[1:5, ])
cat("\n----------\n")

# 13. Row & column means 
row_means <- rowMeans(gene_mat)
col_means <- colMeans(gene_mat)
cat("Do column means match earlier gene means?\n")
print(all.equal(unname(col_means), unname(gene_means)))
cat("Row means (first 5 patients):\n")
print(row_means[1:5])
cat("Column means:\n")
print(col_means)
cat("\n--------\n")

# 14. Create list "study"
study <- list(expr = expr, genes = gene_mat)

# 15. Extract G1 expression for Non-Responders 
cat("G1 expression for Non-Responders:\n")
print(study$expr$G1[study$expr$status == "Non-Responder"])
cat("\n-------\n")

# 16. Threshold gene_mat (<5 → 0) 
changed <- sum(gene_mat < 5)
gene_mat[gene_mat < 5] <- 0
cat("Number of entries changed to zero (<5):", changed, "\n")
cat("\n---------\n")

# 17. Responder vs Non-Responder means post-threshold 
group_means_thresh <- sapply(as.data.frame(gene_mat), function(g) tapply(g, expr$status, mean))
cat("Responder vs Non-Responder means post-threshold:\n")
print(group_means_thresh)
higher_in_resp_thresh <- names(
  which(group_means_thresh["Responder", ] >
          group_means_thresh["Non-Responder", ])
)

cat("Genes higher in Responders after thresholding:\n")
print(higher_in_resp_thresh)

cat("Before thresholding:\n")
print(higher_in_resp)

cat("After thresholding:\n")
print(higher_in_resp_thresh)
cat("\n-------\n")

# 18. Barplot per-patient totals
png("Per-Patient_Total_Expression.png", width = 800, height = 600)

barplot(expr$total,
        names.arg = expr$patient_id,
        col = "Red",
        main = "Per-Patient Total Expression",
        xlab = "Patient",
        ylab = "Total Expression")

dev.off()

# 19. Remove lowest-total patient 
lowest_idx <- which.min(expr$total)
cat("Lowest-total patient:\n")
print(expr[lowest_idx, c("patient_id", "status", "total")])
cat("\nMean total per group before removal:\n")
print(tapply(expr$total, expr$status, mean))
expr_removed <- expr[-lowest_idx, ]

cat("\nMean total per group after removal:\n")
print(tapply(expr_removed$total, expr_removed$status, mean))


# 20. 2-paragraph interpretation 
cat("\nINTERPRETATION:\n\n")

cat("In this pilot study of 10 patients, several genes showed higher average ",
    "expression in Responders compared with Non-Responders. Before thresholding, ",
    "G2, G3, G4, and G7 were higher in Responders. After setting expression values ",
    "below 5 to zero, G3, G4A, and G7 remained higher in Responders, while G2 ",
    "was no longer higher. This suggests that G3, G4A, and G7 showed a more ",
    "consistent pattern of higher expression in Responders in this simulated dataset. ",
    "However, these results describe observed patterns and do not establish that ",
    "these genes cause or predict treatment response.\n\n", sep="")

cat("The patient-level analysis showed that Responders had a higher mean total ",
    "expression than Non-Responders. The lowest-total patient was P3, a Responder, ",
    "with a total expression of 67.43. Before removing P3, the mean total expression ",
    "was 81.98 for Responders and 77.80 for Non-Responders. After removing P3, the ",
    "Responder mean increased to 84.89, while the Non-Responder mean remained 77.80. ",
    "Therefore, removing the lowest-total patient did not change which group had the ",
    "higher mean. Because this is a small simulated pilot dataset with only 10 patients, ",
    "the findings should be considered preliminary and would need to be tested in a ",
    "larger real dataset.\n", sep="")
