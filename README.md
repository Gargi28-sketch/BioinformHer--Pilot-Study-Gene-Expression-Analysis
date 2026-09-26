# README: Pilot Study Gene Expression Analysis

## Project Overview

This project performs a simulated analysis of gene expression data for a pilot study measuring 8 genes across 10 patients. Each patient is labeled as either a **Responder** or **Non-Responder**. The analysis includes computing gene-level and patient-level statistics, visualizing results, and providing interpretations.

## Files

* **analysis.R**: Complete R script performing all steps (1–20) of the analysis, including dataset generation, statistical summaries, variance analysis, thresholding, barplots, and automated interpretation.

## Dataset Generation

Since no dataset was provided, the script generates a synthetic dataset with the following parameters:

* 10 patients (`P1`–`P10`)
* 8 genes (`G1`–`G8`)
* Random Responder / Non-Responder assignment
* Gene expression values simulated from a normal distribution with mean = 10, sd = 3
* Seed is set for reproducibility (`set.seed(123)`)

## Analysis Steps

1. Print dataset structure and gene column summary.
2. Compute per-gene overall means and identify top 2 highest-mean genes.
3. Compute Responder vs Non-Responder mean per gene using `tapply()`.
4. Identify genes with higher expression in Responders.
5. Compute total expression per patient.
6. Label patients as HighExpr or LowExpr based on median total expression and count per group.
7. Identify the patient with the highest total.
8. Compute per-gene variance and top 2 most variable genes.
9. Identify patient × gene combination with the highest single value.
10. Compute per-gene means for Responders and compare to overall means.
11. Rename gene G4 to G4A.
12. Create a gene matrix (`gene_mat`) with rownames as patient IDs.
13. Compute row and column means and verify consistency with previous results.
14. Create a list `study` containing the dataframe and gene matrix.
15. Extract G1 expression for Non-Responders from `study`.
16. Threshold gene matrix: set values < 5 to zero and count changed entries.
17. Recompute Responder vs Non-Responder means post-threshold.
18. Generate a barplot of per-patient total expression.
19. Remove the patient with lowest total and observe group-level changes.
20. Print a 2-paragraph interpretation of results automatically to the console.

## Output

* Dataset preview
* Structure and summaries
* Gene means and variances
* Responder vs Non-Responder comparisons
* Total expression per patient and HighExpr counts
* Top variable genes and genes higher in Responders
* Interpretation paragraphs

  
Full R script is available here- https://github.com/Gargi28-sketch/BioinformHer-R-Class-Project/blob/e76ed00d48b3150ad1c4bbf4fb5fe36a9f6b0c32/Project%20script.R

A barplot of total expression per patient available here-https://github.com/Gargi28-sketch/BioinformHer-R-Class-Project/blob/d59710c21834838ba8cc76732ff482b3a3152e4a/Bar%20Plot.png


## Notes
* All computations use **base R** functions.
---

**Author:** Gargi Durbude

**Contact** gauridilip2001@gmail.com

**Purpose:** Pilot study gene expression analysis and summary for coursework submission.
