# README: Pilot Study Gene Expression Analysis

## Project Overview

This project performs a simulated analysis of gene expression data for a pilot study measuring 8 genes across 10 patients. Each patient is labeled as either a **Responder** or **Non-Responder**. The analysis includes computing gene-level and patient-level statistics, visualizing results, and providing interpretations.

## Files

* **Script.R**: Complete R script performing all steps (1–20) of the analysis, including dataset generation, statistical summaries, variance analysis, thresholding, barplots, and interpretation.

## Dataset Generation

Since no dataset was provided, the script generates a synthetic dataset with the following parameters:

* 10 patients (`P1`–`P10`)
* 8 genes (`G1`–`G8`)
* Random Responder / Non-Responder assignment
* Gene expression values simulated from a normal distribution with mean = 10 and sd = 3
* Seed is set for reproducibility using `set.seed(123)`

## Analysis Steps

1. Print dataset structure and gene column summary.
2. Compute per-gene overall means and identify the top 2 highest-mean genes.
3. Compute Responder vs Non-Responder mean per gene using `tapply()`.
4. Identify genes with higher expression in Responders.
5. Compute total expression per patient.
6. Label patients as HighExpr or LowExpr based on median total expression and count per group.
7. Identify the patient with the highest total expression.
8. Compute per-gene variance and identify the top 2 most variable genes.
9. Identify the patient × gene combination with the highest single expression value.
10. Compute per-gene means for Responders and compare them with overall means.
11. Rename gene G4 to G4A.
12. Create a gene matrix (`gene_mat`) with rownames as patient IDs.
13. Compute row and column means and verify consistency with previous results.
14. Create a list `study` containing the expression data frame and gene matrix.
15. Extract G1 expression for Non-Responders from `study`.
16. Threshold the gene matrix by setting values < 5 to zero and count the changed entries.
17. Recompute Responder vs Non-Responder means after thresholding and compare the conclusions.
18. Generate a barplot of per-patient total expression.
19. Remove the patient with the lowest total expression and observe group-level changes.
20. Print a 2-paragraph interpretation of the results automatically to the console.

## Output

* Dataset preview
* Structure and summaries
* Gene means and variances
* Responder vs Non-Responder comparisons
* Total expression per patient and HighExpr counts
* Top variable genes
* Genes with higher expression in Responders
* Per-patient total expression barplot
* Interpretation of the simulated results

## Notes

* All computations use **base R** functions.
* The dataset is simulated because no dataset was provided.
* The results represent patterns in a small simulated dataset and should not be interpreted as evidence of biological causation or clinical prediction.

**Full R script:** [R Script](https://github.com/Gargi28-sketch/BioinformHer--Pilot-Study-Gene-Expression-Analysis/blob/2d346e4543da1fae510cdb544c42fc1f62142241/Script.R)
**Per-patient total expression barplot:** [https://github.com/Gargi28-sketch/BioinformHer--Pilot-Study-Gene-Expression-Analysis/blob/2d346e4543da1fae510cdb544c42fc1f62142241/Per-Patient_Total_Expression.png]

---

**Author:** Gargi Durbude

**Contact:** gauridilip2001@gmail.com

**Purpose:** Pilot study gene expression analysis and summary for coursework submission.
