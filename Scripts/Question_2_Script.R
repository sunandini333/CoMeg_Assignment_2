setwd("C:/SUNANDINI/Projects/Comeg_ass_2")
# ============================================================
# QUESTION 2: MICROBIAL BETA DIVERSITY ANALYSIS
#
# 1. Bray-Curtis dissimilarity
# 2. Principal Coordinates Analysis (PCoA)
# 3. PERMANOVA
# 4. PERMDISP (multivariate dispersion)
# 5. Alternative PERMANOVA model specification
# 6. Identify strongest metadata association
# ============================================================


# ------------------------------------------------------------
# 1. INSTALL AND LOAD PACKAGES
# ------------------------------------------------------------

packages <- c(
  "tidyverse",
  "vegan",
  "permute",
  "ggplot2"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(tidyverse)
library(vegan)
library(permute)
library(ggplot2)


# ------------------------------------------------------------
# 2. CREATE OUTPUT DIRECTORY
# ------------------------------------------------------------

dir.create(
  "Question2_Results",
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 3. LOAD DATASETS
# ------------------------------------------------------------

metadata <- read.csv(
  "metadata_all.xlsx.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

species_profile <- read.csv(
  "species_profile_norm.xlsx.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Rename first species profile column
colnames(species_profile)[1] <- "sample_id"

# Convert sample IDs to character
metadata$sample_id <- as.character(metadata$sample_id)

species_profile$sample_id <- as.character(
  species_profile$sample_id
)

# Check duplicate sample IDs
if (anyDuplicated(metadata$sample_id)) {
  stop("Duplicate sample IDs found in metadata.")
}

if (anyDuplicated(species_profile$sample_id)) {
  stop("Duplicate sample IDs found in species profile.")
}


# ------------------------------------------------------------
# 4. MATCH SAMPLE IDS
# ------------------------------------------------------------

common_samples <- intersect(
  metadata$sample_id,
  species_profile$sample_id
)

cat(
  "Number of common samples:",
  length(common_samples),
  "\n"
)

# Keep common samples in identical order
metadata <- metadata %>%
  filter(sample_id %in% common_samples) %>%
  arrange(match(sample_id, common_samples))

species_profile <- species_profile %>%
  filter(sample_id %in% common_samples) %>%
  arrange(match(sample_id, common_samples))

stopifnot(
  identical(
    metadata$sample_id,
    species_profile$sample_id
  )
)


# ------------------------------------------------------------
# 5. PREPARE SPECIES ABUNDANCE MATRIX
# ------------------------------------------------------------

species_matrix <- species_profile %>%
  select(-sample_id) %>%
  mutate(
    across(
      everything(),
      ~ as.numeric(as.character(.x))
    )
  ) %>%
  as.matrix()

rownames(species_matrix) <- species_profile$sample_id

# Check invalid values
if (
  anyNA(species_matrix) ||
  any(!is.finite(species_matrix))
) {
  stop("Species matrix contains NA, NaN or Inf.")
}

if (any(species_matrix < 0)) {
  stop("Negative abundance values detected.")
}

# Remove species with zero abundance across all samples
species_matrix <- species_matrix[
  ,
  colSums(species_matrix) > 0,
  drop = FALSE
]

cat(
  "Number of species retained:",
  ncol(species_matrix),
  "\n"
)

# Remove samples with zero total abundance
valid_samples <- rownames(species_matrix)[
  rowSums(species_matrix) > 0
]

species_matrix <- species_matrix[
  valid_samples,
  ,
  drop = FALSE
]

metadata <- metadata %>%
  filter(sample_id %in% valid_samples) %>%
  arrange(match(sample_id, valid_samples))


# ------------------------------------------------------------
# 6. PREPARE METADATA
# ------------------------------------------------------------

# Convert metadata variables to factors
metadata$study_condition <- factor(
  metadata$study_condition
)

metadata$age_category <- factor(
  metadata$age_category
)

metadata$gender <- factor(
  metadata$gender
)

# Retain samples with complete metadata
metadata <- metadata %>%
  filter(
    complete.cases(
      study_condition,
      age_category,
      gender
    )
  )

# Keep species matrix in identical order
species_matrix <- species_matrix[
  metadata$sample_id,
  ,
  drop = FALSE
]

metadata <- droplevels(metadata)

cat("\nFINAL DATA DIMENSIONS:\n")
cat("Samples:", nrow(species_matrix), "\n")
cat("Species:", ncol(species_matrix), "\n")

cat("\nStudy condition distribution:\n")
print(table(metadata$study_condition))

cat("\nAge category distribution:\n")
print(table(metadata$age_category))

cat("\nGender distribution:\n")
print(table(metadata$gender))

write.csv(
  metadata,
  "Question2_Results/Final_Analysis_Metadata.csv",
  row.names = FALSE
)


# ============================================================
# QUESTION 2(i): BRAY-CURTIS AND PCoA
# ============================================================


# ------------------------------------------------------------
# 7. CALCULATE BRAY-CURTIS DISSIMILARITY
# ------------------------------------------------------------

# Bray-Curtis dissimilarity
bray_dist <- vegan::vegdist(
  species_matrix,
  method = "bray"
)

cat("\nBRAY-CURTIS DISTANCE MATRIX:\n")

print(bray_dist)

# Convert distance object into a full matrix
bray_matrix <- as.matrix(bray_dist)

# Save distance matrix
write.csv(
  bray_matrix,
  "Question2_Results/Bray_Curtis_Distance_Matrix.csv"
)


# ------------------------------------------------------------
# 8. PRINCIPAL COORDINATES ANALYSIS (PCoA)
# ------------------------------------------------------------

# Classical multidimensional scaling
# Cailliez correction is applied for negative eigenvalues

pcoa_result <- cmdscale(
  bray_dist,
  k = 2,
  eig = TRUE,
  add = TRUE
)

# Extract coordinates
pcoa_coordinates <- as.data.frame(
  pcoa_result$points
)

colnames(pcoa_coordinates) <- c(
  "PCoA1",
  "PCoA2"
)

rownames(pcoa_coordinates) <- rownames(bray_matrix)

# Calculate percentage variation
positive_eigenvalues <- pcoa_result$eig[
  pcoa_result$eig > 0
]

PCoA1_percent <- (
  positive_eigenvalues[1] /
    sum(positive_eigenvalues)
) * 100

PCoA2_percent <- (
  positive_eigenvalues[2] /
    sum(positive_eigenvalues)
) * 100

cat("\nPCoA VARIANCE EXPLAINED:\n")

cat(
  "PCoA1:",
  round(PCoA1_percent, 2),
  "%\n"
)

cat(
  "PCoA2:",
  round(PCoA2_percent, 2),
  "%\n"
)

cat(
  "Total variation represented by first two axes:",
  round(PCoA1_percent + PCoA2_percent, 2),
  "%\n"
)


# ------------------------------------------------------------
# 9. COMBINE PCoA COORDINATES WITH METADATA
# ------------------------------------------------------------

pcoa_coordinates$sample_id <- rownames(
  pcoa_coordinates
)

pcoa_data <- pcoa_coordinates %>%
  left_join(
    metadata,
    by = "sample_id"
  )

write.csv(
  pcoa_data,
  "Question2_Results/PCoA_Coordinates.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 10. CREATE PCoA PLOT
# ------------------------------------------------------------

# Colour = Study condition
# Shape = Age category

pcoa_plot <- ggplot(
  pcoa_data,
  aes(
    x = PCoA1,
    y = PCoA2,
    color = study_condition,
    shape = age_category
  )
) +
  
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  
  labs(
    title = "PCoA of Microbial Community Composition",
    subtitle = "Bray-Curtis dissimilarity",
    x = paste0(
      "PCoA1 (",
      round(PCoA1_percent, 2),
      "%)"
    ),
    y = paste0(
      "PCoA2 (",
      round(PCoA2_percent, 2),
      "%)"
    ),
    color = "Study Condition",
    shape = "Age Category"
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    legend.position = "right"
  )

print(pcoa_plot)

ggsave(
  "Question2_Results/PCoA_BrayCurtis.png",
  pcoa_plot,
  width = 9,
  height = 7,
  dpi = 300
)


# ============================================================
# QUESTION 2(ii): PRIMARY PERMANOVA
# ============================================================


# ------------------------------------------------------------
# 11. GENERATE A COMMON PERMUTATION MATRIX
# ------------------------------------------------------------

# Use the same permutations for all analyses
# to ensure a fair comparison.

set.seed(123)

n_permutations <- 9999

permutation_matrix <- permute::shuffleSet(
  n = nrow(metadata),
  nset = n_permutations,
  control = permute::how()
)

cat(
  "\nNumber of permutations:",
  nrow(permutation_matrix),
  "\n"
)

# Strategy:
# Unrestricted permutation of sample labels.
# This assumes all samples are independent.


# ------------------------------------------------------------
# 12. PRIMARY PERMANOVA: MARGINAL EFFECTS
# ------------------------------------------------------------

# by = "margin" estimates the marginal association
# of each variable after adjusting for other variables.

permanova_primary <- vegan::adonis2(
  bray_dist ~ study_condition +
    age_category +
    gender,
  data = metadata,
  permutations = permutation_matrix,
  by = "margin"
)

cat("\nPRIMARY PERMANOVA RESULTS:\n")

print(permanova_primary)


# ------------------------------------------------------------
# 13. EXTRACT PERMANOVA RESULTS
# ------------------------------------------------------------

primary_results <- as.data.frame(
  permanova_primary
)

primary_results$Variable <- rownames(
  primary_results
)

primary_results <- primary_results %>%
  filter(
    Variable %in% c(
      "study_condition",
      "age_category",
      "gender"
    )
  ) %>%
  select(
    Variable,
    Df,
    SumOfSqs,
    R2,
    F,
    `Pr(>F)`
  ) %>%
  arrange(desc(R2))

cat("\nPRIMARY PERMANOVA SUMMARY:\n")

print(primary_results)

write.csv(
  primary_results,
  "Question2_Results/Primary_PERMANOVA_Marginal.csv",
  row.names = FALSE
)


# ============================================================
# QUESTION 2(iii): MULTIVARIATE DISPERSION
# ============================================================


# ------------------------------------------------------------
# 14. TEST DISPERSION BETWEEN CONTROL AND DISEASED
# ------------------------------------------------------------

# betadisper measures the distance of samples
# from their group centroid (spatial median here).

dispersion_model <- vegan::betadisper(
  bray_dist,
  group = metadata$study_condition,
  type = "median"
)

# Permutation-based test
dispersion_test <- vegan::permutest(
  dispersion_model,
  permutations = permutation_matrix,
  pairwise = TRUE
)

cat("\nPERMDISP RESULTS:\n")

print(dispersion_test)



# ------------------------------------------------------------
# 15. SAVE DISPERSION RESULTS
# ------------------------------------------------------------

dispersion_output <- data.frame(
  sample_id = metadata$sample_id,
  
  study_condition = as.character(
    dispersion_model$group
  ),
  
  Distance_to_Centroid = as.numeric(
    dispersion_model$distances
  )
)

# Check dimensions
cat("Dispersion output dimensions:\n")
print(dim(dispersion_output))

# Check group distribution
print(
  table(dispersion_output$study_condition)
)

# Save results
write.csv(
  dispersion_output,
  "Question2_Results/Dispersion_Distances.csv",
  row.names = FALSE
)

# Save permutation test results
capture.output(
  print(dispersion_test),
  file = "Question2_Results/Dispersion_Test.txt"
)

cat("Dispersion results saved successfully.\n")


# ------------------------------------------------------------
# 16. VISUALIZE MULTIVARIATE DISPERSION
# ------------------------------------------------------------

dispersion_plot <- ggplot(
  dispersion_output,
  aes(
    x = study_condition,
    y = Distance_to_Centroid,
    fill = study_condition
  )
) +
  
  geom_boxplot(
    alpha = 0.65,
    outlier.shape = NA
  ) +
  
  geom_jitter(
    width = 0.15,
    alpha = 0.5,
    size = 1.5
  ) +
  
  labs(
    title = "Multivariate Dispersion Across Study Conditions",
    subtitle = "Distance from group spatial median",
    x = "Study Condition",
    y = "Bray-Curtis Distance to Group Median"
  ) +
  
  theme_classic(base_size = 13)

print(dispersion_plot)

ggsave(
  "Question2_Results/Dispersion_Plot.png",
  dispersion_plot,
  width = 8,
  height = 6,
  dpi = 300
)


# ============================================================
# QUESTION 2(iv): ALTERNATIVE PERMANOVA
# ============================================================


# ------------------------------------------------------------
# 17. ALTERNATIVE MODEL: DIFFERENT TERM ORDER
# ------------------------------------------------------------

# Primary analysis:
# study_condition -> age_category -> gender
#
# Alternative analysis:
# age_category -> gender -> study_condition
#
# Sequential PERMANOVA assigns explained variation
# according to the order of terms.

permanova_alternative <- vegan::adonis2(
  bray_dist ~ age_category +
    gender +
    study_condition,
  data = metadata,
  permutations = permutation_matrix,
  by = "terms"
)

cat("\nALTERNATIVE SEQUENTIAL PERMANOVA:\n")

print(permanova_alternative)


# ------------------------------------------------------------
# 18. ANOTHER SEQUENTIAL MODEL: CONDITION FIRST
# ------------------------------------------------------------

permanova_sequential_primary <- vegan::adonis2(
  bray_dist ~ study_condition +
    age_category +
    gender,
  data = metadata,
  permutations = permutation_matrix,
  by = "terms"
)

cat("\nPRIMARY ORDER SEQUENTIAL PERMANOVA:\n")

print(permanova_sequential_primary)


# ------------------------------------------------------------
# 19. EXTRACT ALTERNATIVE RESULTS
# ------------------------------------------------------------

alternative_results <- as.data.frame(
  permanova_alternative
)

alternative_results$Variable <- rownames(
  alternative_results
)

alternative_results <- alternative_results %>%
  filter(
    Variable %in% c(
      "study_condition",
      "age_category",
      "gender"
    )
  ) %>%
  select(
    Variable,
    Df,
    SumOfSqs,
    R2,
    F,
    `Pr(>F)`
  )

write.csv(
  alternative_results,
  "Question2_Results/Alternative_PERMANOVA.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 20. COMPARE SEQUENTIAL MODEL ORDERS
# ------------------------------------------------------------

sequential_primary <- as.data.frame(
  permanova_sequential_primary
)

sequential_primary$Variable <- rownames(
  sequential_primary
)

sequential_primary <- sequential_primary %>%
  filter(
    Variable %in% c(
      "study_condition",
      "age_category",
      "gender"
    )
  ) %>%
  select(
    Variable,
    Primary_R2 = R2,
    Primary_F = F,
    Primary_p = `Pr(>F)`
  )

sequential_alternative <- as.data.frame(
  permanova_alternative
)

sequential_alternative$Variable <- rownames(
  sequential_alternative
)

sequential_alternative <- sequential_alternative %>%
  filter(
    Variable %in% c(
      "study_condition",
      "age_category",
      "gender"
    )
  ) %>%
  select(
    Variable,
    Alternative_R2 = R2,
    Alternative_F = F,
    Alternative_p = `Pr(>F)`
  )

model_comparison <- full_join(
  sequential_primary,
  sequential_alternative,
  by = "Variable"
)

cat("\nSEQUENTIAL MODEL COMPARISON:\n")

print(model_comparison)

write.csv(
  model_comparison,
  "Question2_Results/Sequential_Model_Comparison.csv",
  row.names = FALSE
)


# ============================================================
# QUESTION 2(v): STRONGEST METADATA ASSOCIATION
# ============================================================


# ------------------------------------------------------------
# 21. IDENTIFY STRONGEST ASSOCIATION
# ------------------------------------------------------------

# Use marginal PERMANOVA R-squared values
# to compare adjusted effect sizes.

strongest_variable <- primary_results %>%
  arrange(desc(R2)) %>%
  slice(1)

cat("\nSTRONGEST METADATA ASSOCIATION:\n")

print(strongest_variable)

cat(
  "\nVariable:",
  strongest_variable$Variable,
  "\n"
)

cat(
  "Marginal R-squared:",
  round(strongest_variable$R2, 4),
  "\n"
)

cat(
  "Pseudo-F:",
  round(strongest_variable$F, 4),
  "\n"
)

cat(
  "Permutation p-value:",
  strongest_variable$`Pr(>F)`,
  "\n"
)


# ------------------------------------------------------------
# 22. SAVE STRONGEST ASSOCIATION
# ------------------------------------------------------------

write.csv(
  strongest_variable,
  "Question2_Results/Strongest_Metadata_Association.csv",
  row.names = FALSE
)


# ============================================================
# 23. SAVE ALL RESULTS IN ONE TEXT FILE
# ============================================================

capture.output(
  {
    cat("====================================\n")
    cat("QUESTION 2: BETA DIVERSITY ANALYSIS\n")
    cat("====================================\n")
    
    cat("\nSAMPLE SIZE:\n")
    print(nrow(metadata))
    
    cat("\nSPECIES COUNT:\n")
    print(ncol(species_matrix))
    
    cat("\nPCoA VARIATION:\n")
    print(PCoA1_percent)
    print(PCoA2_percent)
    
    cat("\nPRIMARY MARGINAL PERMANOVA:\n")
    print(permanova_primary)
    
    cat("\nPRIMARY SEQUENTIAL PERMANOVA:\n")
    print(permanova_sequential_primary)
    
    cat("\nALTERNATIVE SEQUENTIAL PERMANOVA:\n")
    print(permanova_alternative)
    
    cat("\nPERMDISP:\n")
    print(dispersion_test)
    
    cat("\nMODEL COMPARISON:\n")
    print(model_comparison)
    
    cat("\nSTRONGEST ASSOCIATION:\n")
    print(strongest_variable)
    
  },
  file = "Question2_Results/Complete_Statistical_Results.txt"
)


# ------------------------------------------------------------
# 24. COMPLETION MESSAGE
# ------------------------------------------------------------

cat("\n============================================\n")
cat("QUESTION 2 ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("============================================\n")

cat("\nAll output files are saved in:\n")
cat("Question2_Results/\n")

