setwd("C:/SUNANDINI/Projects/Comeg_ass_2")
# ============================================================
# QUESTION 3: SPECIES PREVALENCE ANALYSIS
#
# i.   Calculate species prevalence
# ii.  Fisher's exact test + BH correction
# iii. Identify top 3 positive and negative differences
# iv.  Visualize prevalence differences
# v.   Compare prevalence and mean abundance
# ============================================================


# ------------------------------------------------------------
# 1. INSTALL AND LOAD PACKAGES
# ------------------------------------------------------------

packages <- c(
  "tidyverse",
  "ggplot2",
  "scales"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(tidyverse)
library(ggplot2)
library(scales)


# ------------------------------------------------------------
# 2. CREATE OUTPUT DIRECTORY
# ------------------------------------------------------------

dir.create(
  "Question3_Results",
  showWarnings = FALSE,
  recursive = TRUE
)


# ------------------------------------------------------------
# 3. LOAD DATA
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

# Rename first species-profile column
colnames(species_profile)[1] <- "sample_id"

# Convert sample IDs to character
metadata$sample_id <- as.character(metadata$sample_id)
species_profile$sample_id <- as.character(species_profile$sample_id)

# Check required metadata columns
required_metadata <- c(
  "sample_id",
  "study_condition"
)

missing_columns <- setdiff(
  required_metadata,
  colnames(metadata)
)

if (length(missing_columns) > 0) {
  stop(
    paste(
      "Missing metadata columns:",
      paste(missing_columns, collapse = ", ")
    )
  )
}

# Check duplicate sample IDs
if (anyDuplicated(metadata$sample_id)) {
  stop("Duplicate sample IDs in metadata.")
}

if (anyDuplicated(species_profile$sample_id)) {
  stop("Duplicate sample IDs in species profile.")
}


# ------------------------------------------------------------
# 4. MATCH SAMPLES
# ------------------------------------------------------------

common_samples <- metadata$sample_id[
  metadata$sample_id %in% species_profile$sample_id
]

metadata <- metadata %>%
  filter(
    sample_id %in% common_samples,
    study_condition %in% c("Control", "Diseased")
  ) %>%
  arrange(match(sample_id, common_samples))

species_profile <- species_profile %>%
  filter(sample_id %in% metadata$sample_id) %>%
  arrange(match(sample_id, metadata$sample_id))

# Verify sample alignment
stopifnot(
  identical(
    metadata$sample_id,
    species_profile$sample_id
  )
)

metadata$study_condition <- factor(
  metadata$study_condition,
  levels = c("Control", "Diseased")
)

cat("\nSAMPLE COUNTS:\n")
print(table(metadata$study_condition))

cat("\nNumber of matched samples:", nrow(metadata), "\n")


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

# Check data validity
if (
  anyNA(species_matrix) ||
  any(!is.finite(species_matrix))
) {
  stop("Species abundance matrix contains NA, NaN, or Inf values.")
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

cat("\nNumber of species retained after removing all-zero species:",
    ncol(species_matrix), "\n")


# ============================================================
# QUESTION 3(i): SPECIES PREVALENCE
# ============================================================


# ------------------------------------------------------------
# 6. CONVERT ABUNDANCE INTO PRESENCE-ABSENCE
# ------------------------------------------------------------

# Presence = 1
# Absence  = 0

presence_absence <- ifelse(
  species_matrix > 0,
  1,
  0
)

rownames(presence_absence) <- rownames(species_matrix)
colnames(presence_absence) <- colnames(species_matrix)

cat("\nPRESENCE-ABSENCE MATRIX DIMENSIONS:\n")
print(dim(presence_absence))


# ------------------------------------------------------------
# 7. DEFINE SAMPLE GROUPS
# ------------------------------------------------------------

control_samples <- metadata$sample_id[
  metadata$study_condition == "Control"
]

diseased_samples <- metadata$sample_id[
  metadata$study_condition == "Diseased"
]

n_control <- length(control_samples)
n_diseased <- length(diseased_samples)

cat("\nControl sample size:", n_control, "\n")
cat("Diseased sample size:", n_diseased, "\n")


# ------------------------------------------------------------
# 8. CALCULATE SPECIES PREVALENCE
# ------------------------------------------------------------

# Number of samples where each species is present

control_count <- colSums(
  presence_absence[control_samples, , drop = FALSE]
)

diseased_count <- colSums(
  presence_absence[diseased_samples, , drop = FALSE]
)

# Prevalence proportions

control_prevalence <- control_count / n_control
diseased_prevalence <- diseased_count / n_diseased

# Prevalence percentages

control_prevalence_percent <- control_prevalence * 100
diseased_prevalence_percent <- diseased_prevalence * 100

# Difference in percentage points
# Positive = higher prevalence in Diseased
# Negative = higher prevalence in Control

prevalence_difference <- (
  diseased_prevalence_percent -
    control_prevalence_percent
)


# ------------------------------------------------------------
# 9. CREATE COMPLETE PREVALENCE TABLE
# ------------------------------------------------------------

prevalence_table <- data.frame(
  
  Species = colnames(species_matrix),
  
  Control_Present = as.numeric(control_count),
  
  Control_Absent = n_control - as.numeric(control_count),
  
  Control_Prevalence = as.numeric(control_prevalence),
  
  Control_Prevalence_Percent =
    as.numeric(control_prevalence_percent),
  
  Diseased_Present = as.numeric(diseased_count),
  
  Diseased_Absent = n_diseased - as.numeric(diseased_count),
  
  Diseased_Prevalence = as.numeric(diseased_prevalence),
  
  Diseased_Prevalence_Percent =
    as.numeric(diseased_prevalence_percent),
  
  Prevalence_Difference_PP =
    as.numeric(prevalence_difference),
  
  Total_Present = as.numeric(
    control_count + diseased_count
  )
)

write.csv(
  prevalence_table,
  "Question3_Results/Complete_Species_Prevalence.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 10. FILTER RARE SPECIES
# ------------------------------------------------------------

# Retain species present in at least 5 samples overall

minimum_prevalence_count <- 5

filtered_prevalence <- prevalence_table %>%
  filter(
    Total_Present >= minimum_prevalence_count
  )

cat("\nNumber of species before filtering:",
    nrow(prevalence_table), "\n")

cat("Number of species retained:",
    nrow(filtered_prevalence), "\n")

cat("Number of species removed:",
    nrow(prevalence_table) - nrow(filtered_prevalence), "\n")

write.csv(
  filtered_prevalence,
  "Question3_Results/Filtered_Species_Prevalence.csv",
  row.names = FALSE
)


# ============================================================
# QUESTION 3(ii): FISHER'S EXACT TEST
# ============================================================


# ------------------------------------------------------------
# 11. FISHER'S EXACT TEST FOR EACH SPECIES
# ------------------------------------------------------------

# Contingency table:
#
#                 Present    Absent
# Control            a         b
# Diseased           c         d
#
# The odds ratio returned by fisher.test() is:
# (a*d)/(b*c)
#
# OR > 1: higher odds of presence in Control
# OR < 1: higher odds of presence in Diseased


fisher_results <- lapply(
  seq_len(nrow(filtered_prevalence)),
  function(i) {
    
    row <- filtered_prevalence[i, ]
    
    contingency_table <- matrix(
      c(
        row$Control_Present,
        row$Control_Absent,
        row$Diseased_Present,
        row$Diseased_Absent
      ),
      nrow = 2,
      byrow = TRUE
    )
    
    test <- fisher.test(
      contingency_table,
      alternative = "two.sided"
    )
    
    data.frame(
      
      Species = row$Species,
      
      Odds_Ratio = unname(test$estimate),
      
      CI_Lower = unname(test$conf.int[1]),
      
      CI_Upper = unname(test$conf.int[2]),
      
      P_Value = test$p.value
    )
  }
)

fisher_results <- bind_rows(fisher_results)


# ------------------------------------------------------------
# 12. MULTIPLE TESTING CORRECTION
# ------------------------------------------------------------

# Benjamini-Hochberg correction

fisher_results <- fisher_results %>%
  mutate(
    
    Adjusted_P_Value = p.adjust(
      P_Value,
      method = "BH"
    ),
    
    Significance = case_when(
      Adjusted_P_Value < 0.001 ~ "***",
      Adjusted_P_Value < 0.01 ~ "**",
      Adjusted_P_Value < 0.05 ~ "*",
      TRUE ~ "NS"
    ),
    
    Significant = Adjusted_P_Value < 0.05
  )


# ------------------------------------------------------------
# 13. COMBINE PREVALENCE AND STATISTICAL RESULTS
# ------------------------------------------------------------

final_prevalence_results <- filtered_prevalence %>%
  left_join(
    fisher_results,
    by = "Species"
  ) %>%
  arrange(Adjusted_P_Value)

significant_species <- final_prevalence_results %>%
  filter(Significant) %>%
  arrange(Adjusted_P_Value)

cat("\nNumber of significant species (FDR < 0.05):",
    nrow(significant_species), "\n")

write.csv(
  final_prevalence_results,
  "Question3_Results/Species_Prevalence_Statistical_Results.csv",
  row.names = FALSE
)

write.csv(
  significant_species,
  "Question3_Results/Significant_Species.csv",
  row.names = FALSE
)


# ============================================================
# QUESTION 3(iii): TOP POSITIVE AND NEGATIVE SPECIES
# ============================================================


# ------------------------------------------------------------
# 14. RANK SPECIES BY PREVALENCE DIFFERENCE
# ------------------------------------------------------------

# Positive difference = higher prevalence in Diseased
# Negative difference = higher prevalence in Control

# Top 3 positive among all retained species

top_positive_all <- final_prevalence_results %>%
  arrange(desc(Prevalence_Difference_PP)) %>%
  slice_head(n = 3)

# Top 3 negative among all retained species

top_negative_all <- final_prevalence_results %>%
  arrange(Prevalence_Difference_PP) %>%
  slice_head(n = 3)

# Top 3 positive among statistically significant species

top_positive_significant <- significant_species %>%
  filter(Prevalence_Difference_PP > 0) %>%
  arrange(desc(Prevalence_Difference_PP)) %>%
  slice_head(n = 3)

# Top 3 negative among statistically significant species

top_negative_significant <- significant_species %>%
  filter(Prevalence_Difference_PP < 0) %>%
  arrange(Prevalence_Difference_PP) %>%
  slice_head(n = 3)


cat("\nTOP 3 POSITIVE PREVALENCE DIFFERENCES:\n")
print(top_positive_all)

cat("\nTOP 3 NEGATIVE PREVALENCE DIFFERENCES:\n")
print(top_negative_all)

cat("\nTOP 3 POSITIVE SIGNIFICANT SPECIES:\n")
print(top_positive_significant)

cat("\nTOP 3 NEGATIVE SIGNIFICANT SPECIES:\n")
print(top_negative_significant)


# Save ranking tables

write.csv(
  top_positive_all,
  "Question3_Results/Top_3_Positive_All.csv",
  row.names = FALSE
)

write.csv(
  top_negative_all,
  "Question3_Results/Top_3_Negative_All.csv",
  row.names = FALSE
)

write.csv(
  top_positive_significant,
  "Question3_Results/Top_3_Positive_Significant.csv",
  row.names = FALSE
)

write.csv(
  top_negative_significant,
  "Question3_Results/Top_3_Negative_Significant.csv",
  row.names = FALSE
)


# ============================================================
# QUESTION 3(iv): PREVALENCE VISUALIZATION
# ============================================================


# ------------------------------------------------------------
# 15. SELECT SIX SPECIES FOR VISUALIZATION
# ------------------------------------------------------------

# Top 3 positive and top 3 negative species
# based on prevalence difference among retained species.

selected_species <- bind_rows(
  top_positive_all,
  top_negative_all
) %>%
  distinct(Species, .keep_all = TRUE) %>%
  mutate(
    Direction = ifelse(
      Prevalence_Difference_PP > 0,
      "Higher in Diseased",
      "Higher in Control"
    ),
    
    Species_Label = ifelse(
      Significant,
      paste0(Species, " *"),
      Species
    )
  ) %>%
  arrange(Prevalence_Difference_PP)

selected_species$Species_Label <- factor(
  selected_species$Species_Label,
  levels = selected_species$Species_Label
)


# ------------------------------------------------------------
# 16. CONVERT TO LONG FORMAT
# ------------------------------------------------------------

plot_data <- selected_species %>%
  select(
    Species_Label,
    Control_Prevalence_Percent,
    Diseased_Prevalence_Percent,
    Significant
  ) %>%
  pivot_longer(
    cols = c(
      Control_Prevalence_Percent,
      Diseased_Prevalence_Percent
    ),
    names_to = "Condition",
    values_to = "Prevalence"
  ) %>%
  mutate(
    
    # Use case_when instead of recode to avoid namespace conflicts
    
    Condition = case_when(
      Condition == "Control_Prevalence_Percent" ~ "Control",
      Condition == "Diseased_Prevalence_Percent" ~ "Diseased",
      TRUE ~ Condition
    ),
    
    Condition = factor(
      Condition,
      levels = c("Control", "Diseased")
    )
  )


# ------------------------------------------------------------
# 17. CREATE GROUPED PREVALENCE BAR PLOT
# ------------------------------------------------------------

prevalence_plot <- ggplot(
  plot_data,
  aes(
    x = Species_Label,
    y = Prevalence,
    fill = Condition
  )
) +
  
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  
  coord_flip() +
  
  scale_y_continuous(
    labels = scales::label_number(suffix = "%"),
    limits = c(0, 105)
  ) +
  
  labs(
    title = "Species Prevalence in Control and Diseased Groups",
    subtitle = "* Indicates statistically significant prevalence difference (FDR < 0.05)",
    x = "Microbial Species",
    y = "Prevalence (%)",
    fill = "Study Condition"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    axis.text.y = element_text(
      size = 9
    ),
    legend.position = "top"
  )

print(prevalence_plot)

ggsave(
  "Question3_Results/Top_6_Species_Prevalence.png",
  prevalence_plot,
  width = 11,
  height = 7,
  dpi = 300
)


# ------------------------------------------------------------
# 18. VISUALIZE PREVALENCE DIFFERENCE
# ------------------------------------------------------------

# Select the 20 largest positive and 20 largest negative differences

top_difference_plot_data <- bind_rows(
  final_prevalence_results %>%
    arrange(desc(Prevalence_Difference_PP)) %>%
    slice_head(n = 20),
  
  final_prevalence_results %>%
    arrange(Prevalence_Difference_PP) %>%
    slice_head(n = 20)
) %>%
  distinct(Species, .keep_all = TRUE) %>%
  mutate(
    Direction = ifelse(
      Prevalence_Difference_PP > 0,
      "Higher in Diseased",
      "Higher in Control"
    )
  ) %>%
  arrange(Prevalence_Difference_PP)

top_difference_plot_data$Species <- factor(
  top_difference_plot_data$Species,
  levels = top_difference_plot_data$Species
)

difference_plot <- ggplot(
  top_difference_plot_data,
  aes(
    x = Species,
    y = Prevalence_Difference_PP,
    fill = Direction
  )
) +
  
  geom_col(width = 0.75) +
  
  coord_flip() +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  
  labs(
    title = "Largest Differences in Species Prevalence",
    subtitle = "Positive values indicate higher prevalence in Diseased samples",
    x = "Microbial Species",
    y = "Prevalence Difference (percentage points)",
    fill = "Direction"
  ) +
  
  theme_classic(base_size = 11) +
  
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    legend.position = "top"
  )

print(difference_plot)

ggsave(
  "Question3_Results/Top_40_Prevalence_Differences.png",
  difference_plot,
  width = 11,
  height = 10,
  dpi = 300
)


# ============================================================
# QUESTION 3(v): PREVALENCE VS MEAN ABUNDANCE
# ============================================================


# ------------------------------------------------------------
# 19. CALCULATE MEAN ABUNDANCE IN EACH GROUP
# ------------------------------------------------------------

# Mean abundance across ALL samples, including zeros.
# This reflects both species prevalence and abundance.

control_matrix <- species_matrix[
  control_samples,
  ,
  drop = FALSE
]

diseased_matrix <- species_matrix[
  diseased_samples,
  ,
  drop = FALSE
]

control_mean_abundance <- colMeans(control_matrix)

diseased_mean_abundance <- colMeans(diseased_matrix)


# ------------------------------------------------------------
# 20. CALCULATE MEAN ABUNDANCE AMONG POSITIVE SAMPLES
# ------------------------------------------------------------

# This calculates abundance only among samples
# in which the species is detected.

mean_positive_abundance <- function(x) {
  
  apply(
    x,
    2,
    function(z) {
      
      positive_values <- z[z > 0]
      
      if (length(positive_values) == 0) {
        return(NA_real_)
      }
      
      mean(positive_values)
    }
  )
}

control_positive_mean <- mean_positive_abundance(
  control_matrix
)

diseased_positive_mean <- mean_positive_abundance(
  diseased_matrix
)


# ------------------------------------------------------------
# 21. CREATE ABUNDANCE COMPARISON TABLE
# ------------------------------------------------------------

abundance_comparison <- data.frame(
  
  Species = colnames(species_matrix),
  
  Control_Mean_Abundance =
    as.numeric(control_mean_abundance),
  
  Diseased_Mean_Abundance =
    as.numeric(diseased_mean_abundance),
  
  Mean_Abundance_Difference =
    as.numeric(
      diseased_mean_abundance -
        control_mean_abundance
    ),
  
  Control_Positive_Mean_Abundance =
    as.numeric(control_positive_mean),
  
  Diseased_Positive_Mean_Abundance =
    as.numeric(diseased_positive_mean),
  
  Positive_Mean_Abundance_Difference =
    as.numeric(
      diseased_positive_mean -
        control_positive_mean
    )
)


# ------------------------------------------------------------
# 22. COMBINE PREVALENCE AND MEAN ABUNDANCE
# ------------------------------------------------------------

prevalence_abundance_comparison <- final_prevalence_results %>%
  left_join(
    abundance_comparison,
    by = "Species"
  ) %>%
  arrange(Adjusted_P_Value)

write.csv(
  abundance_comparison,
  "Question3_Results/Mean_Abundance_Comparison_All_Species.csv",
  row.names = FALSE
)

write.csv(
  prevalence_abundance_comparison,
  "Question3_Results/Prevalence_vs_Mean_Abundance.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 23. CORRELATION BETWEEN PREVALENCE AND MEAN ABUNDANCE
# ------------------------------------------------------------

# Compare:
# 1. Prevalence difference vs mean abundance difference
# 2. Prevalence difference vs positive-sample mean abundance difference

correlation_data <- prevalence_abundance_comparison %>%
  filter(
    is.finite(Prevalence_Difference_PP),
    is.finite(Mean_Abundance_Difference),
    is.finite(Positive_Mean_Abundance_Difference)
  )

prevalence_abundance_correlation <- cor.test(
  correlation_data$Prevalence_Difference_PP,
  correlation_data$Mean_Abundance_Difference,
  method = "spearman",
  exact = FALSE
)

positive_abundance_correlation <- cor.test(
  correlation_data$Prevalence_Difference_PP,
  correlation_data$Positive_Mean_Abundance_Difference,
  method = "spearman",
  exact = FALSE
)

cat("\nCORRELATION: PREVALENCE DIFFERENCE VS MEAN ABUNDANCE DIFFERENCE\n")
print(prevalence_abundance_correlation)

cat("\nCORRELATION: PREVALENCE DIFFERENCE VS POSITIVE-SAMPLE MEAN ABUNDANCE DIFFERENCE\n")
print(positive_abundance_correlation)


# ------------------------------------------------------------
# 24. SAVE CORRELATION RESULTS
# ------------------------------------------------------------

correlation_summary <- data.frame(
  
  Comparison = c(
    "Prevalence difference vs mean abundance difference",
    "Prevalence difference vs positive-sample mean abundance difference"
  ),
  
  Spearman_Rho = c(
    unname(prevalence_abundance_correlation$estimate),
    unname(positive_abundance_correlation$estimate)
  ),
  
  P_Value = c(
    prevalence_abundance_correlation$p.value,
    positive_abundance_correlation$p.value
  )
)

write.csv(
  correlation_summary,
  "Question3_Results/Prevalence_Abundance_Correlation.csv",
  row.names = FALSE
)

print(correlation_summary)


# ------------------------------------------------------------
# 25. CREATE PREVALENCE VS MEAN ABUNDANCE SCATTER PLOT
# ------------------------------------------------------------

abundance_plot <- ggplot(
  prevalence_abundance_comparison,
  aes(
    x = Prevalence_Difference_PP,
    y = Mean_Abundance_Difference,
    color = Significant
  )
) +
  
  geom_point(
    alpha = 0.7,
    size = 2
  ) +
  
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black"
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  
  labs(
    title = "Species Prevalence Difference vs Mean Abundance Difference",
    subtitle = "Mean abundance calculated across all samples, including zeros",
    x = "Prevalence Difference (Diseased - Control, percentage points)",
    y = "Mean Abundance Difference (Diseased - Control)",
    color = "FDR Significant"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    legend.position = "top"
  )

print(abundance_plot)

ggsave(
  "Question3_Results/Prevalence_vs_Mean_Abundance.png",
  abundance_plot,
  width = 9,
  height = 7,
  dpi = 300
)


# ============================================================
# 26. SAVE COMPLETE QUESTION 3 RESULTS
# ============================================================

capture.output(
  {
    
    cat("============================================\n")
    cat("QUESTION 3: SPECIES PREVALENCE ANALYSIS\n")
    cat("============================================\n")
    
    cat("\nSAMPLE COUNTS:\n")
    print(table(metadata$study_condition))
    
    cat("\nNUMBER OF SPECIES BEFORE FILTERING:\n")
    print(nrow(prevalence_table))
    
    cat("\nNUMBER OF SPECIES AFTER FILTERING:\n")
    print(nrow(filtered_prevalence))
    
    cat("\nNUMBER OF SIGNIFICANT SPECIES (FDR < 0.05):\n")
    print(nrow(significant_species))
    
    cat("\nTOP 3 POSITIVE SPECIES:\n")
    print(top_positive_all)
    
    cat("\nTOP 3 NEGATIVE SPECIES:\n")
    print(top_negative_all)
    
    cat("\nTOP 3 POSITIVE SIGNIFICANT SPECIES:\n")
    print(top_positive_significant)
    
    cat("\nTOP 3 NEGATIVE SIGNIFICANT SPECIES:\n")
    print(top_negative_significant)
    
    cat("\nPREVALENCE VS MEAN ABUNDANCE CORRELATION:\n")
    print(prevalence_abundance_correlation)
    
    cat("\nPREVALENCE VS POSITIVE-SAMPLE MEAN ABUNDANCE CORRELATION:\n")
    print(positive_abundance_correlation)
    
  },
  file = "Question3_Results/Complete_Question3_Results.txt"
)


# ------------------------------------------------------------
# 27. COMPLETION MESSAGE
# ------------------------------------------------------------

cat("\n============================================\n")
cat("QUESTION 3 ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("============================================\n")

cat("\nAll results saved in:\n")
cat("Question3_Results/\n")
