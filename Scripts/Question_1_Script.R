setwd("C:/SUNANDINI/Projects/Comeg_ass_2")


# ============================================================
# QUESTION 1: MICROBIAL ALPHA DIVERSITY ANALYSIS
# Shannon Diversity and Pielou's Evenness
# ============================================================

# ------------------------------------------------------------
# 1. INSTALL AND LOAD REQUIRED PACKAGES
# ------------------------------------------------------------

packages <- c(
  "tidyverse",
  "vegan",
  "rstatix",
  "car",
  "ggpubr"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(tidyverse)
library(vegan)
library(rstatix)
library(car)
library(ggpubr)

# ------------------------------------------------------------
# 2. READ THE DATA
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

# Output directory
dir.create("Question1_Results", showWarnings = FALSE)

# ------------------------------------------------------------
# 3. MATCH SAMPLE IDS BETWEEN THE TWO DATASETS
# ------------------------------------------------------------

# Rename first column of species profile as sample_id
colnames(species_profile)[1] <- "sample_id"

# Ensure sample IDs are character
metadata$sample_id <- as.character(metadata$sample_id)
species_profile$sample_id <- as.character(species_profile$sample_id)

# Check duplicate sample IDs
stopifnot(!anyDuplicated(metadata$sample_id))
stopifnot(!anyDuplicated(species_profile$sample_id))

# Identify common samples
common_samples <- intersect(
  metadata$sample_id,
  species_profile$sample_id
)

cat("Common samples:", length(common_samples), "\n")

# Retain common samples and maintain identical ordering
metadata <- metadata %>%
  filter(sample_id %in% common_samples) %>%
  arrange(match(sample_id, common_samples))

species_profile <- species_profile %>%
  filter(sample_id %in% common_samples) %>%
  arrange(match(sample_id, common_samples))

# Check sample ordering
stopifnot(
  identical(metadata$sample_id, species_profile$sample_id)
)

# ------------------------------------------------------------
# 4. PREPARE SPECIES ABUNDANCE MATRIX
# ------------------------------------------------------------

species_matrix <- species_profile %>%
  select(-sample_id) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

rownames(species_matrix) <- species_profile$sample_id

# Check abundance values
if (anyNA(species_matrix) || any(!is.finite(species_matrix))) {
  stop("Species matrix contains NA, NaN or Inf values.")
}

if (any(species_matrix < 0)) {
  stop("Negative abundance values detected.")
}

# Remove species that are zero across all samples
species_matrix <- species_matrix[
  ,
  colSums(species_matrix) > 0,
  drop = FALSE
]

cat("Number of species retained:", ncol(species_matrix), "\n")

# Remove samples with zero total abundance
sample_totals <- rowSums(species_matrix)

if (any(sample_totals == 0)) {
  warning("Samples with zero total abundance are being removed.")
  
  valid_samples <- names(sample_totals[sample_totals > 0])
  
  species_matrix <- species_matrix[valid_samples, , drop = FALSE]
  
  metadata <- metadata %>%
    filter(sample_id %in% valid_samples) %>%
    arrange(match(sample_id, valid_samples))
}

# ------------------------------------------------------------
# 5. CALCULATE SHANNON DIVERSITY
# ------------------------------------------------------------

# Shannon index
shannon <- vegan::diversity(
  species_matrix,
  index = "shannon"
)

# Species richness (number of species with abundance > 0)
richness <- rowSums(species_matrix > 0)

# Pielou's evenness
# J = H / log(S)
pielou <- ifelse(
  richness > 1,
  shannon / log(richness),
  NA_real_
)

# ------------------------------------------------------------
# 6. CREATE SAMPLE-LEVEL DIVERSITY TABLE
# ------------------------------------------------------------

diversity_table <- metadata %>%
  mutate(
    Shannon = as.numeric(shannon[sample_id]),
    Richness = as.numeric(richness[sample_id]),
    Pielou_Evenness = as.numeric(pielou[sample_id])
  ) %>%
  filter(
    !is.na(Shannon),
    !is.na(Pielou_Evenness),
    !is.na(study_condition),
    !is.na(age_category)
  ) %>%
  mutate(
    study_condition = factor(
      study_condition,
      levels = c("Control", "Diseased")
    ),
    age_category = factor(
      age_category,
      levels = c("young", "middle age", "senior")
    )
  )

# Save sample-level table
write.csv(
  diversity_table,
  "Question1_Results/Sample_Level_Diversity.csv",
  row.names = FALSE
)

cat("\nSample-level diversity summary:\n")
print(head(diversity_table))

# ------------------------------------------------------------
# 7. DESCRIPTIVE STATISTICS
# ------------------------------------------------------------

group_summary <- diversity_table %>%
  group_by(study_condition, age_category) %>%
  summarise(
    N = n(),
    
    Shannon_Mean = mean(Shannon),
    Shannon_SD = sd(Shannon),
    Shannon_Median = median(Shannon),
    Shannon_IQR = IQR(Shannon),
    
    Pielou_Mean = mean(Pielou_Evenness),
    Pielou_SD = sd(Pielou_Evenness),
    Pielou_Median = median(Pielou_Evenness),
    
    .groups = "drop"
  )

print(group_summary)

write.csv(
  group_summary,
  "Question1_Results/Group_Summary.csv",
  row.names = FALSE
)

# ============================================================
# QUESTION 1(ii): STATISTICAL COMPARISONS
# ============================================================

# ------------------------------------------------------------
# 8. CONTROL VS DISEASED
# ------------------------------------------------------------

# Wilcoxon rank-sum test
condition_test <- wilcox.test(
  Shannon ~ study_condition,
  data = diversity_table,
  exact = FALSE
)

cat("\nCONTROL VS DISEASED:\n")
print(condition_test)

# Group-wise sample sizes and medians
condition_summary <- diversity_table %>%
  group_by(study_condition) %>%
  summarise(
    N = n(),
    Median_Shannon = median(Shannon),
    Mean_Shannon = mean(Shannon),
    IQR_Shannon = IQR(Shannon),
    .groups = "drop"
  )

print(condition_summary)

# ------------------------------------------------------------
# 9. SHANNON DIVERSITY ACROSS AGE CATEGORIES
# ------------------------------------------------------------

# Kruskal-Wallis test
age_test <- kruskal.test(
  Shannon ~ age_category,
  data = diversity_table
)

cat("\nACROSS AGE CATEGORIES:\n")
print(age_test)

# Dunn's post-hoc test with multiple testing correction
age_posthoc <- diversity_table %>%
  dunn_test(
    Shannon ~ age_category,
    p.adjust.method = "BH"
  )

cat("\nAGE CATEGORY POST-HOC COMPARISONS:\n")
print(age_posthoc)

write.csv(
  age_posthoc,
  "Question1_Results/Age_Posthoc_Dunn.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. AGE CATEGORY COMPARISON WITHIN EACH CONDITION
# ------------------------------------------------------------

age_within_condition <- diversity_table %>%
  group_by(study_condition) %>%
  group_modify(~ {
    dunn_test(
      .x,
      Shannon ~ age_category,
      p.adjust.method = "BH"
    )
  }) %>%
  ungroup()

print(age_within_condition)

write.csv(
  age_within_condition,
  "Question1_Results/Age_Within_Condition_Posthoc.csv",
  row.names = FALSE
)

# ============================================================
# QUESTION 1(iii): AGE × CONDITION INTERACTION
# ============================================================

# ------------------------------------------------------------
# 11. TWO-WAY FACTORIAL ANOVA
# ------------------------------------------------------------

# Type III ANOVA
# Sum contrasts are used for Type III tests.

options(contrasts = c("contr.sum", "contr.poly"))

interaction_model <- lm(
  Shannon ~ study_condition * age_category,
  data = diversity_table
)

interaction_results <- car::Anova(
  interaction_model,
  type = 3
)

cat("\nTWO-WAY ANOVA RESULTS:\n")
print(interaction_results)

# Extract interaction p-value
interaction_p <- interaction_results[
  "study_condition:age_category",
  "Pr(>F)"
]

cat(
  "\nAge × Condition interaction p-value:",
  interaction_p,
  "\n"
)

# ------------------------------------------------------------
# 12. CHECK MODEL ASSUMPTIONS
# ------------------------------------------------------------

# Normality of residuals
residual_normality <- shapiro.test(
  residuals(interaction_model)
)

cat("\nSHAPIRO-WILK TEST:\n")
print(residual_normality)

# Homogeneity of variance
variance_test <- car::leveneTest(
  Shannon ~ interaction(study_condition, age_category),
  data = diversity_table
)

cat("\nLEVENE'S TEST:\n")
print(variance_test)

# Save ANOVA table
anova_output <- as.data.frame(interaction_results)

write.csv(
  anova_output,
  "Question1_Results/Two_Way_ANOVA.csv"
)

# ------------------------------------------------------------
# 13. INTERACTION PLOT
# ------------------------------------------------------------

interaction_plot <- ggplot(
  diversity_table,
  aes(
    x = age_category,
    y = Shannon,
    color = study_condition,
    group = study_condition
  )
) +
  stat_summary(
    fun = mean,
    geom = "line",
    linewidth = 1
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 3
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.15
  ) +
  labs(
    title = "Age and Study Condition Interaction",
    subtitle = paste(
      "Age × Condition interaction p-value =",
      signif(interaction_p, 4)
    ),
    x = "Age Category",
    y = "Mean Shannon Diversity",
    color = "Study Condition"
  ) +
  theme_classic(base_size = 13)

print(interaction_plot)

ggsave(
  "Question1_Results/Interaction_Plot.png",
  interaction_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# ============================================================
# QUESTION 1(iv): SAMPLE-LEVEL DIVERSITY FIGURE
# ============================================================

# ------------------------------------------------------------
# 14. CREATE GROUP LABELS WITH SAMPLE SIZES
# ------------------------------------------------------------

diversity_table <- diversity_table %>%
  mutate(
    Group = interaction(
      study_condition,
      age_category,
      sep = " - "
    )
  )

group_labels <- diversity_table %>%
  group_by(study_condition, age_category) %>%
  summarise(
    N = n(),
    .groups = "drop"
  ) %>%
  mutate(
    x_label = paste0(age_category, "\n(n=", N, ")")
  )

# ------------------------------------------------------------
# 15. SAMPLE-LEVEL SHANNON DIVERSITY FIGURE
# ------------------------------------------------------------

# Calculate global statistical p-values
p_condition <- condition_test$p.value
p_age <- age_test$p.value

diversity_plot <- ggplot(
  diversity_table,
  aes(
    x = age_category,
    y = Shannon,
    fill = study_condition
  )
) +
  
  # Boxplots
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.65
  ) +
  
  # Individual sample points
  geom_point(
    aes(color = study_condition),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    size = 1.8,
    alpha = 0.55
  ) +
  
  scale_x_discrete(
    labels = setNames(
      group_labels$x_label,
      group_labels$age_category
    )
  ) +
  
  labs(
    title = "Microbial Alpha Diversity Across Age and Study Conditions",
    subtitle = paste0(
      "Control vs Diseased: Wilcoxon p = ",
      signif(p_condition, 3),
      " | Age: Kruskal-Wallis p = ",
      signif(p_age, 3)
    ),
    x = "Age Category (sample size shown below)",
    y = "Shannon Diversity Index",
    fill = "Study Condition",
    color = "Study Condition"
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    legend.position = "top",
    axis.text.x = element_text(
      face = "bold"
    )
  )

print(diversity_plot)

ggsave(
  "Question1_Results/Shannon_Diversity_Plot.png",
  diversity_plot,
  width = 10,
  height = 7,
  dpi = 300
)

# ------------------------------------------------------------
# 16. PIELOU'S EVENNESS FIGURE
# ------------------------------------------------------------

pielou_plot <- ggplot(
  diversity_table,
  aes(
    x = age_category,
    y = Pielou_Evenness,
    fill = study_condition
  )
) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA
  ) +
  geom_point(
    aes(color = study_condition),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    alpha = 0.5,
    size = 1.5
  ) +
  labs(
    title = "Pielou's Evenness Across Study Conditions and Age",
    x = "Age Category",
    y = "Pielou's Evenness",
    fill = "Study Condition",
    color = "Study Condition"
  ) +
  theme_classic(base_size = 13)

print(pielou_plot)

ggsave(
  "Question1_Results/Pielou_Evenness_Plot.png",
  pielou_plot,
  width = 10,
  height = 7,
  dpi = 300
)

# ============================================================
# QUESTION 1(v): HIGHEST AND LOWEST MEDIAN DIVERSITY
# ============================================================

# ------------------------------------------------------------
# 17. IDENTIFY HIGHEST AND LOWEST MEDIAN GROUPS
# ------------------------------------------------------------

median_comparison <- diversity_table %>%
  group_by(study_condition, age_category) %>%
  summarise(
    N = n(),
    Median_Shannon = median(Shannon),
    Mean_Shannon = mean(Shannon),
    IQR_Shannon = IQR(Shannon),
    .groups = "drop"
  ) %>%
  arrange(Median_Shannon)

cat("\nGROUPS RANKED BY MEDIAN SHANNON DIVERSITY:\n")
print(median_comparison)

lowest_group <- median_comparison %>%
  slice_min(
    Median_Shannon,
    n = 1,
    with_ties = TRUE
  )

highest_group <- median_comparison %>%
  slice_max(
    Median_Shannon,
    n = 1,
    with_ties = TRUE
  )

cat("\nLOWEST MEDIAN SHANNON DIVERSITY:\n")
print(lowest_group)

cat("\nHIGHEST MEDIAN SHANNON DIVERSITY:\n")
print(highest_group)

write.csv(
  median_comparison,
  "Question1_Results/Median_Diversity_Ranking.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 18. STATISTICAL COMPARISON OF HIGHEST AND LOWEST GROUPS
# ------------------------------------------------------------

if (nrow(lowest_group) == 1 && nrow(highest_group) == 1) {
  
  high_condition <- as.character(highest_group$study_condition)
  high_age <- as.character(highest_group$age_category)
  
  low_condition <- as.character(lowest_group$study_condition)
  low_age <- as.character(lowest_group$age_category)
  
  high_low_data <- diversity_table %>%
    filter(
      (study_condition == high_condition &
         age_category == high_age) |
        (study_condition == low_condition &
           age_category == low_age)
    ) %>%
    mutate(
      comparison_group = ifelse(
        study_condition == high_condition &
          age_category == high_age,
        "Highest median group",
        "Lowest median group"
      )
    )
  
  high_low_test <- wilcox.test(
    Shannon ~ comparison_group,
    data = high_low_data,
    exact = FALSE
  )
  
  cat("\nHIGHEST VS LOWEST MEDIAN GROUP:\n")
  print(high_low_test)
  
} else {
  
  cat(
    "\nTied highest/lowest median groups detected.",
    "A unique two-group comparison cannot be selected automatically.\n"
  )
}

# ============================================================
# 19. SAVE ALL STATISTICAL RESULTS
# ============================================================

capture.output(
  {
    cat("CONTROL VS DISEASED\n")
    print(condition_test)
    
    cat("\nAGE CATEGORY KRUSKAL-WALLIS\n")
    print(age_test)
    
    cat("\nAGE POST-HOC DUNN TEST\n")
    print(age_posthoc)
    
    cat("\nTWO-WAY ANOVA\n")
    print(interaction_results)
    
    cat("\nSHAPIRO-WILK\n")
    print(residual_normality)
    
    cat("\nLEVENE'S TEST\n")
    print(variance_test)
    
    cat("\nHIGHEST AND LOWEST MEDIAN GROUPS\n")
    print(highest_group)
    print(lowest_group)
    
    if (exists("high_low_test")) {
      cat("\nHIGHEST VS LOWEST GROUP TEST\n")
      print(high_low_test)
    }
  },
  file = "Question1_Results/Statistical_Results.txt"
)

cat("\n========================================\n")
cat("QUESTION 1 ANALYSIS COMPLETED!\n")
cat("Check the Question1_Results folder.\n")
cat("========================================\n")

