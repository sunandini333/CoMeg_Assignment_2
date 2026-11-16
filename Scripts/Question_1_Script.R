setwd("C:/SUNANDINI/Projects/Comeg_ass_2")

# Ques1: microbial alpha diversity analysis

packages <- c(
  "tidyverse",
  "vegan",
  "rstatix",
  "car",
  "ggpubr"
)

# preparing and loading env

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

#Loading data
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

#out dir
dir.create("Question1_Results", showWarnings = FALSE)

colnames(species_profile)[1] <- "sample_id"

metadata$sample_id <- as.character(metadata$sample_id)
species_profile$sample_id <- as.character(species_profile$sample_id)

stopifnot(!anyDuplicated(metadata$sample_id))
stopifnot(!anyDuplicated(species_profile$sample_id))

common_samples <- intersect(
  metadata$sample_id,
  species_profile$sample_id
)

cat("Common samples:", length(common_samples), "\n")

metadata <- metadata %>%
  filter(sample_id %in% common_samples) %>%
  arrange(match(sample_id, common_samples))

species_profile <- species_profile %>%
  filter(sample_id %in% common_samples) %>%
  arrange(match(sample_id, common_samples))

stopifnot(
  identical(metadata$sample_id, species_profile$sample_id)
)

# 4. preparing species abundance matrix
species_matrix <- species_profile %>%
  select(-sample_id) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

rownames(species_matrix) <- species_profile$sample_id

#Checking values

if (anyNA(species_matrix) || any(!is.finite(species_matrix))) {
  stop("Species matrix contains NA, NaN or Inf values.")
}
if (any(species_matrix < 0)) {
  stop("Negative abundance values detected.")
}

# Remove zero count species
species_matrix <- species_matrix[
  ,
  colSums(species_matrix) > 0,
  drop = FALSE
]
cat("Number of species retained:", ncol(species_matrix), "\n")

sample_totals <- rowSums(species_matrix)

if (any(sample_totals == 0)) {
  warning("Samples with zero total abundance are being removed.")
  
  valid_samples <- names(sample_totals[sample_totals > 0])
  
  species_matrix <- species_matrix[valid_samples, , drop = FALSE]  
  metadata <- metadata %>%
    filter(sample_id %in% valid_samples) %>%
    arrange(match(sample_id, valid_samples))
}

# shannon diversity

# shannon index
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

# sample level diversity table

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

# Saving tjhe table
write.csv(
  diversity_table,
  "Question1_Results/Sample_Level_Diversity.csv",
  row.names = FALSE
)

cat("\nSample-level diversity summary:\n")
print(head(diversity_table))

#descriptive stas

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

# Ques1 (ii): statstical comparisons

# 8. CONTROL VS DISEASED

# Wilcoxon test
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

#11. age category comparisons

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

#Ques1 (iii): age and condition interaction

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

# checking assumptions
# Normality of residuals
residual_normality <- shapiro.test(
  residuals(interaction_model)
)

cat("\nSHAPIRO-WILK TEST:\n")
print(residual_normality)

variance_test <- car::leveneTest(
  Shannon ~ interaction(study_condition, age_category),
  data = diversity_table
)

cat("\nLEVENE'S TEST:\n")
print(variance_test)

anova_output <- as.data.frame(interaction_results)

write.csv(
  anova_output,
  "Question1_Results/Two_Way_ANOVA.csv"
)

#PLots

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

# Ques1 (iv): Sample level diversity plots
# 14. group tables

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

# 15. shannon diversity figure

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
  
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.65
  ) +
  
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

# 16. pielou's evenness

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

#Ques1 (v): highest and lowest median diversity


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

# 18. statistical comparison 
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

# saving results

capture.output(
  {
    cat("control vs diseased\n")
    print(condition_test)
    
    cat("\nage category kruskal-wallis\n")
    print(age_test)
    
    cat("\nage post-hoc dunn test\n")
    print(age_posthoc)
    
    cat("\ntwo-way anova\n")
    print(interaction_results)
    
    cat("\nshapiro-wilk test\n")
    print(residual_normality)
    
    cat("\nlevene's test\n")
    print(variance_test)
    
    cat("\nhighest and lowest median groups\n")
    print(highest_group)
    print(lowest_group)
    
    if (exists("high_low_test")) {
      cat("\nhighest vs lowest group test\n")
      print(high_low_test)
    }
  },
  file = "Question1_Results/Statistical_Results.txt"
)

cat("\n========================================\n")
cat("Question 1 analysis done!!!\n")
cat("========================================\n")

