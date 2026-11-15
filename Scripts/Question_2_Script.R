setwd("/home/chila/IIITD/Sem_01/COMEG/Assignment_02")

library(readxl)
library(vegan)
library(ggplot2)

# Load data
species_raw <- read_excel("species_profile_norm.xlsx")
meta <- read_excel("metadata_all.xlsx")

# Prepare species abundance matrix
species_mat <- as.matrix(species_raw[, -1])
storage.mode(species_mat) <- "numeric"
rownames(species_mat) <- species_raw[[1]]

# Match metadata to the sample order
meta <- as.data.frame(meta)
rownames(meta) <- meta$sample_id
meta <- meta[rownames(species_mat), ]

stopifnot(identical(rownames(species_mat), meta$sample_id))

# Set metadata factors
meta$study_condition <- factor(
  meta$study_condition,
  levels = c("Control", "Diseased")
)

meta$age_category <- factor(
  meta$age_category,
  levels = c("young", "middle age", "senior")
)

meta$gender <- factor(meta$gender)

# Check for missing values
stopifnot(!anyNA(species_mat))
stopifnot(!anyNA(meta[, c("study_condition", "age_category", "gender")]))

# i. Bray-Curtis dissimilarity and PCoA

bc_dist <- vegdist(species_mat, method = "bray")

pcoa_res <- cmdscale(bc_dist, k = 2, eig = TRUE)

# Calculate the percentage of variation explained
pos_eig <- pcoa_res$eig[pcoa_res$eig > 0]
var_explained <- pos_eig / sum(pos_eig) * 100

cat(sprintf("PCo1: %.2f%%\n", var_explained[1]))
cat(sprintf("PCo2: %.2f%%\n", var_explained[2]))
cat(sprintf("Combined: %.2f%%\n", sum(var_explained[1:2])))

# Prepare plotting data
pcoa_df <- data.frame(
  sample_id = rownames(pcoa_res$points),
  PCo1 = pcoa_res$points[, 1],
  PCo2 = pcoa_res$points[, 2],
  study_condition = meta$study_condition,
  age_category = meta$age_category,
  gender = meta$gender
)

# PCoA plot
pcoa_plot <- ggplot(
  pcoa_df,
  aes(
    x = PCo1,
    y = PCo2,
    color = study_condition,
    shape = age_category
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey70", linewidth = 0.4) +
  geom_point(size = 2.8, alpha = 0.85) +
  scale_color_manual(
    values = c("Control" = "#1f77b4", "Diseased" = "#d62728")
  ) +
  scale_shape_manual(
    values = c("young" = 16, "middle age" = 15, "senior" = 17)
  ) +
  labs(
    title = "PCoA of Bray-Curtis Dissimilarity Across Individuals",
    x = sprintf("PCo1 (%.2f%%)", var_explained[1]),
    y = sprintf("PCo2 (%.2f%%)", var_explained[2]),
    color = "Study Condition",
    shape = "Age Category"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank()
  )

print(pcoa_plot)

ggsave(
  "pcoa_bray_curtis.png",
  plot = pcoa_plot,
  width = 8,
  height = 6,
  dpi = 300
)


# ii. PERMANOVA

set.seed(123)

# Sequential tests
permanova_seq <- adonis2(
  bc_dist ~ study_condition + age_category + gender,
  data = meta,
  permutations = 9999,
  by = "terms"
)

print(permanova_seq)

# Marginal tests
permanova_margin <- adonis2(
  bc_dist ~ study_condition + age_category + gender,
  data = meta,
  permutations = 9999,
  by = "margin"
)

print(permanova_margin)


# iii. Multivariate dispersion

set.seed(123)

disp_condition <- betadisper(
  bc_dist,
  group = meta$study_condition,
  type = "median"
)

print(disp_condition)

# Compare distances between groups
disp_mean <- tapply(
  disp_condition$distances,
  meta$study_condition,
  mean
)

disp_sd <- tapply(
  disp_condition$distances,
  meta$study_condition,
  sd
)

print(disp_mean)
print(disp_sd)

# Test differences in dispersion
anova(disp_condition)

set.seed(123)

permutest_condition <- permutest(
  disp_condition,
  permutations = 9999
)

print(permutest_condition)

# Plot dispersion
disp_df <- data.frame(
  study_condition = disp_condition$group,
  Distance_to_Centroid = disp_condition$distances
)

bp1 <- ggplot(
  disp_df,
  aes(
    x = study_condition,
    y = Distance_to_Centroid,
    fill = study_condition
  )
) +
  geom_boxplot(alpha = 0.65, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
  scale_fill_manual(
    values = c("Control" = "#1f77b4", "Diseased" = "#d62728")
  ) +
  labs(
    title = "Multivariate Dispersion by Study Condition",
    x = "Study Condition",
    y = "Distance to Spatial Median"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")

print(bp1)

ggsave(
  "multivariate_dispersion_study.png",
  plot = bp1,
  width = 8,
  height = 6,
  dpi = 300
)


# iv. Alternative model specifications
# 
# Primary model
set.seed(123)

perm_primary <- adonis2(
  bc_dist ~ study_condition + age_category + gender,
  data = meta,
  permutations = 9999,
  by = "terms"
)

# Alternative model: Age category entered first
set.seed(123)

perm_alt_order <- adonis2(
  bc_dist ~ age_category + study_condition + gender,
  data = meta,
  permutations = 9999,
  by = "terms"
)

print(perm_alt_order)

# Alternative model: Include an interaction between disease status and age
set.seed(123)

perm_interaction <- adonis2(
  bc_dist ~ study_condition * age_category + gender,
  data = meta,
  permutations = 9999,
  by = "terms"
)

print(perm_interaction)


# v. Compare variance explained by metadata variables

es_primary <- as.data.frame(
  perm_primary[1:3, c("Df", "SumOfSqs", "R2", "F", "Pr(>F)")]
)

es_primary$Variable <- rownames(es_primary)
es_primary$R2_percent <- es_primary$R2 * 100

# Rank variables by their R2 values
es_ranked <- es_primary[
  order(-es_primary$R2),
  c("Variable", "Df", "R2", "R2_percent", "F", "Pr(>F)")
]

print(es_ranked, row.names = FALSE)

# Plot R2 values
es_plot <- ggplot(
  es_ranked,
  aes(
    x = reorder(Variable, R2_percent),
    y = R2_percent,
    fill = Variable
  )
) +
  geom_col(width = 0.55, show.legend = FALSE) +
  geom_text(
    aes(label = sprintf("%.2f%% (p = %.4f)", R2_percent, `Pr(>F)`)),
    vjust = -0.4,
    size = 3.8
  ) +
  scale_y_continuous(
    limits = c(0, max(es_ranked$R2_percent) * 1.2),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Variance Explained by Metadata Variables",
    x = "Metadata Variable",
    y = "Variance Explained (R2, %)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(es_plot)
