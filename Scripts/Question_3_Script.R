## CoMeg Assignment 2 - Q3
## Differential prevalence of species between Control and Diseased gut samples
##
## run from the repo root (the folder that has Question/ and Scripts/)
## input  : Question/metadata_all.xlsx.csv, Question/species_profile_norm.xlsx.csv
## output : Question3_Results/

suppressPackageStartupMessages({
  for (p in c("tidyverse", "scales")) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  }
  library(tidyverse)
  library(scales)
})

res_dir <- "Question3_Results"
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# small helper so every write.csv call isn't repeated in full
save_tab <- function(df, fname) {
  write.csv(df, file.path(res_dir, fname), row.names = FALSE)
}


#### ---- reading + sanity checks ---- ####

meta <- read.csv("Question/metadata_all.xlsx.csv",
                 check.names = FALSE, stringsAsFactors = FALSE)
abund <- read.csv("Question/species_profile_norm.xlsx.csv",
                  check.names = FALSE, stringsAsFactors = FALSE)

names(abund)[1] <- "sample_id"     # first col of the profile has the IDs
meta$sample_id  <- as.character(meta$sample_id)
abund$sample_id <- as.character(abund$sample_id)

need <- c("sample_id", "study_condition")
if (!all(need %in% names(meta))) {
  stop("metadata is missing: ",
       paste(setdiff(need, names(meta)), collapse = ", "))
}
if (anyDuplicated(meta$sample_id))  stop("repeated sample_id found in metadata")
if (anyDuplicated(abund$sample_id)) stop("repeated sample_id found in species table")

# keep only samples present in both files and labelled Control / Diseased,
# in the metadata order, then line the abundance rows up the same way
shared <- meta$sample_id[meta$sample_id %in% abund$sample_id]

meta <- meta %>%
  filter(sample_id %in% shared,
         study_condition %in% c("Control", "Diseased")) %>%
  arrange(match(sample_id, shared))

abund <- abund %>%
  filter(sample_id %in% meta$sample_id) %>%
  arrange(match(sample_id, meta$sample_id))

stopifnot(identical(meta$sample_id, abund$sample_id))
meta$study_condition <- factor(meta$study_condition,
                               levels = c("Control", "Diseased"))

message("samples kept: ", nrow(meta))
print(table(meta$study_condition))

# numeric matrix, samples x species
X <- abund %>%
  select(-sample_id) %>%
  mutate(across(everything(), ~ as.numeric(as.character(.x)))) %>%
  as.matrix()
rownames(X) <- abund$sample_id

if (anyNA(X) || any(!is.finite(X))) stop("abundance table has NA/NaN/Inf")
if (any(X < 0)) stop("abundance table has negative entries")

X <- X[, colSums(X) > 0, drop = FALSE]   # drop species never observed
message("species with non-zero abundance somewhere: ", ncol(X))


#### ---- (i) prevalence per group ---- ####

detected <- (X > 0) * 1          # 1 = detected, 0 = not detected
dimnames(detected) <- dimnames(X)
cat("presence/absence matrix:", dim(detected), "\n")

ids_ctrl <- meta$sample_id[meta$study_condition == "Control"]
ids_dis  <- meta$sample_id[meta$study_condition == "Diseased"]
nC <- length(ids_ctrl)
nD <- length(ids_dis)
cat("Control n =", nC, "| Diseased n =", nD, "\n")

hitsC <- colSums(detected[ids_ctrl, , drop = FALSE])
hitsD <- colSums(detected[ids_dis,  , drop = FALSE])

# difference is Diseased minus Control (percentage points),
# so > 0 means the species shows up more often in Diseased
prev_all <- tibble(
  Species                     = colnames(X),
  Control_Present             = as.numeric(hitsC),
  Control_Absent              = nC - as.numeric(hitsC),
  Control_Prevalence          = as.numeric(hitsC) / nC,
  Control_Prevalence_Percent  = 100 * as.numeric(hitsC) / nC,
  Diseased_Present            = as.numeric(hitsD),
  Diseased_Absent             = nD - as.numeric(hitsD),
  Diseased_Prevalence         = as.numeric(hitsD) / nD,
  Diseased_Prevalence_Percent = 100 * as.numeric(hitsD) / nD
) %>%
  mutate(
    Prevalence_Difference_PP = Diseased_Prevalence_Percent - Control_Prevalence_Percent,
    Total_Present            = Control_Present + Diseased_Present
  )

save_tab(prev_all, "Complete_Species_Prevalence.csv")

# very rare species can't give a meaningful test -> need >= 5 detections overall
min_hits <- 5
prev_kept <- filter(prev_all, Total_Present >= min_hits)

cat("species before filter:", nrow(prev_all),
    "| kept:", nrow(prev_kept),
    "| dropped:", nrow(prev_all) - nrow(prev_kept), "\n")

save_tab(prev_kept, "Filtered_Species_Prevalence.csv")


#### ---- (ii) Fisher's exact test, BH-adjusted ---- ####

# 2x2 per species, rows = group, cols = present/absent
#            present  absent
#  Control      a       b
#  Diseased     c       d
# fisher.test gives OR = (a*d)/(b*c); OR > 1 => presence more likely in Control

run_fisher <- function(cP, cA, dP, dA) {
  ft <- fisher.test(matrix(c(cP, cA, dP, dA), nrow = 2, byrow = TRUE),
                    alternative = "two.sided")
  c(Odds_Ratio = unname(ft$estimate),
    CI_Lower   = ft$conf.int[1],
    CI_Upper   = ft$conf.int[2],
    P_Value    = ft$p.value)
}

fisher_out <- pmap(
  list(prev_kept$Control_Present, prev_kept$Control_Absent,
       prev_kept$Diseased_Present, prev_kept$Diseased_Absent),
  run_fisher
)

fisher_tab <- as.data.frame(do.call(rbind, fisher_out))
fisher_tab <- cbind(Species = prev_kept$Species, fisher_tab)

star_code <- function(q) {
  ifelse(q < 0.001, "***",
         ifelse(q < 0.01, "**",
                ifelse(q < 0.05, "*", "NS")))
}

fisher_tab$Adjusted_P_Value <- p.adjust(fisher_tab$P_Value, method = "BH")
fisher_tab$Significance     <- star_code(fisher_tab$Adjusted_P_Value)
fisher_tab$Significant      <- fisher_tab$Adjusted_P_Value < 0.05

prev_stats <- prev_kept %>%
  left_join(fisher_tab, by = "Species") %>%
  arrange(Adjusted_P_Value)

sig_hits <- prev_stats %>% filter(Significant)   # already sorted by q-value

cat("species significant at FDR < 0.05:", nrow(sig_hits), "\n")

save_tab(prev_stats, "Species_Prevalence_Statistical_Results.csv")
save_tab(sig_hits,   "Significant_Species.csv")


#### ---- (iii) biggest shifts in each direction ---- ####

# pick k species with the largest (dir = "up") or most negative (dir = "down")
# prevalence difference
pick_extremes <- function(df, dir = c("up", "down"), k = 3) {
  dir <- match.arg(dir)
  if (dir == "up") {
    df %>% arrange(desc(Prevalence_Difference_PP)) %>% slice_head(n = k)
  } else {
    df %>% arrange(Prevalence_Difference_PP) %>% slice_head(n = k)
  }
}

# over every species that passed the filter...
up3   <- pick_extremes(prev_stats, "up")
down3 <- pick_extremes(prev_stats, "down")

# ...and only among the FDR-significant ones
up3_sig   <- pick_extremes(filter(sig_hits, Prevalence_Difference_PP > 0), "up")
down3_sig <- pick_extremes(filter(sig_hits, Prevalence_Difference_PP < 0), "down")

top_lists <- list(
  "Top_3_Positive_All"         = up3,
  "Top_3_Negative_All"         = down3,
  "Top_3_Positive_Significant" = up3_sig,
  "Top_3_Negative_Significant" = down3_sig
)

for (nm in names(top_lists)) {
  cat("\n>>", gsub("_", " ", nm), "\n")
  print(top_lists[[nm]])
  save_tab(top_lists[[nm]], paste0(nm, ".csv"))
}


#### ---- (iv) plots ---- ####

## grouped bars for the 6 species above (3 up + 3 down), * = FDR significant
six <- bind_rows(up3, down3) %>%
  distinct(Species, .keep_all = TRUE) %>%
  mutate(Direction = if_else(Prevalence_Difference_PP > 0,
                             "Higher in Diseased", "Higher in Control"),
         Species_Label = if_else(Significant, paste0(Species, " *"), Species)) %>%
  arrange(Prevalence_Difference_PP)

six$Species_Label <- factor(six$Species_Label, levels = six$Species_Label)

six_long <- six %>%
  select(Species_Label, Significant,
         Control  = Control_Prevalence_Percent,
         Diseased = Diseased_Prevalence_Percent) %>%
  pivot_longer(c(Control, Diseased),
               names_to = "Condition", values_to = "Prevalence") %>%
  mutate(Condition = factor(Condition, levels = c("Control", "Diseased")))

p_six <- ggplot(six_long, aes(Species_Label, Prevalence, fill = Condition)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  coord_flip() +
  scale_y_continuous(labels = label_number(suffix = "%"), limits = c(0, 105)) +
  labs(title    = "Species Prevalence in Control and Diseased Groups",
       subtitle = "* Indicates statistically significant prevalence difference (FDR < 0.05)",
       x = "Microbial Species", y = "Prevalence (%)", fill = "Study Condition") +
  theme_classic(base_size = 12) +
  theme(plot.title  = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 9),
        legend.position = "top")

print(p_six)
ggsave(file.path(res_dir, "Top_6_Species_Prevalence.png"),
       p_six, width = 11, height = 7, dpi = 300)

## diverging bars: 20 most positive + 20 most negative differences
forty <- bind_rows(pick_extremes(prev_stats, "up",   k = 20),
                   pick_extremes(prev_stats, "down", k = 20)) %>%
  distinct(Species, .keep_all = TRUE) %>%
  mutate(Direction = if_else(Prevalence_Difference_PP > 0,
                             "Higher in Diseased", "Higher in Control")) %>%
  arrange(Prevalence_Difference_PP)

forty$Species <- factor(forty$Species, levels = forty$Species)

p_forty <- ggplot(forty, aes(Species, Prevalence_Difference_PP, fill = Direction)) +
  geom_col(width = 0.75) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title    = "Largest Differences in Species Prevalence",
       subtitle = "Positive values indicate higher prevalence in Diseased samples",
       x = "Microbial Species",
       y = "Prevalence Difference (percentage points)",
       fill = "Direction") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")

print(p_forty)
ggsave(file.path(res_dir, "Top_40_Prevalence_Differences.png"),
       p_forty, width = 11, height = 10, dpi = 300)


#### ---- (v) does prevalence track mean abundance? ---- ####

X_ctrl <- X[ids_ctrl, , drop = FALSE]
X_dis  <- X[ids_dis,  , drop = FALSE]

# two versions of "mean abundance":
#   overall  -> zeros included, so it partly reflects prevalence itself
#   detected -> only samples where the species is actually present
mean_if_detected <- function(M) {
  apply(M, 2, function(v) {
    v <- v[v > 0]
    if (length(v)) mean(v) else NA_real_
  })
}

muC  <- colMeans(X_ctrl)
muD  <- colMeans(X_dis)
muC0 <- mean_if_detected(X_ctrl)
muD0 <- mean_if_detected(X_dis)

abund_tab <- data.frame(
  Species                            = colnames(X),
  Control_Mean_Abundance             = as.numeric(muC),
  Diseased_Mean_Abundance            = as.numeric(muD),
  Mean_Abundance_Difference          = as.numeric(muD - muC),
  Control_Positive_Mean_Abundance    = as.numeric(muC0),
  Diseased_Positive_Mean_Abundance   = as.numeric(muD0),
  Positive_Mean_Abundance_Difference = as.numeric(muD0 - muC0)
)

prev_vs_abund <- prev_stats %>%
  left_join(abund_tab, by = "Species") %>%
  arrange(Adjusted_P_Value)

save_tab(abund_tab,     "Mean_Abundance_Comparison_All_Species.csv")
save_tab(prev_vs_abund, "Prevalence_vs_Mean_Abundance.csv")

# species missing from one group give NA for the detected-only mean;
# both correlations are run on the same complete-case set
cc <- prev_vs_abund %>%
  filter(if_all(c(Prevalence_Difference_PP,
                  Mean_Abundance_Difference,
                  Positive_Mean_Abundance_Difference), is.finite))

sp_cor <- function(y) {
  cor.test(cc$Prevalence_Difference_PP, y, method = "spearman", exact = FALSE)
}

cor_overall  <- sp_cor(cc$Mean_Abundance_Difference)
cor_detected <- sp_cor(cc$Positive_Mean_Abundance_Difference)

cat("\nSpearman: prevalence diff ~ mean abundance diff (zeros included)\n")
print(cor_overall)
cat("\nSpearman: prevalence diff ~ mean abundance diff (detected samples only)\n")
print(cor_detected)

cor_tab <- data.frame(
  Comparison = c("Prevalence difference vs mean abundance difference",
                 "Prevalence difference vs positive-sample mean abundance difference"),
  Spearman_Rho = c(unname(cor_overall$estimate), unname(cor_detected$estimate)),
  P_Value      = c(cor_overall$p.value, cor_detected$p.value)
)
save_tab(cor_tab, "Prevalence_Abundance_Correlation.csv")
print(cor_tab)

p_scatter <- ggplot(prev_vs_abund,
                    aes(Prevalence_Difference_PP, Mean_Abundance_Difference,
                        colour = Significant)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, colour = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title    = "Species Prevalence Difference vs Mean Abundance Difference",
       subtitle = "Mean abundance calculated across all samples, including zeros",
       x = "Prevalence Difference (Diseased - Control, percentage points)",
       y = "Mean Abundance Difference (Diseased - Control)",
       colour = "FDR Significant") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")

print(p_scatter)
ggsave(file.path(res_dir, "Prevalence_vs_Mean_Abundance.png"),
       p_scatter, width = 9, height = 7, dpi = 300)


#### ---- text summary ---- ####

sink(file.path(res_dir, "Complete_Question3_Results.txt"))
cat("Q3 - species prevalence: Control vs Diseased\n")
cat("--------------------------------------------\n")
cat("\nGroup sizes:\n");                         print(table(meta$study_condition))
cat("\nSpecies tested before filtering:\n");     print(nrow(prev_all))
cat("\nSpecies kept (>=", min_hits, "detections):\n"); print(nrow(prev_kept))
cat("\nSignificant species (BH FDR < 0.05):\n"); print(nrow(sig_hits))
for (nm in names(top_lists)) {
  cat("\n", gsub("_", " ", nm), ":\n", sep = "")
  print(top_lists[[nm]])
}
cat("\nSpearman - prevalence vs mean abundance (zeros included):\n");  print(cor_overall)
cat("\nSpearman - prevalence vs mean abundance (detected only):\n");   print(cor_detected)
sink()

cat("\nQ3 finished - outputs written to", res_dir, "\n")
