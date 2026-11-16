# Comparative Analysis of Gut Microbial Diversity, Prevalence, and Co-occurrence Patterns

## CoMeg Assignment 2

---

## 1. Project Overview

This project focuses on the statistical and ecological analysis of gut microbial communities to investigate differences between Control and Diseased individuals.

The analysis uses microbial species abundance profiles and sample-level metadata to investigate microbial diversity, community composition, species prevalence, and microbial co-occurrence patterns.

Four major analyses were performed:

1. **Alpha Diversity Analysis:** Evaluation of within-sample microbial diversity using Shannon diversity and Pielou's evenness.
2. **Beta Diversity Analysis:** Investigation of between-sample microbial community differences using Bray–Curtis dissimilarity, Principal Coordinates Analysis (PCoA), and PERMANOVA.
3. **Differential Species Prevalence Analysis:** Identification of microbial species with significantly different prevalence between Control and Diseased groups.
4. **Microbial Co-occurrence Network Analysis:** Identification of significant microbial associations and network hub species using Spearman correlation.

All four analyses were implemented using R.

---

## 2. Dataset Description

Two datasets were used for the analysis.

| File | Description |
|---|---|
| `metadata_all.xlsx.csv` | Sample-level metadata containing sample IDs, gender, study condition, age category, and age |
| `species_profile_norm.xlsx.csv` | Normalized microbial species abundance matrix |

The original assignment document is included as `Assignment2.pdf`.

### Dataset Characteristics

| Parameter | Value |
|---|---:|
| Total matched samples | 265 |
| Total microbial species | 777 |
| Control samples | 133 |
| Diseased samples | 132 |

Sample identifiers were matched between the metadata and species abundance datasets before downstream analysis.

---

## 3. Repository Structure

```text
CoMeg_Assignment_2/
│
├── Question/
│   ├── Assignment2.pdf
│   ├── metadata_all.xlsx.csv
│   └── species_profile_norm.xlsx.csv
│
├── Question1_Results/
│   ├── Age_Posthoc_Dunn.csv
│   ├── Age_Within_Condition_Posthoc.csv
│   ├── Group_Summary.csv
│   ├── Interaction_Plot.png
│   ├── Median_Diversity_Ranking.csv
│   ├── Pielou_Evenness_Plot.png
│   ├── Sample_Level_Diversity.csv
│   ├── Shannon_Diversity_Plot.png
│   ├── Statistical_Results.txt
│   └── Two_Way_ANOVA.csv
│
├── Question2_Results/
│   ├── Alternative_PERMANOVA.csv
│   ├── Bray_Curtis_Distance_Matrix.csv
│   ├── Complete_Statistical_Results.txt
│   ├── Dispersion_Distances.csv
│   ├── Dispersion_Plot.png
│   ├── Dispersion_Test.txt
│   ├── Final_Analysis_Metadata.csv
│   ├── PCoA_BrayCurtis.png
│   ├── PCoA_Coordinates.csv
│   ├── Primary_PERMANOVA_Marginal.csv
│   ├── Sequential_Model_Comparison.csv
│   └── Strongest_Metadata_Association.csv
│
├── Question3_Results/
│   ├── Complete_Question3_Results.txt
│   ├── Complete_Species_Prevalence.csv
│   ├── Complete_Statistical_Results.txt
│   ├── Filtered_Species_Prevalence.csv
│   ├── Mean_Abundance_Comparison_All_Species.csv
│   ├── Prevalence_Abundance_Correlation.csv
│   ├── Prevalence_And_Abundance_Comparison.csv
│   ├── Prevalence_Difference.png
│   ├── Prevalence_vs_Mean_Abundance.csv
│   ├── Prevalence_vs_Mean_Abundance.png
│   ├── Selected_Species_Biological_Interpretation.csv
│   ├── Selected_Species_Comparison.png
│   ├── Significant_Species.csv
│   ├── Species_Prevalence_Statistical_Results.csv
│   ├── Top_3_Negative_All.csv
│   ├── Top_3_Negative_Significant.csv
│   ├── Top_3_Positive_All.csv
│   ├── Top_3_Positive_Significant.csv
│   ├── Top_40_Prevalence_Differences.png
│   └── Top_6_Species_Prevalence.png
│
├── Question4_Results/
│   ├── All_Spearman_Correlations.csv
│   ├── Alternative_All_Correlations.csv
│   ├── Alternative_Selected_Edges.csv
│   ├── Complete_Statistical_Results.txt
│   ├── Hub_Species.csv
│   ├── Microbial_Cooccurrence_Network.graphml
│   ├── Microbial_Cooccurrence_Network.png
│   ├── Negative_Associations.csv
│   ├── Network_Edge_Overlap.csv
│   ├── Network_Node_Degree.csv
│   ├── Positive_Associations.csv
│   ├── Retained_Species.csv
│   ├── Robust_Associations.csv
│   ├── Robustness_Comparison.csv
│   ├── Selected_Species_Associations.csv
│   ├── Sensitivity_Analysis.png
│   ├── Sensitivity_Analysis_Summary.csv
│   ├── Spearman_Correlation_Matrix.csv
│   ├── Species_Prevalence.csv
│   └── Top_10_Hub_Species.png
│
├── Scripts/
│   ├── Question_1_Script.R
│   ├── Question_2_Script.R
│   ├── Question_3_Script.R
│   └── Question_4_Script.R
│
└── README.md
```

### Directory Description

| Directory | Description |
|---|---|
| `Question/` | Original assignment PDF and input datasets |
| `Question1_results/` | Alpha diversity analysis results, statistical tests, and visualizations |
| `Question2_results/` | Beta diversity analysis results, Bray–Curtis distance matrix, PCoA, PERMANOVA, and PERMDISP outputs |
| `Question3_results/` | Species prevalence calculations, statistical comparisons, significant species, and abundance analysis |
| `Question4_results/` | Microbial co-occurrence network, correlation results, hub species, and sensitivity analysis |
| `Scripts/` | R scripts used to perform all four analyses |
| `README.md` | Project documentation, methodology, results, and execution instructions |

---

## 4. Software and Dependencies

### Programming Language

- R

### Required R Packages

The analyses use the following R packages:

```r
install.packages(c(
  "tidyverse",
  "vegan",
  "rstatix",
  "car",
  "ggpubr",
  "ggplot2",
  "dplyr",
  "tidyr",
  "igraph",
  "ggrepel",
  "reshape2",
  "pheatmap"
))
```

Load the required packages using:

```r
library(tidyverse)
library(vegan)
library(rstatix)
library(car)
library(ggpubr)
library(igraph)
library(ggrepel)
```

Individual scripts may require additional packages depending on their implementation.

---

# 5. Analysis and Methodology

## Question 1: Alpha Diversity Analysis

### Objective

To investigate within-sample microbial diversity and determine whether alpha diversity differs between Control and Diseased individuals and across different age categories.

### Methods

Two alpha diversity indices were calculated.

**1. Shannon Diversity Index**

The Shannon index measures microbial diversity by considering both species richness and relative abundance.

$$
H'=-\sum_{i=1}^{S}p_i\ln(p_i)
$$

Where:

- $p_i$ represents the relative abundance of species $i$.
- $S$ represents the number of observed species.

A higher Shannon index indicates greater microbial diversity.

**2. Pielou's Evenness Index**

Pielou's evenness measures how evenly microbial abundance is distributed across species.

$$
J=\frac{H'}{\ln(S)}
$$

A value closer to 1 indicates a more even microbial community.

### Statistical Analysis

| Analysis | Statistical Method |
|---|---|
| Control vs Diseased | Wilcoxon rank-sum test |
| Diversity across age categories | Kruskal–Wallis test |
| Pairwise age comparisons | Dunn's post-hoc test with BH correction |
| Age × Condition interaction | Two-way factorial ANOVA |
| Residual normality | Shapiro–Wilk test |
| Homogeneity of variance | Levene's test |

### Results

- Shannon diversity differed significantly between Control and Diseased samples.
- Control samples had a median Shannon diversity of **3.63**, compared with **3.44** in Diseased samples.
- Wilcoxon test p-value: **0.0001279**.
- No significant overall difference was observed across age categories (Kruskal–Wallis p = 0.37).
- The age × condition interaction was not significant at the 5% level (p = 0.07538).

The highest median Shannon diversity was observed in the Control–Senior group (3.64), while the lowest was observed in the Diseased–Middle age group (3.42).

### Visualizations

- Shannon diversity boxplots.
- Pielou's evenness boxplots.
- Age and study-condition interaction plot.

### Output Directory

`Question1_results/`

---

## Question 2: Beta Diversity Analysis

### Objective

To investigate differences in microbial community composition between samples and determine whether community composition is associated with study condition, age category, and gender.

### Methods

**1. Bray–Curtis Dissimilarity**

Bray–Curtis dissimilarity was calculated from the normalised species abundance profile (265 samples × 777 species) using `vegdist(method = "bray")`.

**2. Principal Coordinates Analysis (PCoA)**

PCoA was performed on the Bray–Curtis matrix using `cmdscale(k = 2, eig = TRUE)`. Variation explained was calculated as each eigenvalue divided by the sum of the positive eigenvalues (no negative-eigenvalue correction). Samples are coloured by study condition and shaped by age category.

**3. PERMANOVA**

Performed with `adonis2()` from the `vegan` package (9,999 free permutations, `set.seed(123)`), using the model `distance ~ study_condition + age_category + gender`.

- Primary analysis: marginal tests (`by = "margin"`), where each term is assessed after adjusting for the others.
- Sequential tests (`by = "terms"`) were also run.

**4. PERMDISP**

Multivariate dispersion between Control and Diseased samples was assessed with `betadisper()` (spatial median) and tested with `permutest()` (9,999 permutations).

**5. Alternative model specifications**

Sensitivity of the PERMANOVA results was checked with:

- Condition entered first (`study_condition + age_category + gender`)
- Age category entered first (`age_category + study_condition + gender`)
- Interaction model (`study_condition * age_category + gender`)

### PCoA Results

| Axis | Variation Explained |
|---|---:|
| PCo1 | 13.63% |
| PCo2 | 11.15% |
| Combined | 24.77% |

### PERMANOVA Results (marginal)

| Variable | Df | R² | F-value | p-value |
|---|---:|---:|---:|---:|
| Study condition | 1 | 0.01514 | 4.0865 | 0.0002 |
| Age category | 2 | 0.01724 | 2.3274 | 0.0001 |
| Gender | 1 | 0.00473 | 1.2760 | 0.1907 |

Study condition and age category were significantly associated with community composition; gender was not. Each significant variable explained under 2% of the total variation.

### Alternative Model Results

| Term | Condition first (R² / p) | Age first (R² / p) | With interaction (R² / p) |
|---|---:|---:|---:|
| Study condition | 0.01531 / 0.0001 | 0.01514 / 0.0001 | 0.01531 / 0.0001 |
| Age category | 0.01693 / 0.0004 | 0.01710 / 0.0003 | 0.01693 / 0.0004 |
| Gender | 0.00473 / 0.1901 | 0.00473 / 0.1901 | 0.00473 / 0.1907 |
| Condition × age | – | – | 0.00635 / 0.6878 |

Conclusions were unchanged by term order, and no condition × age interaction was detected.

### PERMDISP Results

- Mean distance to median: Control 0.391, Diseased 0.431
- F-value: **15.007**
- p-value: **0.0002**

Diseased samples were significantly more dispersed than Controls. The PERMANOVA result for study condition should therefore be interpreted cautiously, as differences in dispersion may contribute to it.

### Visualizations

- `pcoa_bray_curtis.png`: Bray–Curtis PCoA plot
- `multivariate_dispersion_study.png`: dispersion plot
- Variance explained (R²) bar chart, generated in the script

### Output Directory

`Question2_results/`

---

## Question 3: Differential Species Prevalence Analysis

### Objective

To identify microbial species that differ in prevalence between Control and Diseased samples and investigate the relationship between species prevalence and mean abundance.

### Methods

**1. Presence–Absence Transformation**

Each abundance value was converted to detected (1) or not detected (0):

$$
PA_{ij}=\begin{cases}1 & \text{if abundance}_{ij}>0\\ 0 & \text{otherwise}\end{cases}
$$

**2. Prevalence Calculation**

For each species, prevalence within a group was calculated as:

$$
\text{Prevalence}_{g}=\frac{\text{number of samples in group } g \text{ where the species is detected}}{\text{total samples in group } g}\times 100
$$

The prevalence difference was defined as Diseased − Control (in percentage points), so positive values indicate higher prevalence in Diseased samples.

**3. Prevalence Filtering**

Species detected in at least five of the 265 samples were retained for statistical testing.

- Initial species: **777**
- Retained species: **392**

**4. Statistical Testing**

For each retained species, a 2 × 2 table (Control/Diseased × present/absent) was tested with a two-sided Fisher's exact test. An odds ratio above 1 indicates higher odds of detection in Control samples.

P-values were adjusted using the Benjamini–Hochberg procedure, and species with adjusted p < 0.05 were considered significant.

**5. Correlation Analysis**

Mean abundance was calculated in two ways: across all samples in a group (zeros included) and only across samples where the species was detected. Spearman correlation was used to test whether prevalence differences track each of these mean abundance differences. Species not detected in one of the groups (no detected-only mean) were excluded, leaving 388 species for both correlations.

### Results

- Total species analyzed before filtering: **777**
- Species retained after prevalence filtering: **392**
- Statistically significant species after BH correction: **90**

The tables below list the largest differences among **significant** species. Among all 392 retained species (regardless of significance), the largest positive differences were *Lactobacillus_fermentum* (+14.49 pp), *Granulicatella_adiacens* (+14.05 pp) and *Scardovia_wiggsiae* (+12.42 pp); these three, plus the three negative species below, are shown in `Top_6_Species_Prevalence.png`, where significant species are marked with `*`. The negative top 3 are the same in both rankings.

### Significant Positive Prevalence Differences

| Species | Control Prevalence | Diseased Prevalence | Difference | Adjusted p-value |
|---|---:|---:|---:|---:|
| Lactobacillus_fermentum | 12.03% | 26.52% | +14.49 pp | 0.01772 |
| Streptococcus_suis | 0.00% | 6.06% | +6.06 pp | 0.01875 |
| Lactobacillus_casei | 0.00% | 5.30% | +5.30 pp | 0.03393 |

### Significant Negative Prevalence Differences

| Species | Control Prevalence | Diseased Prevalence | Difference | Adjusted p-value |
|---|---:|---:|---:|---:|
| Clostridium_aminobutyricum | 60.90% | 24.24% | -36.66 pp | 6.49 × 10⁻⁷ |
| Cardiobacterium_valvarum | 89.47% | 58.33% | -31.14 pp | 9.71 × 10⁻⁷ |
| Treponema_pedis | 63.16% | 33.33% | -29.82 pp | 8.04 × 10⁻⁵ |

### Correlation Results

| Comparison | Spearman rho | p-value |
|---|---:|---:|
| Prevalence difference vs mean abundance difference (zeros included) | 0.4677 | 1.77 × 10⁻²² |
| Prevalence difference vs mean abundance difference (detected samples only) | 0.1913 | 1.50 × 10⁻⁴ |

Prevalence differences were moderately and positively associated with mean abundance differences when zeros were included. The association was much weaker when abundance was considered only in samples where the species was detected, indicating that much of the overall abundance difference is driven by how often a species is present rather than how abundant it is when present.

### Visualizations

- `Top_6_Species_Prevalence.png`: Control vs Diseased prevalence for the top 3 positive and top 3 negative species.
- `Top_40_Prevalence_Differences.png`: 20 largest positive and 20 largest negative prevalence differences.
- `Prevalence_vs_Mean_Abundance.png`: prevalence difference vs mean abundance difference.
- - `Prevalence_Difference.png`: prevalence difference (pp) for the part (iii) species and the overall top 3 each way, coloured by FDR significance and labelled with pp and q-value.
- `Selected_Species_Comparison.png`: prevalence vs mean abundance for *Clostridium_aminobutyricum* (largest shift) and *Streptococcus_sinensis* (prevalence and abundance disagree).

### Output Directory

`Question3_Results/`

---

## Question 4: Microbial Co-occurrence Network Analysis

### Objective

To identify statistically significant associations between microbial species and construct a microbial co-occurrence network to investigate highly connected species.

### Methods

**1. Species Prevalence Filtering**

Species were filtered based on their prevalence across samples.

The primary analysis used a minimum prevalence threshold of 10%.

- Initial species: **777**
- Retained species: **267**

**2. Spearman Correlation**

Pairwise Spearman correlation coefficients were calculated between retained species.

Spearman correlation measures the strength and direction of monotonic associations between two variables.

**3. Multiple Testing Correction**

P-values from pairwise correlations were adjusted using the Benjamini–Hochberg procedure.

**4. Network Construction**

An association was included as an edge if it satisfied both conditions:

An association was included as an edge if it satisfied both of the following criteria:

- Absolute correlation coefficient (|rho|) >= 0.60
- Adjusted p-value < 0.05

An undirected network was constructed using the `igraph` package.

- Nodes represent microbial species.
- Edges represent statistically significant associations.
- Positive and negative associations represent the direction of correlation.

**5. Network Analysis**

The following network measures were calculated:

- Degree centrality
- Betweenness centrality
- Eigenvector centrality
- Network density
- Connected components
- Largest connected component

**6. Sensitivity Analysis**

Network robustness was investigated by varying:

- Species prevalence thresholds: 5%, 10%, and 15%.
- Absolute correlation thresholds: 0.50, 0.60, and 0.70.

Jaccard similarity was used to compare edge sets against the primary network.

### Primary Network Results

| Metric | Result |
|---|---:|
| Retained species | 267 |
| Species pairs tested | 35,511 |
| Significant associations | 117 |
| Positive associations | 117 |
| Negative associations | 0 |
| Network nodes | 67 |
| Network edges | 117 |
| Network density | 0.0529 |
| Connected components | 15 |
| Largest component size | 26 |
| Mean degree | 3.49 |

All selected associations in the primary network were positive.

### Top 10 Hub Species

| Rank | Species | Degree | Betweenness |
|---:|---|---:|---:|
| 1 | Treponema_medium | 19 | 0.05561 |
| 2 | Fretibacterium_fastidiosum | 14 | 0.01852 |
| 3 | Treponema_denticola | 13 | 0.01679 |
| 4 | Filifactor_alocis | 13 | 0.01531 |
| 5 | Prevotella_salivae | 8 | 0.00462 |
| 6 | Selenomonas_sputigena | 7 | 0.03186 |
| 7 | Treponema_socranskii | 7 | 0.01688 |
| 8 | Selenomonas_dianae | 7 | 0.01688 |
| 9 | Lachnoanaerobaculum_orale | 7 | 0.00408 |
| 10 | Veillonella_atypica | 7 | 0.00089 |

`Treponema_medium` had the highest degree centrality, with 19 connections, making it the most connected species in the primary network.

### Sensitivity Analysis

At the correlation threshold of 0.60:

- The 5% prevalence network had a Jaccard similarity of 1.00 against the primary network.
- The 10% prevalence network had a Jaccard similarity of 1.00.
- The 15% prevalence network had a Jaccard similarity of approximately 0.9915.

This indicates that the selected associations at the 0.60 correlation threshold were highly consistent across the tested prevalence thresholds.

The number of detected edges changed substantially when the correlation threshold was relaxed to 0.50 or increased to 0.70.

### Visualizations

- Microbial co-occurrence network.
- Top 10 hub species.
- Sensitivity analysis summary.

The network was also exported in GraphML format for further exploration in network visualization tools such as Gephi.

### Output Directory

`Question4_results/`

---

## 6. Overall Findings

The four analyses provide complementary insights into the microbial community.

| Analysis | Main Finding |
|---|---|
| Alpha diversity | Control samples showed significantly higher Shannon diversity than Diseased samples. |
| Beta diversity | Study condition and age category were significantly associated with community composition in PERMANOVA. |
| Species prevalence | 90 species showed statistically significant prevalence differences after BH correction. |
| Co-occurrence network | A network containing 67 nodes and 117 significant positive associations was identified. |

Together, these results indicate differences in microbial diversity, community composition, species prevalence, and microbial association patterns between the study groups.

---

## 7. Limitations

1. The analysis is based on the available sample size and the species included in the provided abundance matrix.
2. Statistical significance does not necessarily indicate biological or clinical importance.
3. The significant PERMDISP result requires cautious interpretation of the PERMANOVA findings.
4. The alpha diversity ANOVA residuals did not satisfy the normality assumption.
5. Correlation-based microbial networks identify statistical associations and do not establish direct biological interactions.
6. Microbial abundance data may be compositional, which can influence correlation-based association estimates.
7. Network hub species should not automatically be interpreted as keystone species.
8. The results describe associations and do not establish causal relationships.

---

## 8. How to Run the Project

### Step 1: Clone the Repository

```bash
git clone https://github.com/sunandini333/CoMeg_Assignment_2.git
```

### Step 2: Navigate to the Project Directory

```bash
cd CoMeg_Assignment_2
```

### Step 3: Install Required Packages

Open R or RStudio and install the required packages listed in Section 4.

### Step 4: Execute the Scripts

Run the scripts in the following order:

```r
source("Scripts/Question_1_Script.R")

source("Scripts/Question_2_Script.R")

source("Scripts/Question_3_Script.R")

source("Scripts/Question_4_Script.R")
```

### Step 5: Input Dataset Paths

The input datasets are stored inside the `Question/` directory.

Use the following relative paths when running the scripts from the repository root:

```r
metadata <- read.csv(
  "Question/metadata_all.xlsx.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

species_profile <- read.csv(
  "Question/species_profile_norm.xlsx.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
```

Ensure that the input file paths and output directory paths in the scripts match the repository structure.

Each script generates its corresponding result directory containing statistical outputs, tables, and visualizations.

---

## 9. Conclusion

This project presents a comprehensive statistical investigation of gut microbial communities using alpha diversity, beta diversity, differential prevalence, and co-occurrence network analysis.

The combination of ecological diversity indices, multivariate statistical testing, prevalence-based comparisons, and correlation-based network analysis provides a broader understanding of microbial community patterns across Control and Diseased individuals.

The project demonstrates the application of R-based statistical and computational approaches for microbiome data analysis.

---

## 10. Author

**CoMeg Assignment 2**

Computational analysis and statistical interpretation performed using R.
