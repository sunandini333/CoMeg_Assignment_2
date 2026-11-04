setwd("C:/SUNANDINI/Projects/Comeg_ass_2")
# ============================================================
# QUESTION 4: MICROBIAL CO-OCCURRENCE NETWORK ANALYSIS
#
# i.   Species filtering based on prevalence
# ii.  Pairwise Spearman correlation + BH correction
# iii. Co-occurrence network construction
# iv.  Hub species identification
# v.   Robustness / sensitivity analysis
# ============================================================


# ------------------------------------------------------------
# 1. INSTALL AND LOAD PACKAGES
# ------------------------------------------------------------

packages <- c(
  "tidyverse",
  "ggplot2",
  "igraph",
  "ggraph",
  "ggrepel",
  "scales"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(tidyverse)
library(ggplot2)
library(igraph)
library(ggraph)
library(ggrepel)
library(scales)


# ------------------------------------------------------------
# 2. CREATE OUTPUT DIRECTORY
# ------------------------------------------------------------

output_dir <- "Question4_Results"

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
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

# Rename first species-profile column
colnames(species_profile)[1] <- "sample_id"

metadata$sample_id <- as.character(metadata$sample_id)
species_profile$sample_id <- as.character(species_profile$sample_id)

# Check duplicate sample IDs
if (
  anyDuplicated(metadata$sample_id) ||
  anyDuplicated(species_profile$sample_id)
) {
  stop("Duplicate sample IDs detected.")
}

# Check missing sample IDs
if (
  anyNA(metadata$sample_id) ||
  anyNA(species_profile$sample_id) ||
  any(metadata$sample_id == "") ||
  any(species_profile$sample_id == "")
) {
  stop("Missing or empty sample IDs detected.")
}


# ------------------------------------------------------------
# 4. MATCH SAMPLES
# ------------------------------------------------------------

common_samples <- intersect(
  metadata$sample_id,
  species_profile$sample_id
)

if (length(common_samples) == 0) {
  stop("No common sample IDs found.")
}

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

cat("\nNUMBER OF MATCHED SAMPLES:\n")
print(nrow(metadata))


# ------------------------------------------------------------
# 5. PREPARE SPECIES ABUNDANCE MATRIX
# ------------------------------------------------------------

species_data <- species_profile %>%
  select(-sample_id)

# Convert all species columns to numeric
species_data <- species_data %>%
  mutate(
    across(
      everything(),
      ~ suppressWarnings(
        as.numeric(as.character(.x))
      )
    )
  )

species_matrix <- as.matrix(species_data)

rownames(species_matrix) <- species_profile$sample_id

# Validate abundance matrix
if (
  anyNA(species_matrix) ||
  any(!is.finite(species_matrix))
) {
  stop(
    "NA, NaN or Inf values detected after numeric conversion. ",
    "Check the original species abundance file."
  )
}

if (any(species_matrix < 0)) {
  stop("Negative abundance values detected.")
}

# Remove species absent in all samples
species_matrix <- species_matrix[
  ,
  colSums(species_matrix) > 0,
  drop = FALSE
]

# Remove samples with zero total abundance
valid_samples <- rownames(species_matrix)[
  rowSums(species_matrix) > 0
]

species_matrix <- species_matrix[
  valid_samples,
  ,
  drop = FALSE
]

# Align metadata
metadata <- metadata %>%
  filter(sample_id %in% valid_samples) %>%
  arrange(match(sample_id, valid_samples))

stopifnot(
  identical(
    metadata$sample_id,
    rownames(species_matrix)
  )
)

if (ncol(species_matrix) < 2) {
  stop("Fewer than two species remain.")
}

cat("\nINITIAL MATRIX DIMENSIONS:\n")
print(dim(species_matrix))


# ============================================================
# QUESTION 4(i): SPECIES FILTERING
# ============================================================


# ------------------------------------------------------------
# 6. CALCULATE SPECIES PREVALENCE
# ------------------------------------------------------------

# Prevalence = proportion of samples with non-zero abundance

species_prevalence <- colMeans(
  species_matrix > 0
)

prevalence_table <- data.frame(
  Species = colnames(species_matrix),
  
  Prevalence = as.numeric(species_prevalence),
  
  Prevalence_Percent = as.numeric(
    species_prevalence * 100
  ),
  
  Mean_Abundance = as.numeric(
    colMeans(species_matrix)
  ),
  
  Variance = as.numeric(
    apply(species_matrix, 2, var)
  )
)

write.csv(
  prevalence_table,
  file.path(output_dir, "Species_Prevalence.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 7. PRIMARY SPECIES FILTERING
# ------------------------------------------------------------

primary_prevalence_threshold <- 0.10

primary_species <- prevalence_table %>%
  filter(
    Prevalence >= primary_prevalence_threshold,
    Variance > 0
  ) %>%
  pull(Species)

if (length(primary_species) < 2) {
  stop("Fewer than two species passed the prevalence filter.")
}

filtered_matrix <- species_matrix[
  ,
  primary_species,
  drop = FALSE
]

cat("\nSPECIES FILTERING SUMMARY:\n")

cat(
  "Initial species:",
  ncol(species_matrix),
  "\n"
)

cat(
  "Species retained:",
  ncol(filtered_matrix),
  "\n"
)

cat(
  "Species removed:",
  ncol(species_matrix) - ncol(filtered_matrix),
  "\n"
)

cat(
  "Prevalence threshold:",
  primary_prevalence_threshold * 100,
  "%\n"
)

write.csv(
  data.frame(Species = colnames(filtered_matrix)),
  file.path(output_dir, "Retained_Species.csv"),
  row.names = FALSE
)


# ============================================================
# QUESTION 4(ii): SPEARMAN CORRELATION
# ============================================================


# ------------------------------------------------------------
# 8. PAIRWISE SPEARMAN CORRELATION FUNCTION
# ------------------------------------------------------------

calculate_spearman <- function(
    abundance_matrix,
    correlation_threshold = 0.60,
    adjusted_p_threshold = 0.05
) {
  
  abundance_matrix <- as.matrix(abundance_matrix)
  
  # Remove constant columns
  variable_columns <- apply(
    abundance_matrix,
    2,
    function(x) length(unique(x)) > 1
  )
  
  abundance_matrix <- abundance_matrix[
    ,
    variable_columns,
    drop = FALSE
  ]
  
  if (ncol(abundance_matrix) < 2) {
    stop("Fewer than two variable species remain.")
  }
  
  # Calculate Spearman correlation matrix
  correlation_matrix <- cor(
    abundance_matrix,
    method = "spearman",
    use = "pairwise.complete.obs"
  )
  
  # All unique species pairs
  species_pairs <- combn(
    colnames(abundance_matrix),
    2,
    simplify = FALSE
  )
  
  # Calculate pairwise correlations
  results <- lapply(
    species_pairs,
    function(pair) {
      
      species1 <- pair[1]
      species2 <- pair[2]
      
      x <- abundance_matrix[, species1]
      y <- abundance_matrix[, species2]
      
      complete <- complete.cases(x, y)
      
      x <- x[complete]
      y <- y[complete]
      
      if (
        length(x) < 3 ||
        length(unique(x)) < 2 ||
        length(unique(y)) < 2
      ) {
        
        return(
          data.frame(
            Species1 = species1,
            Species2 = species2,
            Spearman_Rho = NA_real_,
            P_Value = NA_real_,
            N_Samples = length(x)
          )
        )
      }
      
      test <- suppressWarnings(
        cor.test(
          x,
          y,
          method = "spearman",
          exact = FALSE
        )
      )
      
      data.frame(
        Species1 = species1,
        Species2 = species2,
        Spearman_Rho = unname(test$estimate),
        P_Value = test$p.value,
        N_Samples = length(x)
      )
    }
  )
  
  correlation_results <- bind_rows(results)
  
  # Remove invalid tests
  correlation_results <- correlation_results %>%
    filter(
      is.finite(Spearman_Rho),
      is.finite(P_Value)
    )
  
  if (nrow(correlation_results) == 0) {
    stop("No valid pairwise correlation tests were obtained.")
  }
  
  # BH correction
  correlation_results <- correlation_results %>%
    mutate(
      Adjusted_P_Value = p.adjust(
        P_Value,
        method = "BH"
      ),
      
      Correlation_Strength = case_when(
        abs(Spearman_Rho) >= 0.80 ~ "Very strong",
        abs(Spearman_Rho) >= 0.60 ~ "Strong",
        abs(Spearman_Rho) >= 0.40 ~ "Moderate",
        TRUE ~ "Weak"
      ),
      
      Association = case_when(
        Spearman_Rho > 0 ~ "Positive",
        Spearman_Rho < 0 ~ "Negative",
        TRUE ~ "Zero"
      ),
      
      Significant = (
        Adjusted_P_Value < adjusted_p_threshold
      )
    )
  
  # Select associations
  selected_edges <- correlation_results %>%
    filter(
      abs(Spearman_Rho) >= correlation_threshold,
      Adjusted_P_Value < adjusted_p_threshold
    )
  
  return(
    list(
      correlation_matrix = correlation_matrix,
      all_results = correlation_results,
      selected_edges = selected_edges
    )
  )
}


# ------------------------------------------------------------
# 9. RUN PRIMARY CORRELATION ANALYSIS
# ------------------------------------------------------------

primary_correlation_threshold <- 0.60
primary_adjusted_p_threshold <- 0.05

primary_results <- calculate_spearman(
  abundance_matrix = filtered_matrix,
  correlation_threshold = primary_correlation_threshold,
  adjusted_p_threshold = primary_adjusted_p_threshold
)

# Save correlation matrix
write.csv(
  primary_results$correlation_matrix,
  file.path(output_dir, "Spearman_Correlation_Matrix.csv")
)

# Save all correlations
write.csv(
  primary_results$all_results,
  file.path(output_dir, "All_Spearman_Correlations.csv"),
  row.names = FALSE
)

primary_edges <- primary_results$selected_edges

cat("\nPRIMARY CORRELATION RESULTS:\n")

cat(
  "Total species pairs tested:",
  nrow(primary_results$all_results),
  "\n"
)

cat(
  "Selected associations:",
  nrow(primary_edges),
  "\n"
)

cat(
  "Positive associations:",
  sum(primary_edges$Spearman_Rho > 0),
  "\n"
)

cat(
  "Negative associations:",
  sum(primary_edges$Spearman_Rho < 0),
  "\n"
)

write.csv(
  primary_edges,
  file.path(output_dir, "Selected_Species_Associations.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 10. SEPARATE POSITIVE AND NEGATIVE ASSOCIATIONS
# ------------------------------------------------------------

positive_edges <- primary_edges %>%
  filter(Spearman_Rho > 0) %>%
  arrange(desc(Spearman_Rho))

negative_edges <- primary_edges %>%
  filter(Spearman_Rho < 0) %>%
  arrange(Spearman_Rho)

write.csv(
  positive_edges,
  file.path(output_dir, "Positive_Associations.csv"),
  row.names = FALSE
)

write.csv(
  negative_edges,
  file.path(output_dir, "Negative_Associations.csv"),
  row.names = FALSE
)


# ============================================================
# QUESTION 4(iii): NETWORK CONSTRUCTION
# ============================================================


# ------------------------------------------------------------
# 11. NETWORK CONSTRUCTION FUNCTION
# ------------------------------------------------------------

construct_network <- function(edges) {
  
  if (is.null(edges) || nrow(edges) == 0) {
    return(NULL)
  }
  
  network_nodes <- unique(
    c(
      edges$Species1,
      edges$Species2
    )
  )
  
  node_data <- data.frame(
    name = network_nodes
  )
  
  edge_data <- edges %>%
    transmute(
      from = Species1,
      to = Species2,
      rho = Spearman_Rho,
      adjusted_p = Adjusted_P_Value,
      association = Association
    )
  
  microbial_network <- igraph::graph_from_data_frame(
    d = edge_data,
    directed = FALSE,
    vertices = node_data
  )
  
  # IMPORTANT:
  # Add degree as a vertex attribute.
  # This is required for ggraph's size = degree mapping.
  
  igraph::V(microbial_network)$degree <- igraph::degree(
    microbial_network,
    mode = "all"
  )
  
  return(microbial_network)
}


# ------------------------------------------------------------
# 12. CREATE PRIMARY NETWORK
# ------------------------------------------------------------

microbial_network <- construct_network(primary_edges)

if (is.null(microbial_network)) {
  
  warning(
    "No associations passed the selected thresholds. ",
    "The primary network is empty."
  )
  
} else {
  
  cat("\nNETWORK SUMMARY:\n")
  
  cat(
    "Number of nodes:",
    igraph::vcount(microbial_network),
    "\n"
  )
  
  cat(
    "Number of edges:",
    igraph::ecount(microbial_network),
    "\n"
  )
  
  cat(
    "Positive edges:",
    sum(igraph::E(microbial_network)$rho > 0),
    "\n"
  )
  
  cat(
    "Negative edges:",
    sum(igraph::E(microbial_network)$rho < 0),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # 13. NETWORK VISUALIZATION
  # ----------------------------------------------------------
  
  # Verify degree attribute before plotting
  stopifnot(
    "degree" %in% igraph::vertex_attr_names(microbial_network)
  )
  
  network_plot <- ggraph(
    microbial_network,
    layout = "fr"
  ) +
    
    geom_edge_link(
      aes(
        edge_colour = association,
        edge_width = abs(rho)
      ),
      alpha = 0.65
    ) +
    
    geom_node_point(
      aes(
        size = degree
      ),
      color = "steelblue",
      alpha = 0.85
    ) +
    
    geom_node_text(
      aes(
        label = name
      ),
      repel = TRUE,
      size = 3
    ) +
    
    scale_edge_colour_manual(
      values = c(
        "Positive" = "forestgreen",
        "Negative" = "firebrick"
      ),
      drop = FALSE
    ) +
    
    scale_edge_width(
      range = c(0.4, 2)
    ) +
    
    scale_size(
      range = c(3, 10)
    ) +
    
    labs(
      title = "Microbial Species Co-occurrence Network",
      subtitle = paste0(
        "Spearman |rho| >= ",
        primary_correlation_threshold,
        "; BH-adjusted p < ",
        primary_adjusted_p_threshold
      ),
      edge_colour = "Association",
      edge_width = "|Spearman rho|",
      size = "Node degree"
    ) +
    
    theme_void() +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 15
      ),
      legend.position = "right"
    )
  
  print(network_plot)
  
  ggsave(
    file.path(
      output_dir,
      "Microbial_Cooccurrence_Network.png"
    ),
    network_plot,
    width = 12,
    height = 10,
    dpi = 300
  )
  
  # Export network as GraphML
  igraph::write_graph(
    microbial_network,
    file.path(
      output_dir,
      "Microbial_Cooccurrence_Network.graphml"
    ),
    format = "graphml"
  )
}


# ============================================================
# QUESTION 4(iv): HUB SPECIES IDENTIFICATION
# ============================================================


# ------------------------------------------------------------
# 14. CALCULATE NETWORK METRICS
# ------------------------------------------------------------

if (!is.null(microbial_network)) {
  
  node_degree_table <- data.frame(
    
    Species = igraph::V(microbial_network)$name,
    
    Degree = igraph::degree(
      microbial_network,
      mode = "all"
    ),
    
    Betweenness = igraph::betweenness(
      microbial_network,
      directed = FALSE,
      normalized = TRUE
    ),
    
    Eigenvector_Centrality = igraph::eigen_centrality(
      microbial_network,
      directed = FALSE,
      scale = TRUE
    )$vector
    
  )
  
  # Calculate positive degree
  node_degree_table$Positive_Degree <- vapply(
    node_degree_table$Species,
    function(sp) {
      
      sum(
        (
          primary_edges$Species1 == sp |
            primary_edges$Species2 == sp
        ) &
          primary_edges$Spearman_Rho > 0
      )
    },
    numeric(1)
  )
  
  # Calculate negative degree
  node_degree_table$Negative_Degree <- vapply(
    node_degree_table$Species,
    function(sp) {
      
      sum(
        (
          primary_edges$Species1 == sp |
            primary_edges$Species2 == sp
        ) &
          primary_edges$Spearman_Rho < 0
      )
    },
    numeric(1)
  )
  
  node_degree_table <- node_degree_table %>%
    arrange(
      desc(Degree),
      desc(Betweenness)
    )
  
  write.csv(
    node_degree_table,
    file.path(
      output_dir,
      "Network_Node_Degree.csv"
    ),
    row.names = FALSE
  )
  
  
  # ----------------------------------------------------------
  # 15. IDENTIFY HUB SPECIES
  # ----------------------------------------------------------
  
  hub_species <- node_degree_table %>%
    slice_head(n = min(10, nrow(node_degree_table)))
  
  write.csv(
    hub_species,
    file.path(
      output_dir,
      "Hub_Species.csv"
    ),
    row.names = FALSE
  )
  
  cat("\nTOP HUB SPECIES:\n")
  print(hub_species)
  
  
  # Hub visualization
  hub_plot <- ggplot(
    hub_species,
    aes(
      x = reorder(Species, Degree),
      y = Degree
    )
  ) +
    
    geom_col(fill = "steelblue") +
    
    coord_flip() +
    
    labs(
      title = "Top 10 Hub Species",
      x = "Species",
      y = "Degree"
    ) +
    
    theme_minimal() +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 14
      )
    )
  
  print(hub_plot)
  
  ggsave(
    file.path(
      output_dir,
      "Top_10_Hub_Species.png"
    ),
    hub_plot,
    width = 9,
    height = 6,
    dpi = 300
  )
  
}


# ============================================================
# QUESTION 4(v): ROBUSTNESS / SENSITIVITY ANALYSIS
# ============================================================


# ------------------------------------------------------------
# 16. NETWORK SUMMARY FUNCTION
# ------------------------------------------------------------

network_summary <- function(edges) {
  
  network <- construct_network(edges)
  
  if (is.null(network)) {
    
    return(
      data.frame(
        Nodes = 0,
        Edges = 0,
        Density = 0,
        Connected_Components = 0,
        Largest_Component_Size = 0,
        Mean_Degree = 0
      )
    )
  }
  
  components_result <- igraph::components(network)
  
  data.frame(
    Nodes = igraph::vcount(network),
    
    Edges = igraph::ecount(network),
    
    Density = igraph::edge_density(
      network,
      loops = FALSE
    ),
    
    Connected_Components = components_result$no,
    
    Largest_Component_Size = max(
      components_result$csize
    ),
    
    Mean_Degree = mean(
      igraph::degree(network)
    )
  )
}


# ------------------------------------------------------------
# 17. SENSITIVITY ANALYSIS SETTINGS
# ------------------------------------------------------------

prevalence_thresholds <- c(
  0.05,
  0.10,
  0.15
)

correlation_thresholds <- c(
  0.50,
  0.60,
  0.70
)

sensitivity_results <- list()

sensitivity_counter <- 1


# ------------------------------------------------------------
# 18. RUN SENSITIVITY ANALYSIS
# ------------------------------------------------------------

# Create edge identifiers for comparison
make_edge_ids <- function(edges) {
  
  if (is.null(edges) || nrow(edges) == 0) {
    return(character(0))
  }
  
  apply(
    edges[, c("Species1", "Species2"), drop = FALSE],
    1,
    function(x) {
      paste(sort(x), collapse = " | ")
    }
  )
}

primary_edge_ids <- make_edge_ids(primary_edges)


for (prev_threshold in prevalence_thresholds) {
  
  cat(
    "\nRunning prevalence threshold:",
    prev_threshold * 100,
    "%\n"
  )
  
  # Select species
  selected_species <- prevalence_table %>%
    filter(
      Prevalence >= prev_threshold,
      Variance > 0
    ) %>%
    pull(Species)
  
  if (length(selected_species) < 2) {
    
    cat("Fewer than two species retained. Skipping.\n")
    
    next
  }
  
  current_matrix <- species_matrix[
    ,
    selected_species,
    drop = FALSE
  ]
  
  # Calculate pairwise correlations once for this filter
  current_correlation_results <- calculate_spearman(
    abundance_matrix = current_matrix,
    correlation_threshold = 0,
    adjusted_p_threshold = 0.05
  )
  
  for (rho_threshold in correlation_thresholds) {
    
    current_edges <- current_correlation_results$all_results %>%
      filter(
        abs(Spearman_Rho) >= rho_threshold,
        Adjusted_P_Value < 0.05
      )
    
    # Network summary
    current_network_summary <- network_summary(
      current_edges
    )
    
    # Jaccard similarity
    current_edge_ids <- make_edge_ids(current_edges)
    
    union_edges <- union(
      current_edge_ids,
      primary_edge_ids
    )
    
    intersection_edges <- intersect(
      current_edge_ids,
      primary_edge_ids
    )
    
    jaccard_similarity <- if (
      length(union_edges) == 0
    ) {
      NA_real_
    } else {
      length(intersection_edges) /
        length(union_edges)
    }
    
    sensitivity_results[[sensitivity_counter]] <-
      data.frame(
        
        Prevalence_Threshold = prev_threshold,
        
        Correlation_Threshold = rho_threshold,
        
        Species_Retained = length(selected_species),
        
        Significant_Edges = nrow(current_edges),
        
        Positive_Edges = sum(
          current_edges$Spearman_Rho > 0
        ),
        
        Negative_Edges = sum(
          current_edges$Spearman_Rho < 0
        ),
        
        Jaccard_Similarity_vs_Primary =
          jaccard_similarity
        
      ) %>%
      bind_cols(current_network_summary)
    
    sensitivity_counter <- sensitivity_counter + 1
  }
}


# ------------------------------------------------------------
# 19. SAVE SENSITIVITY RESULTS
# ------------------------------------------------------------

if (length(sensitivity_results) == 0) {
  
  sensitivity_table <- data.frame()
  
  warning("No sensitivity-analysis results were generated.")
  
} else {
  
  sensitivity_table <- bind_rows(
    sensitivity_results
  )
  
  write.csv(
    sensitivity_table,
    file.path(
      output_dir,
      "Sensitivity_Analysis_Summary.csv"
    ),
    row.names = FALSE
  )
  
  cat("\nSENSITIVITY ANALYSIS SUMMARY:\n")
  print(sensitivity_table)
  
}


# ------------------------------------------------------------
# 20. SENSITIVITY ANALYSIS VISUALIZATION
# ------------------------------------------------------------

if (nrow(sensitivity_table) > 0) {
  
  sensitivity_plot <- ggplot(
    sensitivity_table,
    aes(
      x = factor(Correlation_Threshold),
      y = Significant_Edges,
      fill = factor(Prevalence_Threshold)
    )
  ) +
    
    geom_col(
      position = "dodge"
    ) +
    
    labs(
      title = "Network Sensitivity Analysis",
      x = "Absolute Spearman Correlation Threshold",
      y = "Number of Significant Edges",
      fill = "Prevalence Threshold"
    ) +
    
    theme_minimal() +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 14
      )
    )
  
  print(sensitivity_plot)
  
  ggsave(
    file.path(
      output_dir,
      "Sensitivity_Analysis.png"
    ),
    sensitivity_plot,
    width = 9,
    height = 6,
    dpi = 300
  )
  
}


# ------------------------------------------------------------
# 21. NETWORK EDGE OVERLAP
# ------------------------------------------------------------

if (nrow(sensitivity_table) > 0) {
  
  edge_overlap_table <- sensitivity_table %>%
    select(
      Prevalence_Threshold,
      Correlation_Threshold,
      Significant_Edges,
      Jaccard_Similarity_vs_Primary
    )
  
  write.csv(
    edge_overlap_table,
    file.path(
      output_dir,
      "Network_Edge_Overlap.csv"
    ),
    row.names = FALSE
  )
  
}


# ============================================================
# 22. FINAL SUMMARY
# ============================================================

cat("\n============================================\n")
cat("QUESTION 4 ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("============================================\n")

cat(
  "\nMatched samples:",
  nrow(species_matrix),
  "\n"
)

cat(
  "Initial species:",
  ncol(species_matrix),
  "\n"
)

cat(
  "Species retained in primary analysis:",
  ncol(filtered_matrix),
  "\n"
)

cat(
  "Species pairs tested:",
  nrow(primary_results$all_results),
  "\n"
)

cat(
  "Selected primary associations:",
  nrow(primary_edges),
  "\n"
)

if (!is.null(microbial_network)) {
  
  cat(
    "Network nodes:",
    igraph::vcount(microbial_network),
    "\n"
  )
  
  cat(
    "Network edges:",
    igraph::ecount(microbial_network),
    "\n"
  )
  
}

cat(
  "\nResults saved in:",
  output_dir,
  "\n"
)

cat("\n============================================\n")
