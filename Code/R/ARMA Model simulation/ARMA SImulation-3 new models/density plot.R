
# --- Packages ---
# install.packages(c("readxl","dplyr","purrr","stringr","tibble","ggplot2","scales","tidyr"))
library(readxl)
library(dplyr)
library(purrr)
library(stringr)
library(tibble)
library(ggplot2)
library(scales)
library(tidyr)    # <-- Needed for drop_na()

# --- Inputs ---
data_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/ARMA Time series results n1000_n5000/HC_Results_all_Models_n1000_n5000.xlsx"
sheet_suffix <- "_n1000"   # change to "_n5000" if you want the n=5000 sheets
metrics <- c("H_Shannon","C_Shannon")

cases <- list(
  Case1 = c("AR1_M1","AR1_M2","AR2_M1",
            "MA1_M1","MA1_M2","MA2_M1",
            "ARMA11_M1","ARMA11_M2","ARMA22_M1"),
  Case2 = c("AR1_M3","AR1_M4","AR2_M4",
            "MA1_M3","MA1_M4","MA2_M4",
            "ARMA11_M3","ARMA11_M4","ARMA22_M4"),
  Case3 = c("AR2_M2","AR2_M3",
            "MA2_M2","MA2_M3",
            "ARMA22_M2","ARMA22_M3")
)

# --- Output folder ---
out_dir <- file.path(getwd(), "ERCA_Output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- Helpers ---

# Read one model's metric column from its sheet, align on Rep
read_model_metric <- function(path, model, metric, sheet_suffix, available_sheets) {
  sheet_candidates <- c(paste0(model, sheet_suffix), model)
  sheet_name <- dplyr::first(sheet_candidates[sheet_candidates %in% available_sheets])
  if (is.na(sheet_name)) {
    warning(sprintf("Sheet not found for model '%s' (tried: %s). Skipping.",
                    model, paste(sheet_candidates, collapse = ", ")))
    return(NULL)
  }
  
  df <- readxl::read_excel(path, sheet = sheet_name)
  names(df) <- names(df) |> as.character() |> stringr::str_trim()
  
  # Case-insensitive match for Rep and metric
  nm_lower <- setNames(names(df), tolower(names(df)))
  rep_col    <- nm_lower[["rep"]]
  metric_col <- nm_lower[[tolower(metric)]]
  if (is.null(rep_col) || is.null(metric_col)) {
    warning(sprintf("Columns 'Rep' or '%s' not found on sheet '%s'. Skipping.", metric, sheet_name))
    return(NULL)
  }
  
  df |>
    dplyr::select(Rep = !!rep_col, value = !!metric_col) |>
    dplyr::mutate(Rep = suppressWarnings(as.numeric(Rep))) |>
    tidyr::drop_na(Rep, value) |>
    dplyr::arrange(Rep) |>
    dplyr::rename(!!model := value)
}

# Build X matrix (Rep-aligned) for a case & metric
build_case_matrix <- function(path, models, metric, sheet_suffix, available_sheets) {
  pieces <- purrr::map(models, ~ read_model_metric(path, .x, metric, sheet_suffix, available_sheets)) |> purrr::compact()
  if (length(pieces) == 0) return(NULL)
  
  # Inner join on Rep to ensure common replicates across all models
  X <- purrr::reduce(pieces, ~ dplyr::inner_join(.x, .y, by = "Rep"))
  X <- X |> dplyr::select(-Rep) |> tidyr::drop_na()
  if (ncol(X) < 2 || nrow(X) < 2) return(NULL)
  return(X)
}

# ERCA: eigenvalues & condition indices from correlation matrix
erca_from_matrix <- function(X) {
  cor_mat <- stats::cor(X, use = "pairwise.complete.obs")
  eig <- eigen(cor_mat, symmetric = TRUE)
  eigenvalues <- eig$values
  lambda_max <- max(eigenvalues, na.rm = TRUE)
  cond_index <- sqrt(lambda_max / eigenvalues)
  
  tbl_desc <- tibble(
    component = paste0("Comp_", seq_along(eigenvalues)),
    eigenvalue = eigenvalues
  ) |>
    dplyr::arrange(dplyr::desc(eigenvalue)) |>
    dplyr::mutate(
      proportion = eigenvalue / sum(eigenvalue),
      cumulative = cumsum(proportion),
      condition_index = sqrt(lambda_max / eigenvalue)
    )
  
  tbl_asc <- tbl_desc |> dplyr::arrange(eigenvalue)
  
  list(
    cor = cor_mat,
    eigenvalues = eigenvalues,
    cond_index = cond_index,
    table_desc = tbl_desc,
    table_asc  = tbl_asc,
    condition_number = sqrt(lambda_max / min(eigenvalues))
  )
}

# Plot scree and condition-index bars
plot_scree <- function(tbl_desc, title, file_path) {
  p <- ggplot(tbl_desc, aes(x = reorder(component, -eigenvalue), y = eigenvalue)) +
    geom_col(fill = "#2F7ED8") +
    geom_text(aes(label = round(eigenvalue, 2)), vjust = -0.3, size = 3) +
    labs(title = title, x = "Component (sorted by eigenvalue ↓)", y = "Eigenvalue") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file_path, p, width = 8, height = 5, dpi = 200)
}

plot_condition_indices <- function(tbl_asc, title, file_path) {
  p <- ggplot(tbl_asc, aes(x = reorder(component, eigenvalue), y = condition_index)) +
    geom_col(fill = "#D64F2F") +
    geom_hline(yintercept = 10, linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = 30, linetype = "dashed", color = "grey40") +
    geom_text(aes(label = round(condition_index, 2)), vjust = -0.3, size = 3) +
    labs(title = title, x = "Component (sorted by eigenvalue ↑)", y = "Condition Index = sqrt(λmax / λi)") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file_path, p, width = 8, height = 5, dpi = 200)
}

# --- Run ERCA for all cases & metrics ---
sheets <- readxl::excel_sheets(data_path)
summary_rows <- list()

for (case_name in names(cases)) {
  models <- cases[[case_name]]
  
  for (metric in metrics) {
    X <- build_case_matrix(data_path, models, metric, sheet_suffix, sheets)
    if (is.null(X)) {
      message(sprintf("Skipping %s / %s (matrix empty or not enough data).", case_name, metric))
      next
    }
    
    er <- erca_from_matrix(X)
    
    # Save tables
    tbl_path_desc <- file.path(out_dir, sprintf("ERCA_%s_%s_table_desc.csv", case_name, metric))
    tbl_path_asc  <- file.path(out_dir, sprintf("ERCA_%s_%s_table_asc.csv",  case_name, metric))
    write.csv(er$table_desc, tbl_path_desc, row.names = FALSE)
    write.csv(er$table_asc,  tbl_path_asc,  row.names = FALSE)
    
    # Save plots
    scree_path <- file.path(out_dir, sprintf("ERCA_%s_%s_scree.pdf", case_name, metric))
    ci_path    <- file.path(out_dir, sprintf("ERCA_%s_%s_condition_indices.pdf", case_name, metric))
    plot_scree(er$table_desc, sprintf("Scree Plot — %s (%s)", case_name, metric), scree_path)
    plot_condition_indices(er$table_asc, sprintf("Condition Indices — %s (%s)", case_name, metric), ci_path)
    
    # Append summary
    summary_rows[[length(summary_rows)+1]] <- tibble(
      case = case_name,
      metric = metric,
      p = ncol(X),
      n_rows = nrow(X),
      condition_number = er$condition_number,
      max_condition_index = max(er$cond_index),
      min_eigenvalue = min(er$eigenvalues),
      max_eigenvalue = max(er$eigenvalues),
      any_ci_ge_10 = max(er$cond_index) >= 10,
      any_ci_ge_30 = max(er$cond_index) >= 30,
      table_desc_csv = tbl_path_desc,
      table_asc_csv  = tbl_path_asc,
      scree_png = scree_path,
      ci_png    = ci_path
    )
  }
}

summary_df <- dplyr::bind_rows(summary_rows)
summary_path <- file.path(out_dir, "ERCA_summary.csv")
write.csv(summary_df, summary_path, row.names = FALSE)

message("Done. Summary: ", summary_path)
print(summary_df)


#-------------------------------------------------------------------------------------------
# This script performs PCA and clustering on time series metrics from an Excel file.
# It processes multiple cases, computes PCA scores, identifies clusters, and generates plots.
# Outputs include PCA scores, model centroids, cluster assignments, and various visualizations.
#-------------------------------------------------------------------------------------------
library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(tibble)

#--------------------------------------------------------
# Paths
#--------------------------------------------------------
data_path  <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/ARMA Time series results n1000_n5000/HC_Results_all_Models_n1000_n5000.xlsx"
base_out   <- "C:/Users/UserA1/Documents/ERCA_Output"

# Make base output dir
dir.create(base_out, recursive = TRUE, showWarnings = FALSE)

#--------------------------------------------------------
# Case definitions (models per case)
#--------------------------------------------------------
cases <- list(
  Case1 = c("AR1_M1","AR1_M2","AR2_M1",
            "MA1_M1","MA1_M2","MA2_M1",
            "ARMA11_M1","ARMA11_M2","ARMA22_M1"),
  Case2 = c("AR1_M3","AR1_M4","AR2_M4",
            "MA1_M3","MA1_M4","MA2_M4",
            "ARMA11_M3","ARMA11_M4","ARMA22_M4"),
  Case3 = c("AR2_M2","AR2_M3",
            "MA2_M2","MA2_M3",
            "ARMA22_M2","ARMA22_M3")
)

#--------------------------------------------------------
# Sample sizes
#--------------------------------------------------------
sample_sizes <- c(1000, 5000)

#--------------------------------------------------------
# Helper: choose metric columns for PCA
#   - Prefer H_Shannon, C_Shannon if present
#   - Otherwise use all numeric columns except 'n'
#--------------------------------------------------------
select_metric_columns <- function(df) {
  candidates <- c("H_Shannon", "C_Shannon", "H", "C")
  present    <- intersect(candidates, names(df))
  
  if (length(present) >= 2) {
    message("  Using metric columns: ", paste(present, collapse = ", "))
    return(present)
  } else {
    # fallback: all numeric except 'n'
    num_cols <- names(df)[sapply(df, is.numeric)]
    num_cols <- setdiff(num_cols, "n")
    message("  Using all numeric metric columns (fallback): ",
            paste(num_cols, collapse = ", "))
    return(num_cols)
  }
}

#--------------------------------------------------------
# Main loops: over cases and sample sizes
#--------------------------------------------------------
set.seed(123)  # reproducible clustering

for (case_name in names(cases)) {
  
  models_case <- cases[[case_name]]
  
  for (n_val in sample_sizes) {
    
    message("===== Processing ", case_name, ", n = ", n_val, " =====")
    
    #--------------------------------------------------------
    # Read Excel sheets for this case & n
    #--------------------------------------------------------
    df_list <- lapply(models_case, function(m) {
      sheet_name <- paste0(m, "_n", n_val)
      message("  Reading sheet: ", sheet_name)
      read_excel(data_path, sheet = sheet_name) |>
        mutate(Model = m, n = n_val)
    })
    
    df_case <- bind_rows(df_list)
    
    # Decide on metric columns for PCA
    metric_cols <- select_metric_columns(df_case)
    
    # Drop rows with non-finite values in any metric
    df_case <- df_case |>
      filter(across(all_of(metric_cols), ~ is.finite(.x)))
    
    if (nrow(df_case) == 0) {
      warning("No data for ", case_name, " n = ", n_val, ", skipping.")
      next
    }
    
    #--------------------------------------------------------
    # Classify Families (AR / MA / ARMA) and Type (Model)
    #--------------------------------------------------------
    df_case <- df_case |>
      mutate(
        Family = case_when(
          str_starts(Model, "ARMA") ~ "ARMA",
          str_starts(Model, "AR")   ~ "AR",
          str_starts(Model, "MA")   ~ "MA",
          TRUE                      ~ "Other"
        ),
        Type = factor(Model, levels = models_case)
      )
    
    #--------------------------------------------------------
    # PCA on metrics
    #--------------------------------------------------------
    X <- as.matrix(df_case[, metric_cols])
    X_scaled <- scale(X)   # centre & scale each metric
    
    pca <- prcomp(X_scaled, center = FALSE, scale. = FALSE)
    
    # PCA scores (PC1, PC2, ...)
    scores <- as.data.frame(pca$x)
    colnames(scores) <- paste0("PC", seq_len(ncol(scores)))
    
    # Keep at least PC1, PC2
    # Combine with labels
    pca_scores <- df_case |>
      select(Model, Family, Type, n) |>
      bind_cols(scores)
    
    #--------------------------------------------------------
    # Variance explained (for scree plot)
    #--------------------------------------------------------
    var_explained <- pca$sdev^2
    var_explained <- var_explained / sum(var_explained)
    scree_df <- tibble(
      PC = factor(paste0("PC", seq_along(var_explained)),
                  levels = paste0("PC", seq_along(var_explained))),
      VarExplained = var_explained
    )
    
    #--------------------------------------------------------
    # Clustering in PCA space (k-means on first 2 PCs)
    #   You can change k_clusters as you like.
    #   Here we use number of distinct models in this case.
    #--------------------------------------------------------
    k_clusters <- length(unique(pca_scores$Type))
    message("  Clustering with k = ", k_clusters, " (number of models in this case)")
    
    km <- kmeans(pca_scores[, c("PC1", "PC2")], centers = k_clusters, nstart = 25)
    
    pca_scores$Cluster <- factor(km$cluster)
    
    # Cluster assignment per model Type (majority vote)
    cluster_by_type <- pca_scores |>
      group_by(Type) |>
      count(Cluster, name = "n_in_cluster") |>
      slice_max(order_by = n_in_cluster, n = 1, with_ties = FALSE) |>
      ungroup()
    
    #--------------------------------------------------------
    # Model centroids in PCA space
    #--------------------------------------------------------
    centroids <- pca_scores |>
      group_by(Type, Model, Family) |>
      summarise(
        PC1_mean = mean(PC1),
        PC2_mean = mean(PC2),
        .groups = "drop"
      )
    
    #--------------------------------------------------------
    # Prepare output directory for this Case & n
    #--------------------------------------------------------
    case_out_dir <- file.path(base_out, case_name, paste0("n", n_val))
    dir.create(case_out_dir, recursive = TRUE, showWarnings = FALSE)
    
    #--------------------------------------------------------
    # Save numeric outputs as CSV
    #--------------------------------------------------------
    # 1) PCA scores per observation
    write.csv(
      pca_scores,
      file = file.path(case_out_dir, paste0("PCA_Scores_", case_name, "_n", n_val, ".csv")),
      row.names = FALSE
    )
    
    # 2) Model centroids in PCA space
    write.csv(
      centroids,
      file = file.path(case_out_dir, paste0("PCA_Model_Centroids_", case_name, "_n", n_val, ".csv")),
      row.names = FALSE
    )
    
    # 3) Cluster assignment per model
    write.csv(
      cluster_by_type,
      file = file.path(case_out_dir, paste0("Cluster_Assignments_", case_name, "_n", n_val, ".csv")),
      row.names = FALSE
    )
    
    message("  Saved PCA and clustering CSVs to: ", case_out_dir)
    
    #--------------------------------------------------------
    # PLOTS
    #--------------------------------------------------------
    
    # 1) Scree plot
    p_scree <- ggplot(scree_df, aes(x = PC, y = VarExplained)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = scales::percent(VarExplained, accuracy = 0.1)),
                vjust = -0.3, size = 3) +
      labs(
        title = paste("Scree plot —", case_name, "(n =", n_val, ")"),
        x = "Principal Component",
        y = "Proportion of variance explained"
      ) +
      theme_minimal(base_size = 14)
    
    ggsave(
      filename = file.path(case_out_dir, paste0("Scree_", case_name, "_n", n_val, ".pdf")),
      plot = p_scree,
      width = 7, height = 5
    )
    
    # 2) PC1 vs PC2 coloured by Model Type
    p_pca_type <- ggplot(pca_scores, aes(PC1, PC2, colour = Type, shape = Family)) +
      geom_point(alpha = 0.5, size = 2) +
      geom_point(data = centroids, aes(PC1_mean, PC2_mean, colour = Type),
                 size = 4, shape = 4, stroke = 1.2, inherit.aes = FALSE) +
      labs(
        title = paste("PCA (PC1 vs PC2) by Model Type —", case_name, "(n =", n_val, ")"),
        x = "PC1",
        y = "PC2"
      ) +
      theme_minimal(base_size = 14)
    
    ggsave(
      filename = file.path(case_out_dir, paste0("PCA_PC1_PC2_by_Type_", case_name, "_n", n_val, ".pdf")),
      plot = p_pca_type,
      width = 7, height = 5
    )
    
    # 3) PC1 vs PC2 coloured by cluster
    p_pca_cluster <- ggplot(pca_scores, aes(PC1, PC2, colour = Cluster)) +
      geom_point(alpha = 0.5, size = 2) +
      labs(
        title = paste("PCA (PC1 vs PC2) by k-means Cluster —", case_name, "(n =", n_val, ")"),
        x = "PC1",
        y = "PC2"
      ) +
      theme_minimal(base_size = 14)
    
    ggsave(
      filename = file.path(case_out_dir, paste0("PCA_PC1_PC2_by_Cluster_", case_name, "_n", n_val, ".pdf")),
      plot = p_pca_cluster,
      width = 7, height = 5
    )
    
    message("  Saved plots to: ", case_out_dir, "\n")
  }
}

message("🎉 PCA + clustering completed for all cases and sample sizes!")
