
# --- Packages ---
# install.packages(c("readxl","dplyr","purrr","stringr","tibble","ggplot2","scales"))
library(readxl)
library(dplyr)
library(purrr)
library(stringr)
library(tibble)
library(ggplot2)
library(scales)

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
  # Prefer "<model>_n1000" (or n5000), fallback to "<model>" if needed
  sheet_candidates <- c(paste0(model, sheet_suffix), model)
  sheet_name <- sheet_candidates[sheet_candidates %in% available_sheets] %>% first()
  if (is.na(sheet_name)) {
    warning(sprintf("Sheet not found for model '%s' (tried: %s). Skipping.",
                    model, paste(sheet_candidates, collapse = ", ")))
    return(NULL)
  }
  
  df <- read_excel(path, sheet = sheet_name)
  # Normalize names
  names(df) <- names(df) %>% as.character() %>% str_trim()
  
  # Case-insensitive match for Rep and metric
  lower <- setNames(names(df), tolower(names(df)))
  rep_col    <- lower[["rep"]]
  metric_col <- lower[[tolower(metric)]]
  if (is.null(rep_col) || is.null(metric_col)) {
    warning(sprintf("Columns 'Rep' or '%s' not found on sheet '%s'. Skipping.", metric, sheet_name))
    return(NULL)
  }
  
  df %>%
    select(Rep = !!rep_col, value = !!metric_col) %>%
    mutate(Rep = suppressWarnings(as.numeric(Rep))) %>%
    filter(!is.na(Rep), !is.na(value)) %>%
    arrange(Rep) %>%
    rename(!!model := value)
}

# Build X matrix (Rep-aligned) for a case & metric
build_case_matrix <- function(path, models, metric, sheet_suffix, available_sheets) {
  pieces <- map(models, ~ read_model_metric(path, .x, metric, sheet_suffix, available_sheets))
  pieces <- compact(pieces)  # drop NULLs
  if (length(pieces) == 0) return(NULL)
  
  # Inner join on Rep to ensure common replicates across all models
  X <- reduce(pieces, ~ inner_join(.x, .y, by = "Rep"))
  # Keep only model columns (drop Rep), drop any NA rows
  X <- X %>% select(-Rep) %>% drop_na()
  if (ncol(X) < 2 || nrow(X) < 2) return(NULL)
  return(X)
}

# ERCA: eigenvalues & condition indices from correlation matrix
erca_from_matrix <- function(X) {
  # Correlation across models (columns)
  cor_mat <- cor(X, use = "pairwise.complete.obs")
  eig <- eigen(cor_mat, symmetric = TRUE)
  eigenvalues <- eig$values
  # Condition indices: sqrt(lambda_max / lambda_i)
  lambda_max <- max(eigenvalues, na.rm = TRUE)
  cond_index <- sqrt(lambda_max / eigenvalues)
  
  # Tables (sorted by descending eigenvalue for scree; also provide ascending for CI emphasis)
  tbl_desc <- tibble(
    component = paste0("Comp_", seq_along(eigenvalues)),
    eigenvalue = eigenvalues
  ) %>%
    arrange(desc(eigenvalue)) %>%
    mutate(
      proportion = eigenvalue / sum(eigenvalue),
      cumulative = cumsum(proportion),
      condition_index = sqrt(lambda_max / eigenvalue)
    )
  
  tbl_asc <- tbl_desc %>% arrange(eigenvalue)
  
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
    annotate("text", x = 1, y = 10, label = "10", vjust = -0.6, size = 3, color = "grey40") +
    annotate("text", x = 1, y = 30, label = "30", vjust = -0.6, size = 3, color = "grey40") +
    geom_text(aes(label = round(condition_index, 2)), vjust = -0.3, size = 3) +
    labs(title = title, x = "Component (sorted by eigenvalue ↑)", y = "Condition Index = sqrt(λmax / λi)") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file_path, p, width = 8, height = 5, dpi = 200)
}

# --- Run ERCA for all cases & metrics ---
sheets <- excel_sheets(data_path)
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
    scree_path <- file.path(out_dir, sprintf("ERCA_%s_%s_scree.png", case_name, metric))
    ci_path    <- file.path(out_dir, sprintf("ERCA_%s_%s_condition_indices.png", case_name, metric))
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

summary_df <- bind_rows(summary_rows)
summary_path <- file.path(out_dir, "ERCA_summary.csv")
write.csv(summary_df, summary_path, row.names = FALSE)

message("Done. Summary: ", summary_path)
print(summary_df)
