# Install packages if needed:
required <- c("readxl", "dplyr", "tidyr", "ggplot2", "ggExtra", "hexbin")
newpkgs <- required[!(required %in% installed.packages()[, "Package"])]
if(length(newpkgs)) install.packages(newpkgs)

# Load libraries
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggExtra)
library(hexbin)

# Path to your uploaded Excel file (provided by you)
xlsx_path <- "/mnt/data/Summary results.xlsx"

# ---- 1) Read sheet 2 ----
sheet_index <- 2
df_raw <- readxl::read_excel(xlsx_path, sheet = sheet_index)
message("Loaded sheet ", sheet_index, " — columns found:")
print(colnames(df_raw))

# ---- 2) Auto-detect entropy & complexity column names ----
# We'll look for common patterns. We build a mapping of entropy -> column name,
# and complexity -> column name for each entropy type.

colnames_lower <- tolower(colnames(df_raw))

find_col <- function(patterns) {
  idx <- which(sapply(patterns, function(p) grepl(p, colnames_lower)))
  if (length(idx)==0) return(NA_character_)
  # prefer exact matches or the first matching column
  return(colnames(df_raw)[idx[1]])
}

mapping <- list(
  Shannon_H = find_col(c("^h_shannon$", "h_shannon", "shannon", "h_shannon", "h\\.shannon")),
  Shannon_C = find_col(c("^c_shannon$", "c_shannon", "c shannon", "c\\.shannon", "statistical_complexity", "c_shannon")),
  Renyi_H   = find_col(c("h_renyi", "renyi", "h\\.renyi")),
  Renyi_C   = find_col(c("c_renyi", "c renyi", "c\\.renyi")),
  Tsallis_H = find_col(c("h_tsallis", "tsallis", "h\\.tsallis")),
  Tsallis_C = find_col(c("c_tsallis", "c tsallis", "c\\.tsallis")),
  Fisher_H  = find_col(c("h_fisher", "fisher", "hfisher", "h\\.fisher")),
  Fisher_C  = find_col(c("c_fisher", "c fisher", "c\\.fisher"))
)

message("Detected column mapping (NA = not found):")
print(mapping)

# ---- 3) Detect model / class columns ----
# Common column names: "Model", "model", "Series", "SeriesType", "Type"
model_col <- find_col(c("^model$", "model", "series", "type"))
if (is.na(model_col)) {
  stop("Cannot find a model/series column. Please ensure your sheet contains a column named 'Model' or similar.")
}
message("Using model column: ", model_col)

# Detect group/class column (M1/M2/M3/M4). If not present, try to parse from Model.
group_col <- find_col(c("class", "group", "modelclass", "model_class", "group_label"))
if (is.na(group_col)) {
  # try to extract M1..M4 from model_col values
  df_tmp <- df_raw %>% mutate(.model_raw = as.character(.data[[model_col]]))
  df_tmp <- df_tmp %>% mutate(Group = ifelse(grepl("M[1234]", .model_raw, ignore.case = TRUE),
                                             toupper(sub(".*(M[1234]).*", "\\1", .model_raw, ignore.case = TRUE)),
                                             NA_character_))
  # if generated groups are mostly non-NA, use them
  if (sum(!is.na(df_tmp$Group)) / nrow(df_tmp) > 0.2) {
    df_raw$Group <- df_tmp$Group
    group_col <- "Group"
    message("Inferred group (M1..M4) from Model column; new column 'Group' created.")
  } else {
    # fallback: use Model itself as group
    df_raw$Group <- as.character(df_raw[[model_col]])
    group_col <- "Group"
    message("No separate group detected. Using Model values as group identifiers.")
  }
} else {
  message("Using provided group column: ", group_col)
}

# ---- 4) Detect coefficient columns to filter "all positive coefficients" ----
coef_candidates <- grep("coef|coeff|ar\\b|ma\\b|ar1|ar2|ma1|ma2", colnames_lower, value = TRUE)
coef_cols <- intersect(colnames(df_raw), coef_candidates) # names preserved
if (length(coef_cols) == 0) {
  message("No coefficient columns detected automatically. The script will proceed without filtering by positive coefficients.")
  has_coef <- FALSE
} else {
  message("Detected coefficient columns (for positivity check): ", paste(coef_cols, collapse = ", "))
  has_coef <- TRUE
}

# ---- 5) Prepare working dataframe ----
# keep only rows with at least one entropy/complexity pair detected
entropy_pairs <- list(
  Shannon = c(mapping$Shannon_H, mapping$Shannon_C),
  Renyi   = c(mapping$Renyi_H, mapping$Renyi_C),
  Tsallis = c(mapping$Tsallis_H, mapping$Tsallis_C),
  Fisher  = c(mapping$Fisher_H, mapping$Fisher_C)
)

# helper: valid pair?
valid_pairs <- lapply(entropy_pairs, function(x) all(!is.na(x)))
message("Valid entropy-complexity pairs found:")
print(valid_pairs)

# ---- 6) Filter dataset to "positive coefficients" if possible ----
df <- df_raw
if (has_coef) {
  # consider positive if ALL detected coefficient columns > 0 for a row
  df <- df %>% rowwise() %>%
    mutate(.allcoefpos = all(c_across(all_of(coef_cols)) > 0, na.rm = TRUE)) %>%
    ungroup()
  df_pos <- df %>% filter(.allcoefpos == TRUE)
  message("Rows with all positive coefficients: ", nrow(df_pos), " / ", nrow(df))
  df_pos$filter_reason <- "all_positive_coeffs"
} else {
  # fallback: do not filter (use full dataset) but mark it
  df$.allcoefpos <- NA
  df_pos <- df
  df_pos$filter_reason <- "no_coef_columns"
  message("Proceeding without coefficient filter. All rows kept.")
}

# ---- 7) Plot function ----
plot_entropy_complexity <- function(data, entropy_col, complexity_col,
                                    entropy_name = "Entropy", complexity_name = "Complexity",
                                    out_basename = "entropy_complexity",
                                    facet_by = group_col,
                                    point_alpha = 0.7, point_size = 1.8) {
  if (is.na(entropy_col) || is.na(complexity_col)) {
    message("Skipping ", entropy_name, " — columns not found.")
    return(NULL)
  }
  # Minimal tidy data
  plot_df <- data %>%
    dplyr::select(all_of(c(entropy_col, complexity_col, model_col, facet_by))) %>%
    rename(Entropy = !!entropy_col, Complexity = !!complexity_col, Model = !!model_col, Group = !!facet_by) %>%
    mutate(Group = as.factor(Group), Model = as.factor(Model))
  
  plot_df <- plot_df %>% filter(!is.na(Entropy) & !is.na(Complexity))
  
  if (nrow(plot_df) == 0) {
    message("No data available for plotting ", entropy_name)
    return(NULL)
  }
  
  p <- ggplot(plot_df, aes(x = Entropy, y = Complexity)) +
    # filled 2d density (contour-like)
    stat_density_2d(aes(fill = ..level..), geom = "polygon", contour = TRUE, alpha = 0.35) +
    geom_point(aes(color = Model, shape = Group), size = point_size, alpha = point_alpha) +
    geom_rug(sides = "bl", alpha = 0.2) +
    scale_fill_viridis_c(option = "magma") +
    theme_minimal(base_size = 12) +
    labs(title = paste0(entropy_name, " – Complexity plane"),
         subtitle = paste0("Data source: sheet ", sheet_index, " (positive coefficients subset if available)"),
         x = paste0(entropy_name, " (normalized)"),
         y = paste0(complexity_name),
         fill = "Density") +
    theme(legend.position = "bottom")
  
  # Facet by group if too many groups will reduce clarity; otherwise return single plot
  if (length(unique(plot_df$Group)) > 1 && length(unique(plot_df$Group)) <= 8) {
    p_f <- p + facet_wrap(~ Group, scales = "free")
  } else {
    p_f <- p
  }
  
  # Hexbin variant (useful for dense regions)
  p_hex <- ggplot(plot_df, aes(x = Entropy, y = Complexity)) +
    stat_binhex(bins = 40, alpha = 0.8) +
    geom_point(aes(color = Model, shape = Group), size = 0.9, alpha = 0.6) +
    theme_minimal() +
    labs(title = paste0(entropy_name, " – Complexity (hexbin)"),
         x = paste0(entropy_name), y = paste0(complexity_name))
  
  # Save plots
  ggplot2::ggsave(filename = paste0(out_basename, "_", entropy_name, "_points_density.png"),
                  plot = p_f, width = 10, height = 7, dpi = 300)
  ggplot2::ggsave(filename = paste0(out_basename, "_", entropy_name, "_hexbin.png"),
                  plot = p_hex, width = 10, height = 7, dpi = 300)
  
  message("Saved plots for ", entropy_name, " as PNGs with basename: ", out_basename)
  return(list(scatter = p_f, hexbin = p_hex, data = plot_df))
}

# ---- 8) Run analysis for each entropy type ----
results_list <- list()
for (entropy_name in names(entropy_pairs)) {
  pair <- entropy_pairs[[entropy_name]]
  Hcol <- pair[1]; Ccol <- pair[2]
  if (any(is.na(c(Hcol, Ccol)))) {
    message("Skipping ", entropy_name, ": required columns not found.")
    next
  }
  out <- plot_entropy_complexity(df_pos, entropy_col = Hcol, complexity_col = Ccol,
                                 entropy_name = entropy_name,
                                 complexity_name = paste0(entropy_name, "_Complexity"),
                                 out_basename = "EntropyComplexity")
  results_list[[entropy_name]] <- out
}

# ---- 9) Save a CSV summary of the subset used (positive coeffs or full) ----
summary_outfile <- "EntropyComplexity_Summary_used_rows.csv"
df_pos %>% dplyr::select(any_of(c(model_col, group_col, coef_cols, unlist(mapping)))) %>%
  write.csv(summary_outfile, row.names = FALSE)
message("Saved summary CSV: ", summary_outfile)

# ---- 10) Optional: create small summary statistics CSV per group & entropy ----
summary_stats <- list()
for (entropy_name in names(entropy_pairs)) {
  pair <- entropy_pairs[[entropy_name]]
  Hcol <- pair[1]; Ccol <- pair[2]
  if (any(is.na(c(Hcol, Ccol)))) next
  tmp <- df_pos %>%
    dplyr::filter(!is.na(.data[[Hcol]]) & !is.na(.data[[Ccol]])) %>%
    dplyr::group_by(!!sym(group_col)) %>%
    summarise(
      N = n(),
      Entropy_mean = mean(.data[[Hcol]], na.rm = TRUE),
      Entropy_sd   = sd(.data[[Hcol]], na.rm = TRUE),
      Complexity_mean = mean(.data[[Ccol]], na.rm = TRUE),
      Complexity_sd   = sd(.data[[Ccol]], na.rm = TRUE)
    ) %>% mutate(Entropy = entropy_name)
  summary_stats[[entropy_name]] <- tmp
}
all_stats <- bind_rows(summary_stats)
write.csv(all_stats, "EntropyComplexity_Group_Summaries.csv", row.names = FALSE)
message("Saved group summary stats -> EntropyComplexity_Group_Summaries.csv")

# ---- Done ----
message("All plotting done. Check PNGs in working directory and the summary CSVs.")
