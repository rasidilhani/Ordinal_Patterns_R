#--------------------------------------------------------
# Libraries
#--------------------------------------------------------
library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)
library(StatOrdPattHxC)

#--------------------------------------------------------
# Paths
#--------------------------------------------------------
#data_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/ARMA Time series results/Entropy_Complexity_Results_all_Models.xlsx"
data_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/ARMA Time series results n5000_n10000/HC_Results_all_Models_n5000_n10000.xlsx"
output_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Plots/ARMA plots/ARMA Plots n5000_n10000"
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
# Choose Case and sample size
#--------------------------------------------------------
selected_case <- "Case1"
selected_n    <- 5000

models_case <- cases[[selected_case]]

#--------------------------------------------------------
# Read only sheets for Case1, n = 5000
#--------------------------------------------------------
df_list <- lapply(models_case, function(m) {
  sheet_name <- paste0(m, "_n", selected_n)
  message("Reading sheet: ", sheet_name)
  read_excel(data_path, sheet = sheet_name) |>
    mutate(Model = m, n = selected_n)
})

df_case <- bind_rows(df_list) |>
  rename(
    H = H_Shannon,
    C = C_Shannon
  ) |>
  filter(is.finite(H), is.finite(C))

#--------------------------------------------------------
# Classify model type: AR / MA / ARMA
####-----------------------------------------------
df_case <- df_case |>
  mutate(Type = case_when(
    str_starts(Model, "ARMA") ~ "ARMA",
    str_starts(Model, "AR")   ~ "AR",
    str_starts(Model, "MA")   ~ "MA",
    TRUE                      ~ "Other"
  ))

df_case$Type <- factor(df_case$Type, levels = c("AR","MA","ARMA","Other"))

##---------------------------------------------------------------
# Set plotting limits
##-------------------------------------------------
H_min <- min(df_case$H) - 0.02
H_max <- max(df_case$H) + 0.02
C_min <- min(df_case$C) - 0.02
C_max <- max(df_case$C) + 0.02

#-----------------------------------------------------
# Feasible region
##---------------------------------------------------------
data("LinfLsup")
Linf <- subset(LinfLsup, Side == "Lower" & Dimension == "3")
Lsup <- subset(LinfLsup, Side == "Upper" & Dimension == "3")

Linf_crop <- Linf |> filter(H >= H_min, H <= H_max)
Lsup_crop <- Lsup |> filter(H >= H_min, H <= H_max)

##------------------------------------------------------
# Improved Colors (viridis-inspired)
##----------------------------------------------------------
type_cols <- c(
  "AR"   = "#21918c",
  "MA"   = "#3182bd",
  "ARMA" = "#b2182b",
  "Other" = "grey50"
)

##--------------------------------------------------------------
# Heatmap with Top 3 Improvements
#---------------------------------------------------------------------
p <- ggplot(df_case, aes(H, C)) +
  
  # 1. Feasible region boundaries
  geom_line(data = Linf_crop, aes(H, C),
            colour = "black", linewidth = 0.8, inherit.aes = FALSE) +
  geom_line(data = Lsup_crop, aes(H, C),
            colour = "black", linewidth = 0.8, inherit.aes = FALSE) +
  
  # Background points
  geom_point(alpha = 0.08, size = 0.8, colour = "black") +
  
  ###------------------------------------------------------------
# 2. High-resolution density polygons (bins = 50)
##----------------------------------------------------------------------
stat_density_2d(
  aes(fill = Type, alpha = after_stat(level)),
  geom = "polygon",
  contour = TRUE,
  bins = 15,
  colour = NA
) +
  
  ##--------------------------------------------------------------
# 3. Add contour lines on top (improvement #1)
#-----------------------------------------------------------------
stat_density_2d(
  aes(color = Type),
  contour = TRUE,
  bins = 10,
  linewidth = 0.05
) +
  
  # Apply color scales
  scale_fill_manual(values = type_cols, name = "Model type") +
  scale_color_manual(values = type_cols, guide = "none") +
  
  # Opacity improvement (Option 1)
  scale_alpha(range = c(0.10, 0.55), guide = "none") +
  
  labs(
    title = paste("Heatmap —", selected_case, "(n =", selected_n, ")"),
    x = expression(italic(H)),
    y = expression(italic(C))
  ) +
  
  coord_cartesian(xlim = c(H_min, H_max), ylim = c(C_min, C_max)) +
  
  theme_minimal(base_size = 16) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid       = element_blank(),
    legend.position  = "right"
  )

#----------------------------------------------------------
# Print final plot
#-----------------------------------------------------------
print(p)



#----------------------------------------------------------# Save plot to file
#-----------------------------------------------------------
output_plot_path <- file.path(
  "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Plots/ARMA plots/ARMA Plots n5000_n10000",
  selected_case, "Plots",
  paste0("HC_Heatmap_", selected_case, "_n", selected_n, ".pdf")
)
ggsave(output_plot_path, p, width = 10, height = 8)


#-------------------------------------------------------------------------
# End of example_code_heatmap.R for Case 1
#
#-------------------------------------------------------------------------
# In this script we generate a heatmap of the Entropy-Complexity plane
# for a selected case (Case2, or Case3) and sample size (n = 5000 or n = 10000).
# The heatmap visualizes the density of points corresponding to different
# time series models (AR, MA, ARMA) in the Entropy-Complexity space.
#--------------------------------------------------------
# Libraries
#--------------------------------------------------------
library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)
library(StatOrdPattHxC)

#--------------------------------------------------------
# Paths
#--------------------------------------------------------
data_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/ARMA Time series results n5000_n10000/HC_Results_all_Models_n5000_n10000.xlsx"

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
# Choose Case and sample size
#--------------------------------------------------------
selected_case <- "Case2"
selected_n    <- 5000

models_case <- cases[[selected_case]]

#--------------------------------------------------------
# Read Excel sheets
#--------------------------------------------------------
df_list <- lapply(models_case, function(m) {
  sheet_name <- paste0(m, "_n", selected_n)
  read_excel(data_path, sheet = sheet_name) |>
    mutate(Model = m)
})

df_case <- bind_rows(df_list) |>
  rename(H = H_Shannon, C = C_Shannon) |>
  filter(is.finite(H), is.finite(C))

#--------------------------------------------------------
# Classify model type
#--------------------------------------------------------
df_case <- df_case |>
  mutate(Type = case_when(
    str_starts(Model, "ARMA") ~ "ARMA",
    str_starts(Model, "AR")   ~ "AR",
    str_starts(Model, "MA")   ~ "MA",
    TRUE                      ~ "Other"
  ))

df_case$Type <- factor(df_case$Type, levels = c("AR","MA","ARMA","Other"))

#--------------------------------------------------------
# Feasible region (crop to H-range)
#--------------------------------------------------------
data("LinfLsup")
Linf <- subset(LinfLsup, Side == "Lower" & Dimension == 3)
Lsup <- subset(LinfLsup, Side == "Upper" & Dimension == 3)

H_min <- min(df_case$H)
H_max <- max(df_case$H)

Linf_crop <- Linf |> filter(H >= H_min, H <= H_max)
Lsup_crop <- Lsup |> filter(H >= H_min, H <= H_max)

# Boundary polygon (CROPPED)
boundary_df <- rbind(
  Linf_crop[, c("H","C")],
  Lsup_crop[order(Lsup_crop$H, decreasing = TRUE), c("H","C")],
  Linf_crop[1, c("H","C")]
)
boundary_df <- as.data.frame(boundary_df)

#--------------------------------------------------------
# Dynamic bandwidth for density smoothing
#--------------------------------------------------------
Hx <- diff(range(df_case$H))
Cx <- diff(range(df_case$C))

h_vec <- c(Hx / 5, Cx / 5)   # adaptive smoothing

#--------------------------------------------------------
# Colors
#--------------------------------------------------------
type_cols <- c(
  "AR"   = "#21918c",
  "MA"   = "#3182bd",
  "ARMA" = "#b2182b",
  "Other" = "grey50"
)

#--------------------------------------------------------
# Plot
#--------------------------------------------------------
p <- ggplot(df_case, aes(H, C)) +
  
  # Feasible region boundaries
  geom_line(data = Linf_crop, aes(H, C),
            colour = "black", linewidth = 0.8) +
  geom_line(data = Lsup_crop, aes(H, C),
            colour = "black", linewidth = 0.8) +
  
  # Background points
  geom_point(alpha = 0.07, size = 0.8, colour = "black") +
  
  # Density polygons (dynamic smoothing!)
  stat_density_2d(
    aes(fill = Type, alpha = after_stat(level)),
    geom = "polygon",
    contour = TRUE,
    bins = 15,
    colour = NA,
    h = h_vec      # <<< FIXED FOR ALL CASES
  ) +
  
  # Contour lines
  stat_density_2d(
    aes(color = Type),
    contour = TRUE,
    bins = 10,
    linewidth = 0.25,
    h = h_vec      # <<< FIXED FOR ALL CASES
  ) +
  
  scale_fill_manual(values = type_cols, name = "Model type") +
  scale_color_manual(values = type_cols, guide = "none") +
  scale_alpha(range = c(0.10, 0.55), guide = "none") +
  
  labs(
    title = paste("Heatmap —", selected_case, "(n =", selected_n, ")"),
    x = expression(italic(H)),
    y = expression(italic(C))
  ) +
  
  coord_cartesian(xlim = c(H_min, H_max),
                  ylim = c(min(df_case$C)-0.02, max(df_case$C)+0.02)) +
  
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

print(p)

# Save plot
output_plot_path <- file.path(
  "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Plots/ARMA plots/ARMA Plots n5000_n10000",
  selected_case, "Plots",
  paste0("HC_Heatmap_", selected_case, "_n", selected_n, ".pdf")
)
ggsave(output_plot_path, p, width = 10, height = 8)

#-------------------------------------------------------------------------
# End of example_code_heatmap.R for Case 2
#-------------------------------------------------------------------------
# In this script we generate a heatmap of the Entropy-Complexity plane
# for a selected case (Case3) and sample size (n = 5000 or n = 10000).
# The heatmap visualizes the density of points corresponding to different
# time series models (AR, MA, ARMA) in the Entropy-Complexity space.
#--------------------------------------------------------
# Libraries
#--------------------------------------------------------
library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)
library(StatOrdPattHxC)

#--------------------------------------------------------
# Paths
#--------------------------------------------------------
data_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/ARMA Time series results n5000_n10000/HC_Results_all_Models_n5000_n10000.xlsx"

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
# Choose Case and sample size
#--------------------------------------------------------
selected_case <- "Case3"
selected_n    <- 5000

models_case <- cases[[selected_case]]

#--------------------------------------------------------
# Read Excel sheets
#--------------------------------------------------------
df_list <- lapply(models_case, function(m) {
  sheet_name <- paste0(m, "_n", selected_n)
  read_excel(data_path, sheet = sheet_name) |>
    mutate(Model = m)
})

df_case <- bind_rows(df_list) |>
  rename(H = H_Shannon, C = C_Shannon) |>
  filter(is.finite(H), is.finite(C))

#--------------------------------------------------------
# Classify model type
#--------------------------------------------------------
df_case <- df_case |>
  mutate(Type = case_when(
    str_starts(Model, "ARMA") ~ "ARMA",
    str_starts(Model, "AR")   ~ "AR",
    str_starts(Model, "MA")   ~ "MA",
    TRUE                      ~ "Other"
  ))

df_case$Type <- factor(df_case$Type, levels = c("AR","MA","ARMA","Other"))

#--------------------------------------------------------
# Feasible region (crop to H-range)
#--------------------------------------------------------
data("LinfLsup")
Linf <- subset(LinfLsup, Side == "Lower" & Dimension == 3)
Lsup <- subset(LinfLsup, Side == "Upper" & Dimension == 3)

H_min <- min(df_case$H)
H_max <- max(df_case$H)

Linf_crop <- Linf |> filter(H >= H_min, H <= H_max)
Lsup_crop <- Lsup |> filter(H >= H_min, H <= H_max)

# Boundary polygon (CROPPED)
boundary_df <- rbind(
  Linf_crop[, c("H","C")],
  Lsup_crop[order(Lsup_crop$H, decreasing = TRUE), c("H","C")],
  Linf_crop[1, c("H","C")]
)
boundary_df <- as.data.frame(boundary_df)

#--------------------------------------------------------
# Dynamic bandwidth for density smoothing
#--------------------------------------------------------
Hx <- diff(range(df_case$H))
Cx <- diff(range(df_case$C))

h_vec <- c(Hx / 5, Cx / 5)   # adaptive smoothing

#--------------------------------------------------------
# Colors
#--------------------------------------------------------
type_cols <- c(
  "AR"   = "#21918c",
  "MA"   = "#3182bd",
  "ARMA" = "#b2182b",
  "Other" = "grey50"
)

#--------------------------------------------------------
# Plot
#--------------------------------------------------------
p <- ggplot(df_case, aes(H, C)) +
  
  # Feasible region boundaries
  geom_line(data = Linf_crop, aes(H, C),
            colour = "black", linewidth = 0.8) +
  geom_line(data = Lsup_crop, aes(H, C),
            colour = "black", linewidth = 0.8) +
  
  # Background points
  geom_point(alpha = 0.07, size = 0.8, colour = "black") +
  
  # Density polygons (dynamic smoothing!)
  stat_density_2d(
    aes(fill = Type, alpha = after_stat(level)),
    geom = "polygon",
    contour = TRUE,
    bins = 15,
    colour = NA,
    h = h_vec      # <<< FIXED FOR ALL CASES
  ) +
  
  # Contour lines
  stat_density_2d(
    aes(color = Type),
    contour = TRUE,
    bins = 10,
    linewidth = 0.25,
    h = h_vec      # <<< FIXED FOR ALL CASES
  ) +
  
  scale_fill_manual(values = type_cols, name = "Model type") +
  scale_color_manual(values = type_cols, guide = "none") +
  scale_alpha(range = c(0.10, 0.55), guide = "none") +
  
  labs(
    title = paste("Heatmap —", selected_case, "(n =", selected_n, ")"),
    x = expression(italic(H)),
    y = expression(italic(C))
  ) +
  
  coord_cartesian(xlim = c(H_min, H_max),
                  ylim = c(min(df_case$C)-0.02, max(df_case$C)+0.02)) +
  
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

print(p)

# Save plot
output_plot_path <- file.path(
  "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Plots/ARMA plots/ARMA Plots n5000_n10000",
  selected_case, "Plots",
  paste0("HC_Heatmap_", selected_case, "_n", selected_n, ".pdf")
)
ggsave(output_plot_path, p, width = 10, height = 8)
#-------------------------------------------------------------------------
# End of example_code_heatmap.R for Case 3
#-------------------------------------------------------------------------


# In this script we generate a heatmap of the Entropy-Complexity plane
# for a selected cases (Case1-3) and sample size (n = 5000 and 10000).
# The heatmap visualizes the density of points corresponding to different
# time series models (AR, MA, ARMA) in the Entropy-Complexity space.
#--------------------------------------------------------
#--------------------------------------------------------
# Libraries
#--------------------------------------------------------
library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)
library(StatOrdPattHxC)

#--------------------------------------------------------
# Paths
#--------------------------------------------------------
data_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/ARMA Time series results n5000_n10000/HC_Results_all_Models_n5000_n10000.xlsx"
base_output_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Plots/ARMA plots/ARMA Plots n5000_n10000"

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
sample_sizes <- c(5000, 10000)

#--------------------------------------------------------
# Feasible region (Linf / Lsup) for D = 3
#--------------------------------------------------------
data("LinfLsup")
Linf_all <- subset(LinfLsup, Side == "Lower" & Dimension == "3")
Lsup_all <- subset(LinfLsup, Side == "Upper" & Dimension == "3")

#--------------------------------------------------------
# Colors
#--------------------------------------------------------
type_cols <- c(
  "AR"   = "#21918c",
  "MA"   = "#3182bd",
  "ARMA" = "#b2182b",
  "Other" = "grey50"
)

#--------------------------------------------------------
# Main loops: over cases and sample sizes
#--------------------------------------------------------
for (case_name in names(cases)) {
  
  models_case <- cases[[case_name]]
  
  for (n_val in sample_sizes) {
    
    message("===== Processing ", case_name, ", n = ", n_val, " =====")
    
    #--------------------------------------------------------
    # Read Excel sheets for this case & n
    #--------------------------------------------------------
    df_list <- lapply(models_case, function(m) {
      sheet_name <- paste0(m, "_n", n_val)
      message("Reading sheet: ", sheet_name)
      read_excel(data_path, sheet = sheet_name) |>
        mutate(Model = m, n = n_val)
    })
    
    df_case <- bind_rows(df_list) |>
      rename(
        H = H_Shannon,
        C = C_Shannon
      ) |>
      filter(is.finite(H), is.finite(C))
    
    if (nrow(df_case) == 0) {
      warning("No data for ", case_name, " n = ", n_val, ", skipping.")
      next
    }
    
    #--------------------------------------------------------
    # Classify model type: AR / MA / ARMA
    #--------------------------------------------------------
    df_case <- df_case |>
      mutate(Type = case_when(
        str_starts(Model, "ARMA") ~ "ARMA",
        str_starts(Model, "AR")   ~ "AR",
        str_starts(Model, "MA")   ~ "MA",
        TRUE                      ~ "Other"
      ))
    
    df_case$Type <- factor(df_case$Type, levels = c("AR","MA","ARMA","Other"))
    
    #--------------------------------------------------------
    # Feasible region (crop to H-range of this case)
    #--------------------------------------------------------
    H_min <- min(df_case$H)
    H_max <- max(df_case$H)
    C_min <- min(df_case$C)
    C_max <- max(df_case$C)
    
    Linf_crop <- Linf_all |>
      filter(H >= H_min, H <= H_max)
    Lsup_crop <- Lsup_all |>
      filter(H >= H_min, H <= H_max)
    
    #--------------------------------------------------------
    # Dynamic bandwidth for density smoothing
    #--------------------------------------------------------
    Hx <- diff(range(df_case$H))
    Cx <- diff(range(df_case$C))
    h_vec <- c(Hx / 5, Cx / 5)   # adaptive smoothing
    
    #--------------------------------------------------------
    # Plot
    #--------------------------------------------------------
    p <- ggplot(df_case, aes(H, C)) +
      
      # Feasible region boundaries
      geom_line(data = Linf_crop, aes(H, C),
                colour = "black", linewidth = 0.8) +
      geom_line(data = Lsup_crop, aes(H, C),
                colour = "black", linewidth = 0.8) +
      
      # Background points
      geom_point(alpha = 0.07, size = 0.8, colour = "black") +
      
      # Density polygons (heatmap-like)
      stat_density_2d(
        aes(fill = Type, alpha = after_stat(level)),
        geom = "polygon",
        contour = TRUE,
        bins = 15,
        colour = NA,
        h = h_vec
      ) +
      
      # Contour lines
      stat_density_2d(
        aes(color = Type),
        contour = TRUE,
        bins = 10,
        linewidth = 0.25,
        h = h_vec
      ) +
      
      scale_fill_manual(values = type_cols, name = "Model type") +
      scale_color_manual(values = type_cols, guide = "none") +
      scale_alpha(range = c(0.10, 0.55), guide = "none") +
      
      labs(
        title = paste("Heatmap —", case_name, "(n =", n_val, ")"),
        x = expression(italic(H)),
        y = expression(italic(C))
      ) +
      
      coord_cartesian(
        xlim = c(H_min, H_max),
        ylim = c(C_min - 0.02, C_max + 0.02)
      ) +
      
      theme_minimal(base_size = 16) +
      theme(
        panel.grid      = element_blank(),
        legend.position = "right"
      )
    
    print(p)
    
    #--------------------------------------------------------
    # Save plot
    #--------------------------------------------------------
    case_output_dir <- file.path(base_output_dir, case_name, "Plots")
    dir.create(case_output_dir, recursive = TRUE, showWarnings = FALSE)
    
    output_plot_path <- file.path(
      case_output_dir,
      paste0("HC_Heatmap_", case_name, "_n", n_val, ".pdf")
    )
    
    ggsave(output_plot_path, p, width = 10, height = 8)
    message("Saved: ", output_plot_path, "\n")
  }
}

message("🎉 All heatmaps completed for all cases and sample sizes!")
#-------------------------------------------------------------------------
#End of combined heatmap code
#-------------------------------------------------------------------
# In this script we generate heatmaps of the Entropy-Complexity plane
# for all cases (Case1-3) and sample sizes (n = 5000 and 10000).
# The heatmaps visualize the density of points corresponding to different
# time series models (AR, MA, ARMA) in the Entropy-Complexity space.
#--------------------------------------------------------
# Libraries
#--------------------------------------------------------
library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)
library(StatOrdPattHxC)
library(viridis)

#--------------------------------------------------------
# Paths
#--------------------------------------------------------
data_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/ARMA Time series results n5000_n10000/HC_Results_all_Models_n5000_n10000.xlsx"
base_output_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Plots/ARMA plots/ARMA Plots n5000_n10000"

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
sample_sizes <- c(5000, 10000)

#--------------------------------------------------------
# Feasible region (Linf / Lsup) for D = 3
#--------------------------------------------------------
data("LinfLsup")
Linf_all <- subset(LinfLsup, Side == "Lower" & Dimension == "3")
Lsup_all <- subset(LinfLsup, Side == "Upper" & Dimension == "3")

#--------------------------------------------------------
# Main loops: over cases and sample sizes
#--------------------------------------------------------
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
    
    df_case <- bind_rows(df_list) |>
      rename(
        H = H_Shannon,
        C = C_Shannon
      ) |>
      filter(is.finite(H), is.finite(C))
    
    if (nrow(df_case) == 0) {
      warning("No data for ", case_name, " n = ", n_val, ", skipping.")
      next
    }
    
    #--------------------------------------------------------
    # Classify Families (AR / MA / ARMA) and fine Types (each model)
    #--------------------------------------------------------
    df_case <- df_case |>
      mutate(
        Family = case_when(
          str_starts(Model, "ARMA") ~ "ARMA",
          str_starts(Model, "AR")   ~ "AR",
          str_starts(Model, "MA")   ~ "MA",
          TRUE                      ~ "Other"
        ),
        # Fine-grained type: each model in this case is its own "type"
        Type = factor(Model, levels = models_case)
      )
    
    df_case$Family <- factor(df_case$Family,
                             levels = c("AR","MA","ARMA","Other"))
    
    #--------------------------------------------------------
    # Feasible region (crop to H-range of this case)
    #--------------------------------------------------------
    H_min <- min(df_case$H)
    H_max <- max(df_case$H)
    C_min <- min(df_case$C)
    C_max <- max(df_case$C)
    
    Linf_crop <- Linf_all |>
      filter(H >= H_min, H <= H_max)
    Lsup_crop <- Lsup_all |>
      filter(H >= H_min, H <= H_max)
    
    #--------------------------------------------------------
    # Dynamic bandwidth for density smoothing
    #--------------------------------------------------------
    Hx <- diff(range(df_case$H))
    Cx <- diff(range(df_case$C))
    h_vec <- c(Hx / 5, Cx / 5)   # adaptive smoothing
    
    #--------------------------------------------------------
    # Colors per fine Type (each model gets its own color)
    #--------------------------------------------------------
    type_cols <- setNames(
      viridis(length(models_case)),
      models_case
    )
    
    #--------------------------------------------------------
    # Plot
    #--------------------------------------------------------
    p <- ggplot(df_case, aes(H, C)) +
      
      # Feasible region boundaries
      geom_line(data = Linf_crop, aes(H, C),
                colour = "black", linewidth = 0.8) +
      geom_line(data = Lsup_crop, aes(H, C),
                colour = "black", linewidth = 0.8) +
      
      # Background points
      geom_point(alpha = 0.07, size = 0.8, colour = "black") +
      
      # Density polygons per model Type
      stat_density_2d(
        aes(fill = Type, alpha = after_stat(level)),
        geom = "polygon",
        contour = TRUE,
        bins = 15,
        colour = NA,
        h = h_vec
      ) +
      
      # Contour lines per model Type
      stat_density_2d(
        aes(color = Type),
        contour = TRUE,
        bins = 10,
        linewidth = 0.25,
        h = h_vec
      ) +
      
      scale_fill_manual(values = type_cols, name = "Model (Type)") +
      scale_color_manual(values = type_cols, guide = "none") +
      scale_alpha(range = c(0.10, 0.55), guide = "none") +
      
      labs(
        title = paste("Heatmap by Model Type —", case_name, "(n =", n_val, ")"),
        x = expression(italic(H)),
        y = expression(italic(C))
      ) +
      
      coord_cartesian(
        xlim = c(H_min, H_max),
        ylim = c(C_min - 0.02, C_max + 0.02)
      ) +
      
      theme_minimal(base_size = 16) +
      theme(
        panel.grid      = element_blank(),
        legend.position = "right"
      )
    
    print(p)
    
    #--------------------------------------------------------
    # Save plot
    #--------------------------------------------------------
    case_output_dir <- file.path(base_output_dir, case_name, "Plots")
    dir.create(case_output_dir, recursive = TRUE, showWarnings = FALSE)
    
    output_plot_path <- file.path(
      case_output_dir,
      paste0("HC_Heatmap_ModelTypes_", case_name, "_n", n_val, ".pdf")
    )
    
    ggsave(output_plot_path, p, width = 10, height = 8)
    message("  Saved: ", output_plot_path, "\n")
  }
}

message("🎉 All heatmaps completed for all cases, sample sizes, and model types!")
