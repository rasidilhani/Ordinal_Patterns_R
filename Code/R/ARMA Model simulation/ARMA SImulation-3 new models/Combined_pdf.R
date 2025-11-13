# ==========================================
# 📄 Combine Four PDF Plots (2x2 Grid Layout)
# ==========================================
library(magick)
library(grid)
library(gridExtra)

# ---- 1. Define your directory ----
plot_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/ARMA Time series results/Case1/Plots"
setwd(plot_dir)

# ---- 2. Define the PDF files to combine ----
files <- c(
  "Scatter_Case1_n500.pdf",
  "CentralPoints_Case1_n500.pdf",
  "CentralPoints_TimeSeries_Case1_n500.pdf"
  #"new_centralpoint_timeseries_Case1_n500.pdf"
)

# ---- 3. Check all files exist ----
missing <- files[!file.exists(files)]
if (length(missing) > 0) stop("❌ Missing files:\n", paste(missing, collapse = "\n"))

# ---- 4. Read each PDF as magick images ----
imgs <- lapply(files, function(f) {
  cat("Reading:", f, "\n")
  image_read_pdf(normalizePath(f), density = 300)
})

# ---- 5. Convert magick images to raster grobs ----
grobs <- lapply(imgs, function(im) {
  # Convert to PNG internally to make it raster-friendly
  im_converted <- image_convert(im, format = "png")
  rasterGrob(as.raster(im_converted), interpolate = FALSE)
})

# ---- 6. Arrange 2x2 layout and save as combined PDF ----
output_file <- file.path(plot_dir, "Combined_3Plots_Case1_n500.pdf")

pdf(output_file, width = 16, height = 12)
grid.arrange(
  grobs[[1]], grobs[[2]],
  grobs[[3]], grobs[[4]],
  nrow = 2, ncol = 2,
  top = textGrob("Combined Entropy–Complexity and Time Series Plots",
                 gp = gpar(fontsize = 16, fontface = "bold"))
)
dev.off()

cat("\n✅ Combined 4-in-1 PDF saved successfully to:\n", output_file, "\n")

###################################################################################
# 3 by 1
# ============================================================
# 📄 Combine Three PDF Plots into a 2×2 Grid (Option A)
# ============================================================
# ============================================================
# 📄 Combine Three PDF Plots into a 2×2 Grid (Option A, FIXED)
# ============================================================
library(magick)
library(grid)
library(gridExtra)

# ---- 1. Define your Case and n ----
case_name <- "Case3"
n_val     <- 500

# ---- 2. Directory containing the plots ----
plot_dir <- file.path(
  "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/ARMA Time series results",
  case_name, "Plots"
)

setwd(plot_dir)

# ---- 3. Define the PDF files to combine ----
files <- c(
  paste0("Scatter_", case_name, "_n", n_val, ".pdf"),
  paste0("CentralPoints_", case_name, "_n", n_val, ".pdf"),
  paste0("CentralPoints_TimeSeries_", case_name, "_n", n_val, ".pdf")
)

# ---- 4. Check for missing files ----
missing <- files[!file.exists(files)]
if (length(missing) > 0) {
  stop("❌ Missing files:\n", paste(missing, collapse = "\n"))
}

# ---- 5. Read PDFs correctly — use only PAGE 1 ----
imgs <- lapply(files, function(f) {
  cat("Reading:", f, "\n")
  im_all <- image_read_pdf(normalizePath(f), density = 300)
  
  # ---- THE FIX: extract only the first page ----
  im_page1 <- im_all[1]
  return(im_page1)
})

# ---- 6. Convert to raster grobs ----
grobs <- lapply(imgs, function(im) {
  im_png <- image_convert(im, format = "png")
  rasterGrob(as.raster(im_png), interpolate = FALSE)
})

# ---- 7. Create an empty panel for the 4th cell ----
empty_panel <- grid.rect(gp = gpar(col = NA, fill = NA))

# ---- 8. Output combined PDF ----
output_file <- file.path(plot_dir,
                         paste0("Combined3Plots_",case_name,"_n",n_val,".pdf"))

pdf(output_file, width = 16, height = 12)
grid.arrange(
  grobs[[1]], grobs[[2]],
  grobs[[3]], empty_panel,
  nrow = 2, ncol = 2,
  top = textGrob(
    paste("Combined3Plots_ (", case_name, ", n=", n_val, ")"),
    gp = gpar(fontsize = 16, fontface = "bold")
  )
)
dev.off()

cat("\n✅ Combined 3-in-1 PDF saved successfully to:\n", output_file, "\n")
