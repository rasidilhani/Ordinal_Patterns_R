# ==========================================
# 📄 Combine Four PDF Plots (2x2 Grid Layout)
# ==========================================
library(magick)
library(grid)
library(gridExtra)

# ---- 1. Define your directory ----
plot_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/ARMA Time series results/Case 3/Plots"
setwd(plot_dir)

# ---- 2. Define the PDF files to combine ----
files <- c(
  "HC_Scatter_Case 3_Shannon.pdf",
  "HC_Scatter_CI_Case 3_n1000.pdf",
  "HC_Heatmap_Case 3_n1000.pdf",
  "Combined_HC_TS_grid_case3.pdf"
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
output_file <- file.path(plot_dir, "Combined_4Plots_Case3.pdf")

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


