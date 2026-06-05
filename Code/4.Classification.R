################################################################################
# This code is written by Dazhi Yang
# School of Electrical Engineering and Automation
# Harbin Institute of Technology
# emails: yangdazhi.nus@gmail.com
################################################################################

rm(list = ls(all = TRUE))
libs <- c("dplyr", "lubridate", "sn", "ggplot2", "ComplexHeatmap", "viridis", "raster")
invisible(lapply(libs, library, character.only = TRUE))

################################################################################
# global input
################################################################################
plot.size <- 8
line.size <- 0.2
point.size <- 1.4
legend.size <- 0.4
text.size <- plot.size * 5 / 14
# dir0 <- "/Users/dyang/Dropbox/Working papers/distribution"
dir0 <- "/Users/seryangd/Library/CloudStorage/Dropbox/Working papers/distribution"

# source the necessary functions for fitting mixture models
source(file.path(dir0, "Code/0.Functions.R"))

# get the variables of interest
rvs <- c("kappa", "kt", "kb", "kd")
################################################################################

# read the metadata of sites
setwd(file.path(dir0, "Data"))
locations <- read.csv("location.csv", header = TRUE)

# access the folder contained arranged data from 126 stations
setwd(file.path(dir0, "Data/Results"))
files <- dir(pattern = "\\.RData$", ignore.case = TRUE)

stations <- list()
pb <- txtProgressBar(max = length(files), style = 3) # progress bar
for (i in seq_along(files))
{
  load(files[i])
  data <- data.plot %>%
    filter(group == "3-component TSN")

  for (j in seq_along(rvs))
  {
    tmp <- data %>% filter(quantity == rvs[j])
    assign(rvs[j], list(x = tmp$x, y = tmp$y))
  }

  stn <- substr(files[i], 1, 3)
  stations[[stn]] <- list(kappa = kappa, kt = kt, kb = kb, kd = kd)

  setTxtProgressBar(pb, i)
} # end i
close(pb)

# Create distance matrix for all stations
n_stations <- length(stations)
dist_matrix <- matrix(0, nrow = n_stations, ncol = n_stations)
for (i in 1:(n_stations - 1)) {
  for (j in (i + 1):n_stations) {
    dist_matrix[i, j] <- calculate_station_distance(stations[[i]], stations[[j]], weights = c(0.5, 0.5, 1, 1))
    dist_matrix[j, i] <- dist_matrix[i, j] # Symmetric
  }
}

# Convert to distance object
dist_object <- as.dist(dist_matrix)

# Get the clustering that heatmap uses
heatmap_result <- heatmap(dist_matrix, symm = TRUE)

# First, create a mapping from station order to cluster number
heatmap_order <- heatmap_result$rowInd # This is the order stations appear in heatmap
cluster_assignments <- c(
  rep(1, 14), # First 14 stations in heatmap -> Cluster 1
  rep(2, 28), # Next 28 stations -> Cluster 2
  rep(3, 36), # Next 36 stations -> Cluster 3
  rep(4, 25), # Next 25 stations -> Cluster 4
  rep(5, 23) # Last 23 stations -> Cluster 5
)

# Create a vector that maps station number to cluster
class <- numeric(126)
for (i in 1:126) {
  station_number <- heatmap_order[i] # The actual station number at position i
  class[station_number] <- cluster_assignments[i]
}

################################################################################
# plot the heatmap using Heatmap
################################################################################

rownames(dist_matrix) <- locations$stn
colnames(dist_matrix) <- locations$stn
# Create quantile-based color function
create_quantile_colors <- function(matrix, n_colors = 256) {
  vals <- as.vector(matrix)
  quantiles <- quantile(vals, probs = seq(0, 1, length.out = n_colors), na.rm = TRUE)
  colors <- viridis::viridis(n_colors)
  circlize::colorRamp2(quantiles, colors)
}

# Use it
col_fun <- create_quantile_colors(dist_matrix)

# Create symmetric group labels
group_labels <- c("C4", "C3", "C1", "C5", "C2")
group_colors <- ggthemes::colorblind_pal()(8)[2:6] # Same colors for both sides

# simple version
Heatmap(dist_matrix,
  name = "WD",
  col = col_fun,
  # Disable row and column names
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE
)

# paper version
Heatmap(dist_matrix,
  name = "WD",
  col = col_fun,
  # Change dendrogram line size
  row_dend_width = unit(1.3, "cm"), # Width of row dendrogram
  column_dend_height = unit(1.3, "cm"), # Height of column dendrogram
  row_dend_gp = gpar(lwd = 0.3), # Row dendrogram line width
  column_dend_gp = gpar(lwd = 0.3), # Column dendrogram line width
  # Disable row and column names
  show_row_names = FALSE,
  show_column_names = FALSE,
  # # Top annotation (columns)
  top_annotation = HeatmapAnnotation(
    foo = anno_block(
      gp = gpar(fill = group_colors),
      labels = group_labels,
      height = unit(0.3, "cm"),
      labels_gp = gpar(col = "white", fontsize = plot.size - 1, fontfamily = "Times")
    )
  ),
  # Left annotation (rows) - SYMMETRIC to top
  left_annotation = rowAnnotation(
    foo = anno_block(
      gp = gpar(fill = group_colors),
      labels = group_labels,
      width = unit(0.3, "cm"),
      labels_gp = gpar(col = "white", fontsize = plot.size - 1, fontfamily = "Times")
    )
  ),
  column_split = class,
  row_split = class,
  column_title = NULL,
  row_title = NULL,
  # Same font settings for both axes
  column_names_gp = gpar(fontsize = plot.size - 1, fontfamily = "Times"),
  row_names_gp = gpar(fontsize = plot.size - 1, fontfamily = "Times"),
  # Legend settings
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = plot.size - 1, fontfamily = "Times", fontface = "plain"),
    labels_gp = gpar(fontsize = plot.size - 1, fontfamily = "Times"),
    legend_height = unit(4, "cm")
  )
)

################################################################################
# plot a world map showing clustering results
################################################################################
kgc <- raster(file.path(dir0, "Data/Beck_KG_V1_present_0p083.tif"))
kgc <- as_tibble(rasterToPoints(kgc))
names(kgc) <- c("lon", "lat", "class")
kgc <- kgc %>%
  filter(class > 0)

kgc$class <- factor(kgc$class, levels = c(1:30))
levels(kgc$class) <- c(
  "Af", "Am", "As",
  "BWh", "BWk", "BSh", "BSk",
  "Csa", "Csb", "Csc",
  "Cwa", "Cwb", "Cwc",
  "Cfa", "Cfb", "Cfc",
  "Dsa", "Dsb", "Dsc", "Dsd",
  "Dwa", "Dwb", "Dwc", "Dwd",
  "Dfa", "Dfb", "Dfc", "Dfd",
  "ET", "EF"
)

kp.color30 <- c(
  "#1100FF", "#0F72FA", "#3397E4",
  "#FC0200", "#F4998F", "#F0A000", "#EBC35B",
  "#FDFF03", "#D0C708", "#8E8F0A",
  "#8BFF96", "#5BC861", "#349431",
  "#C3FB48", "#62FD33", "#37C800",
  "#FC00FF", "#CB00C4", "#98329F", "#8F5A92",
  "#9FB5FF", "#4377DC", "#4851AE", "#2D028C",
  "#00FAFF", "#43C2FF", "#057B7B", "#004665",
  "#9F9F9F", "#5F615E"
) # "#69B0E5",

loc <- locations %>%
  mutate(RCC = class)
loc$RCC <- factor(loc$RCC, levels = 1:5, labels = paste0("C", 1:5))

p <- ggplot() +
  geom_raster(data = kgc, aes(lon, lat, fill = class)) +
  geom_point(data = loc, aes(x = lon, y = lat, shape = RCC), size = point.size, stroke = 0.2) +
  # geom_polygon(data = cn, aes(long, lat, group=group), fill = NA, color = "gray50", size = line.size, alpha = 0.7) +
  coord_quickmap(expand = 0) +
  scale_fill_manual(name = "Koppen-Geiger climate classification system", values = kp.color30) +
  scale_shape_manual(name = "Radiation climate classification", values = 1:5) +
  theme_minimal() +
  theme(plot.margin = unit(c(0, 0, 0, -0.1), "lines"), panel.spacing = unit(0, "lines"), panel.grid = element_blank(), text = element_text(family = "Times", size = plot.size), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), strip.text.x = element_text(margin = margin(0, 0, 0, 0, "lines"), size = plot.size), legend.position = "bottom", legend.direction = "horizontal", legend.text = element_text(family = "Times", size = plot.size, color = "black"), legend.title = element_text(family = "Times", size = plot.size, color = "black"), legend.key.height = unit(0.2, "lines"), legend.box.margin = unit(c(-0.7, 0, 0, 0), "lines"), legend.background = element_rect(fill = "transparent", colour = "transparent"), legend.box = "vertical", legend.spacing.y = unit(0, "lines")) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5), shape = guide_legend(title.position = "top", title.hjust = 0.5))

p

setwd(file.path(dir0, "Revision1"))
ggsave(filename = "mapStn.pdf", plot = p, width = 85, height = 90, unit = "mm")


################################################################################
# plot the densities after clustering
################################################################################
setwd(file.path(dir0, "Data/Results"))

# a loop to store the kappa values for plotting
data.line <- NULL
pb <- txtProgressBar(max = length(files), style = 3) # progress bar
for (i in seq_along(files))
{
  load(files[i])
  kappa <- pull(data.hist %>% filter(quantity == "kappa"), "x")
  kt <- pull(data.hist %>% filter(quantity == "kt"), "x")
  kb <- pull(data.hist %>% filter(quantity == "kb"), "x")
  kd <- pull(data.hist %>% filter(quantity == "kd"), "x")

  # Calculate the KDE
  kde_kappa <- density(kappa)
  kde_kt <- density(kt)
  kde_kb <- density(kb)
  kde_kd <- density(kd)

  # Get the KDE value for a specific point
  x_new <- seq(0, 1.65, by = 0.01)
  kappa_new <- approx(kde_kappa$x, kde_kappa$y, xout = x_new)$y
  kt_new <- approx(kde_kt$x, kde_kt$y, xout = x_new)$y
  kb_new <- approx(kde_kb$x, kde_kb$y, xout = x_new)$y
  kd_new <- approx(kde_kd$x, kde_kd$y, xout = x_new)$y

  data.line <- data.line %>%
    bind_rows(., tibble(x = x_new, y = kappa_new, stn = substr(files[i], 1, 3), RCC = class[i], index = "kappa")) %>%
    bind_rows(., tibble(x = x_new, y = kt_new, stn = substr(files[i], 1, 3), RCC = class[i], index = "italic(k)[italic(t)]")) %>%
    bind_rows(., tibble(x = x_new, y = kb_new, stn = substr(files[i], 1, 3), RCC = class[i], index = "italic(k)[italic(b)]")) %>%
    bind_rows(., tibble(x = x_new, y = kd_new, stn = substr(files[i], 1, 3), RCC = class[i], index = "italic(k)[italic(d)]"))

  setTxtProgressBar(pb, i)
}
close(pb)

# distribution plot
p1 <- ggplot() +
  geom_line(aes(x = x, y = y, group = stn), data = data.line %>% filter(index == "kappa"), linewidth = line.size, alpha = 0.5) +
  # scale_linetype_manual(values=c("dashed", "dotdash", "solid")) +
  # scale_color_manual(values = colorblind_pal()(8)[2:4]) +
  facet_wrap(~RCC, nrow = 1) +
  scale_x_continuous(name = expression(paste("Clear-sky index, ", kappa, " [dimensionless]")), breaks = c(0, 0.5, 1, 1.5)) +
  scale_y_continuous(name = expression(paste("KDE [dimensionless]")), limits = c(0, 15)) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2, 0.2, 0, 0.4), "lines"), text = element_text(family = "Times", size = plot.size), axis.text = element_text(size = plot.size), panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "gray80", linewidth = line.size / 2), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), strip.text.x = element_text(size = plot.size, margin = unit(c(0.1, 0, 0.1, 0), "lines")), strip.text.y = element_text(size = plot.size, margin = unit(c(0, 0.1, 0, 0.1), "lines")), panel.spacing = unit(0, "lines"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = plot.size), legend.background = element_rect(fill = "transparent"), legend.key = element_rect(fill = "transparent"), legend.box.margin = unit(c(-0.7, 0, 0, 0.0), "lines"), legend.direction = "horizontal") # legend.position = "inside", legend.position.inside = c(0.3, 0.7),

p1

p2 <- ggplot() +
  geom_line(aes(x = x, y = y, group = stn), data = data.line %>% filter(index == "italic(k)[italic(t)]"), linewidth = line.size, alpha = 0.5) +
  # scale_linetype_manual(values=c("dashed", "dotdash", "solid")) +
  # scale_color_manual(values = colorblind_pal()(8)[2:4]) +
  facet_wrap(~RCC, nrow = 1) +
  scale_x_continuous(name = expression(paste("Clearness index, ", italic(k)[italic(t)], " [dimensionless]")), limits = c(0, 1.2), breaks = c(0, 0.5, 1)) +
  scale_y_continuous(name = expression(paste("KDE [dimensionless]"))) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2, 0.2, 0, 0.4), "lines"), text = element_text(family = "Times", size = plot.size), axis.text = element_text(size = plot.size), panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "gray80", linewidth = line.size / 2), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), strip.text.x = element_text(size = plot.size, margin = unit(c(0.1, 0, 0.1, 0), "lines")), strip.text.y = element_text(size = plot.size, margin = unit(c(0, 0.1, 0, 0.1), "lines")), panel.spacing = unit(0, "lines"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = plot.size), legend.background = element_rect(fill = "transparent"), legend.key = element_rect(fill = "transparent"), legend.box.margin = unit(c(-0.7, 0, 0, 0.0), "lines"), legend.direction = "horizontal") # legend.position = "inside", legend.position.inside = c(0.3, 0.7),

p2

p3 <- ggplot() +
  geom_line(aes(x = x, y = y, group = stn), data = data.line %>% filter(index == "italic(k)[italic(b)]"), linewidth = line.size, alpha = 0.5) +
  # scale_linetype_manual(values=c("dashed", "dotdash", "solid")) +
  # scale_color_manual(values = colorblind_pal()(8)[2:4]) +
  facet_wrap(~RCC, nrow = 1) +
  scale_x_continuous(name = expression(paste("Clearness index, ", italic(k)[italic(t)], " [dimensionless]")), limits = c(0, 0.9), breaks = c(0, 0.4, 0.8)) +
  scale_y_continuous(name = expression(paste("KDE [dimensionless]"))) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2, 0.2, 0, 0.4), "lines"), text = element_text(family = "Times", size = plot.size), axis.text = element_text(size = plot.size), panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "gray80", linewidth = line.size / 2), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), strip.text.x = element_text(size = plot.size, margin = unit(c(0.1, 0, 0.1, 0), "lines")), strip.text.y = element_text(size = plot.size, margin = unit(c(0, 0.1, 0, 0.1), "lines")), panel.spacing = unit(0, "lines"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = plot.size), legend.background = element_rect(fill = "transparent"), legend.key = element_rect(fill = "transparent"), legend.box.margin = unit(c(-0.7, 0, 0, 0.0), "lines"), legend.direction = "horizontal") # legend.position = "inside", legend.position.inside = c(0.3, 0.7),

p3

p4 <- ggplot() +
  geom_line(aes(x = x, y = y, group = stn), data = data.line %>% filter(index == "italic(k)[italic(d)]"), linewidth = line.size, alpha = 0.5) +
  # scale_linetype_manual(values=c("dashed", "dotdash", "solid")) +
  # scale_color_manual(values = colorblind_pal()(8)[2:4]) +
  facet_wrap(~RCC, nrow = 1) +
  scale_x_continuous(name = expression(paste("Diffuse transmittance, ", italic(k)[italic(t)], " [dimensionless]")), limits = c(0, 1.2), breaks = c(0, 0.5, 1)) +
  scale_y_continuous(name = expression(paste("KDE [dimensionless]"))) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2, 0.2, 0, 0.4), "lines"), text = element_text(family = "Times", size = plot.size), axis.text = element_text(size = plot.size), panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "gray80", linewidth = line.size / 2), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), strip.text.x = element_text(size = plot.size, margin = unit(c(0.1, 0, 0.1, 0), "lines")), strip.text.y = element_text(size = plot.size, margin = unit(c(0, 0.1, 0, 0.1), "lines")), panel.spacing = unit(0, "lines"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = plot.size), legend.background = element_rect(fill = "transparent"), legend.key = element_rect(fill = "transparent"), legend.box.margin = unit(c(-0.7, 0, 0, 0.0), "lines"), legend.direction = "horizontal") # legend.position = "inside", legend.position.inside = c(0.3, 0.7),

p4
