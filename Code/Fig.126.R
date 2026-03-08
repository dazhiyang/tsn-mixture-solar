################################################################################
# This code is written by Dazhi Yang
# School of Electrical Engineering and Automation
# Harbin Institute of Technology
# emails: yangdazhi.nus@gmail.com
################################################################################

rm(list = ls(all = TRUE))
libs <- c("dplyr", "lubridate", "mixtools", "sn", "bbmle", "ggplot2", "xtable")
invisible(lapply(libs, library, character.only = TRUE))

################################################################################
# global input
################################################################################
plot.size <- 8
line.size <- 0.2
point.size <- 0.05
legend.size <- 0.4
text.size <- plot.size * 5 / 14
dir0 <- "/Users/dyang/Dropbox/Working papers/distribution"
# dir0 <- "/Users/seryangd/Library/CloudStorage/Dropbox/Working papers/distribution"
dir.data <- "/Volumes/Macintosh Research/Data/Separation/tibbles"

# source the necessary functions for fitting mixture models
source(file.path(dir0, "Code/0.Functions.R"))
################################################################################

# read the metadata of sites
setwd(file.path(dir0, "Data"))
locations <- read.csv("location.csv", header = TRUE)
climates <- read.csv("climatology variables of 126 sites.csv", header = TRUE)

# access the folder contained arranged data from 126 stations
setwd(dir.data)
files <- dir()

# a loop to store the kappa values for plotting
data.plot <- NULL
set.seed(1234) # for reproducibility
pb <- txtProgressBar(max = length(files), style = 3) # progress bar
for (i in seq_along(files))
{
  load(files[i])
  kappa <- data$kc
  kappa <- sample(kappa, size = min(30000, nrow(data)), replace = FALSE)
  kt <- data$kt
  kt <- sample(kt, size = min(30000, nrow(data)), replace = FALSE)

  # Calculate the KDE
  kde_kappa <- density(kappa)
  kde_kt <- density(kt)
  # Get the KDE value for a specific point
  x_new <- seq(0, 1.65, by = 0.01)
  kappa_new <- approx(kde_kappa$x, kde_kappa$y, xout = x_new)$y
  kt_new <- approx(kde_kt$x, kde_kt$y, xout = x_new)$y

  data.plot <- data.plot %>%
    bind_rows(., tibble(x = x_new, y = kappa_new, stn = substr(files[i], 1, 3), KGC = substr(locations$CZ[i], 1, 1), RCC = climates$cluster[i], index = "kappa")) %>%
    bind_rows(., tibble(x = x_new, y = kt_new, stn = substr(files[i], 1, 3), KGC = substr(locations$CZ[i], 1, 1), RCC = climates$cluster[i], index = "italic(k)[italic(t)]"))

  setTxtProgressBar(pb, i)
}
close(pb)

# plotting and tables
data.plot$KGC <- factor(data.plot$KGC, levels = c("A", "B", "C", "D", "E", "O"), labels = c("A (tropical)", "B (arid)", "C (temperate)", "D (continental)", "E (polar)", "Ocean"))
data.plot$RCC <- factor(data.plot$RCC, levels = c("1", "2", "3", "4", "5"), labels = c("C1 (high aerosol)", "C2 (dust & sand)", "C3 (high albedo)", "C4 (clear sky)", "C5 (cloudy)"))

# distribution plot
p1 <- ggplot() +
  geom_line(aes(x = x, y = y, group = stn), data = data.plot %>% filter(index == "kappa"), linewidth = line.size, alpha = 0.5) +
  # scale_linetype_manual(values=c("dashed", "dotdash", "solid")) +
  # scale_color_manual(values = colorblind_pal()(8)[2:4]) +
  facet_wrap(~RCC, nrow = 1) +
  scale_x_continuous(name = expression(paste("Clear-sky index, ", kappa, " [dimensionless]")), breaks = c(0, 0.5, 1, 1.5)) +
  scale_y_continuous(name = expression(paste("KDE [dimensionless]")), limits = c(0, 15)) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2, 0.2, 0, 0.4), "lines"), text = element_text(family = "Times", size = plot.size), axis.text = element_text(size = plot.size), panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "gray80", linewidth = line.size / 2), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), strip.text.x = element_text(size = plot.size, margin = unit(c(0.1, 0, 0.1, 0), "lines")), strip.text.y = element_text(size = plot.size, margin = unit(c(0, 0.1, 0, 0.1), "lines")), panel.spacing = unit(0, "lines"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = plot.size), legend.background = element_rect(fill = "transparent"), legend.key = element_rect(fill = "transparent"), legend.box.margin = unit(c(-0.7, 0, 0, 0.0), "lines"), legend.direction = "horizontal") # legend.position = "inside", legend.position.inside = c(0.3, 0.7),

p2 <- ggplot() +
  geom_line(aes(x = x, y = y, group = stn), data = data.plot %>% filter(index == "italic(k)[italic(t)]"), linewidth = line.size, alpha = 0.5) +
  # scale_linetype_manual(values=c("dashed", "dotdash", "solid")) +
  # scale_color_manual(values = colorblind_pal()(8)[2:4]) +
  facet_wrap(~RCC, nrow = 1) +
  scale_x_continuous(name = expression(paste("Clearness index, ", italic(k)[italic(t)], " [dimensionless]")), limits = c(0, 1.2), breaks = c(0, 0.5, 1)) +
  scale_y_continuous(name = expression(paste("KDE [dimensionless]"))) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2, 0.2, 0, 0.4), "lines"), text = element_text(family = "Times", size = plot.size), axis.text = element_text(size = plot.size), panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "gray80", linewidth = line.size / 2), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), strip.text.x = element_text(size = plot.size, margin = unit(c(0.1, 0, 0.1, 0), "lines")), strip.text.y = element_text(size = plot.size, margin = unit(c(0, 0.1, 0, 0.1), "lines")), panel.spacing = unit(0, "lines"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = plot.size), legend.background = element_rect(fill = "transparent"), legend.key = element_rect(fill = "transparent"), legend.box.margin = unit(c(-0.7, 0, 0, 0.0), "lines"), legend.direction = "horizontal") # legend.position = "inside", legend.position.inside = c(0.3, 0.7),


p <- ggpubr::ggarrange(p1, p2, ncol = 1, align = "v", labels = c("(a)", "(b)"), heights = c(1, 1), font.label = list(size = plot.size, color = "black", face = "plain", family = "Times"))

setwd(file.path(dir0, "tex"))
ggsave(filename = "kappa126.pdf", plot = p, width = 160, height = 80, unit = "mm")
