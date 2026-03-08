################################################################################
# This code is written by Dazhi Yang
# School of Electrical Engineering and Automation
# Harbin Institute of Technology
# emails: yangdazhi.nus@gmail.com
################################################################################

rm(list = ls(all = TRUE))
libs <- c("dplyr", "sn", "ggplot2", "ggthemes")
invisible(lapply(libs, library, character.only = TRUE))

################################################################################
# global input
################################################################################
# source the necessary functions for fitting mixture models
source(file.path(dir0, "Code/0.Functions.R"))

plot.size <- 8
line.size <- 0.3
point.size <- 0.05
legend.size <- 0.4
text.size <- plot.size * 5 / 14
# dir0 <- "/Users/dyang/Dropbox/Working papers/QC"
dir0 <- "/Users/seryangd/Library/CloudStorage/Dropbox/Working papers/distribution"
################################################################################
# Parameters
xi <- 0.5
omega <- 0.3
a <- 0 # lower truncation bound
b <- 1 # upper truncation bound
x <- seq(-0.5, 1.5, by = 0.001)
alpha <- c(-2, 0, 5)

data.plot <- NULL
for (i in seq_along(alpha))
{
  # untruncated version
  tmp <- dsn(x, xi = xi, omega = omega, alpha = alpha[i])
  data.plot <- data.plot %>%
    bind_rows(., tibble(x = x, y = tmp, group = paste0("xi==", xi, "*\', \'*~omega==", omega, "*\', \'*~alpha==", alpha[i]), type = "Untruncated"))
  # truncated version
  tmp <- dtsn(x, xi = xi, omega = omega, alpha = alpha[i], a, b)
  data.plot <- data.plot %>%
    bind_rows(., tibble(x = x, y = tmp, group = paste0("xi==", xi, "*\', \'*~omega==", omega, "*\', \'*~alpha==", alpha[i]), type = "Truncated"))
}

data.plot$group <- factor(data.plot$group, levels = paste0("xi==", xi, "*\', \'*~omega==", omega, "*\', \'*~alpha==", alpha))

p <- ggplot() +
  geom_line(aes(x = x, y = y, linetype = type), data = data.plot, linewidth = line.size) +
  facet_wrap(~group, labeller = label_parsed, nrow = 1) +
  # scale_color_manual(values = colorblind_pal()(8)[1:5], labels = scales::parse_format()) +
  scale_x_continuous(name = expression(paste(italic(x)))) +
  # scale_x_datetime(name = "Year", breaks = seq(ymd_h("2016-01-01 00"), ymd_h("2020-01-01 00"), by = "1 year"), labels = date_format("%Y"), expand = c(0.07,0)) +
  scale_y_continuous(name = expression(paste(italic(f)(italic(x))))) +
  theme_bw() +
  theme(plot.margin = unit(c(0.1, 0.2, 0, 0.2), "lines"), text = element_text(family = "Times", size = plot.size), axis.text = element_text(size = plot.size), panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "gray80", linewidth = line.size / 2), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), strip.text.x = element_text(size = plot.size, margin = unit(c(0.1, 0, 0.1, 0), "lines")), strip.text.y = element_text(size = plot.size, margin = unit(c(0, 0.1, 0, 0.1), "lines")), panel.spacing = unit(0, "lines"), legend.position = "inside", legend.position.inside = c(0.5, 0.75), legend.title = element_blank(), legend.text = element_text(size = plot.size), legend.background = element_rect(fill = "transparent"), legend.key = element_rect(fill = "transparent"), legend.box.margin = unit(c(-0.7, 0, 0, 0.0), "lines"))

p


setwd(file.path(dir0, "tex"))
ggsave(filename = "TSN.pdf", plot = p, width = 160, height = 40, unit = "mm")
