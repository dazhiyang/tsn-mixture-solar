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
plot.size <- 8
line.size <- 0.3
point.size <- 0.05
legend.size <- 0.4
text.size <- plot.size * 5 / 14
# dir0 <- "/Users/dyang/Dropbox/Working papers/QC"
dir0 <- "/Users/seryangd/Library/CloudStorage/Dropbox/Working papers/distribution"
################################################################################

x <- seq(-3, 3, by = 0.01)
alpha <- c(-4, -1, 0, 1, 1.5)

data.plot <- NULL
for (i in seq_along(alpha))
{
  y.tmp <- dsn(x, xi = 0, omega = 1, alpha = alpha[i])
  data.plot <- data.plot %>%
    bind_rows(., tibble(x = x, y = y.tmp, group = paste0("italic(alpha)==", alpha[i])))
}

data.plot$group <- factor(data.plot$group, levels = paste0("italic(alpha)==", alpha))

p1 <- ggplot() +
  geom_line(aes(x = x, y = y, linetype = group), data = data.plot, linewidth = line.size) +
  scale_linetype_manual(values = c("solid", "longdash", "dashed", "dotdash", "dotted"), labels = scales::parse_format()) +
  # scale_color_manual(values = colorblind_pal()(8)[1:5], labels = scales::parse_format()) +
  scale_x_continuous(name = expression(paste(italic(x)))) +
  # scale_x_datetime(name = "Year", breaks = seq(ymd_h("2016-01-01 00"), ymd_h("2020-01-01 00"), by = "1 year"), labels = date_format("%Y"), expand = c(0.07,0)) +
  scale_y_continuous(name = expression(paste(italic(f)(italic(x))))) +
  theme_bw() +
  theme(plot.margin = unit(c(0.1, 0.2, 0, 0.2), "lines"), text = element_text(family = "Times", size = plot.size), axis.text = element_text(size = plot.size), panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "gray80", linewidth = line.size / 2), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), strip.text.x = element_text(size = plot.size, margin = unit(c(0.1, 0, 0.1, 0), "lines")), strip.text.y = element_text(size = plot.size, margin = unit(c(0, 0.1, 0, 0.1), "lines")), panel.spacing = unit(0, "lines"), legend.position = "inside", legend.position.inside = c(0.85, 0.56), legend.title = element_blank(), legend.text = element_text(size = plot.size), legend.background = element_rect(fill = "transparent"), legend.key = element_rect(fill = "transparent"), legend.box.margin = unit(c(-0.7, 0, 0, 0.0), "lines"))

p1

data.plot <- NULL
for (i in seq_along(alpha))
{
  y.tmp <- psn(x, xi = 0, omega = 1, alpha = alpha[i])
  data.plot <- data.plot %>%
    bind_rows(., tibble(x = x, y = y.tmp, group = paste0("italic(alpha)==", alpha[i])))
}

data.plot$group <- factor(data.plot$group, levels = paste0("italic(alpha)==", alpha))

p2 <- ggplot() +
  geom_line(aes(x = x, y = y, linetype = group), data = data.plot, linewidth = line.size) +
  scale_linetype_manual(values = c("solid", "longdash", "dashed", "dotdash", "dotted"), labels = scales::parse_format()) +
  # scale_color_manual(values = colorblind_pal()(8)[1:5], labels = scales::parse_format()) +
  scale_x_continuous(name = expression(paste(italic(x)))) +
  # scale_x_datetime(name = "Year", breaks = seq(ymd_h("2016-01-01 00"), ymd_h("2020-01-01 00"), by = "1 year"), labels = date_format("%Y"), expand = c(0.07,0)) +
  scale_y_continuous(name = expression(paste(italic(F)(italic(x))))) +
  theme_bw() +
  theme(plot.margin = unit(c(0.1, 0.2, 0, 0.2), "lines"), text = element_text(family = "Times", size = plot.size), axis.text = element_text(size = plot.size), panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "gray80", linewidth = line.size / 2), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), strip.text.x = element_text(size = plot.size, margin = unit(c(0.1, 0, 0.1, 0), "lines")), strip.text.y = element_text(size = plot.size, margin = unit(c(0, 0.1, 0, 0.1), "lines")), panel.spacing = unit(0, "lines"), legend.position = "inside", legend.position.inside = c(0.85, 0.4), legend.title = element_blank(), legend.text = element_text(size = plot.size), legend.background = element_rect(fill = "transparent"), legend.key = element_rect(fill = "transparent"), legend.box.margin = unit(c(-0.7, 0, 0, 0.0), "lines"))

p2

p <- ggpubr::ggarrange(p1, p2, ncol = 1, align = "v", labels = c("(a)", "(b)"), heights = c(1, 1), font.label = list(size = plot.size, color = "black", face = "plain", family = "Times"))

setwd(file.path(dir0, "tex"))
ggsave(filename = "SN.pdf", plot = p, width = 80, height = 90, unit = "mm")
