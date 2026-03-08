################################################################################
# This code is written by Dazhi Yang
# School of Electrical Engineering and Automation
# Harbin Institute of Technology
# emails: yangdazhi.nus@gmail.com
################################################################################

rm(list = ls(all = TRUE))
libs <- c("dplyr", "lubridate", "ggplot2", "rworldmap")
invisible(lapply(libs, library, character.only = TRUE))

################################################################################
# global input
################################################################################
plot.size <- 8
line.size <- 0.2
point.size <- 0.05
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
climates <- read.csv("climatology variables of 126 sites.csv", header = TRUE)
locations <- read.csv("location.csv", header = TRUE)
locations$Continent <- coords2continent(cbind(locations$lon, locations$lat))
locations$Continent <- ifelse(is.na(locations$Continent), "Island", as.character(locations$Continent))
locations$KGC <- substr(locations$CZ, 1, 1)
locations$RCC <- climates$cluster

# access the folder contained arranged data from 126 stations
setwd(file.path(dir0, "Data/Results"))
files <- dir()

mixtures <- c("2-component N", "3-component N", "3-component TSN")
metrics <- c("AIC", "BIC", "WD", "KSD")
error.kappa <- error.kt <- error.kb <- error.kd <- array(NA, dim = c(nrow(locations), 3, 4)) # N stn, 3 mixtures, 4 metrics
pb <- txtProgressBar(max = length(files), style = 3) # progress bar
for (i in seq_along(files))
{
  load(files[i])
  data.perf <- data.perf %>%
    mutate(aic = as.numeric(aic)) %>%
    mutate(bic = as.numeric(bic)) %>%
    mutate(was = latex_sci_to_numeric(was))
  dat <- data.matrix(data.perf[, 5:8])
  for (j in seq_along()(mixtures))
  {
    for (k in seq_along(metrics))
    {
      error.kappa[i, j, k] <- dat[j, k]
      error.kt[i, j, k] <- dat[j + 3, k]
      error.kb[i, j, k] <- dat[j + 6, k]
      error.kd[i, j, k] <- dat[j + 9, k]
    } # end k
  } # end j
  setTxtProgressBar(pb, i)
} # end i
close(pb)

#################################################################################
# Box plot by continents
#################################################################################
continents <- c("Africa", "Antarctica", "Asia", "Australia", "Europe", "North America", "South America", "Island")
data.plot.box <- NULL
for (l in seq_along()(continents))
{
  select <- which(locations$Continent == continents[l])
  names <- locations$stn[select]
  for (j in 1:3) # mixtures
  {
    for (k in 1:4) # metrics
    {
      metric.tmp <- metrics[k]
      mixture.tmp <- mixtures[j]
      continent.tmp <- continents[l]
      data.plot.box <- data.plot.box %>%
        bind_rows(., tibble(error = error.kappa[select, j, k], mixture = mixture.tmp, metric = metric.tmp, continent = continent.tmp, index = "kappa")) %>%
        bind_rows(., tibble(error = error.kt[select, j, k], mixture = mixture.tmp, metric = metric.tmp, continent = continent.tmp, index = "kt")) %>%
        bind_rows(., tibble(error = error.kb[select, j, k], mixture = mixture.tmp, metric = metric.tmp, continent = continent.tmp, index = "kb")) %>%
        bind_rows(., tibble(error = error.kd[select, j, k], mixture = mixture.tmp, metric = metric.tmp, continent = continent.tmp, index = "kd"))
    } # end k
  } # end j
} # end l

data.plot.box$index <- factor(data.plot.box$index, levels = c("kappa", "kt", "kb", "kd"), labels = c("Clear~sky~index*\', \'*~kappa", "Clearness~index*\', \'*~italic(k)[italic(t)]", "Beam~transmittance*\', \'*~italic(k)[italic(b)]", "Diffuse~transmittance*\', \'*~italic(k)[italic(d)]"))
data.plot.box$metric <- factor(data.plot.box$metric, levels = c("AIC", "BIC", "WD", "KSD"))

data.plot.box1 <- data.plot.box %>%
  filter(mixture %in% c("3-component N", "3-component TSN")) %>%
  filter(error < 0.049)

p1 <- ggplot(data.plot.box1) +
  geom_hline(yintercept = 0, color = "black", linetype = "dotted") +
  geom_boxplot(aes(x = continent, y = error, fill = mixture), size = line.size, outlier.size = point.size, varwidth = TRUE) +
  scale_fill_manual(values = c("gray80", "gray50")) +
  facet_wrap(~ metric + index, labeller = label_parsed, ncol = 4, scales = "free_x") +
  coord_flip() +
  theme_minimal() +
  theme(plot.margin = unit(c(0.5, 0.2, 0, 0), "lines"), panel.spacing = unit(0.5, "lines"), text = element_text(family = "Times", size = plot.size), axis.title = element_blank(), strip.text.x = element_text(margin = margin(0, 0, 0, 0, "lines")), plot.title = element_text(family = "Times", size = plot.size, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.box.margin = unit(c(-0.7, 0, 0, 0.0), "lines"), legend.direction = "horizontal") # c(0.92,0.15)

p1

# setwd(file.path(dir0, "tex"))
# ggsave(filename = "PerfConti.pdf", plot = p1, width = 160, height = 120, unit = "mm")

#################################################################################
# Box plot by RCC
#################################################################################
locations$RCC <- factor(locations$RCC, levels = c("1", "2", "3", "4", "5"), labels = c("C1 (high aerosol)", "C2 (dust & sand)", "C3 (high albedo)", "C4 (clear sky)", "C5 (cloudy)"))
RCC <- c("C1 (high aerosol)", "C2 (dust & sand)", "C3 (high albedo)", "C4 (clear sky)", "C5 (cloudy)")
data.plot.box <- NULL
for (l in seq_along(RCC))
{
  select <- which(locations$RCC == RCC[l])
  names <- locations$stn[select]
  for (j in 1:3) # mixtures
  {
    for (k in 1:4) # metrics
    {
      metric.tmp <- metrics[k]
      mixture.tmp <- mixtures[j]
      RCC.tmp <- RCC[l]
      data.plot.box <- data.plot.box %>%
        bind_rows(., tibble(error = error.kappa[select, j, k], mixture = mixture.tmp, metric = metric.tmp, RCC = RCC.tmp, index = "kappa")) %>%
        bind_rows(., tibble(error = error.kt[select, j, k], mixture = mixture.tmp, metric = metric.tmp, RCC = RCC.tmp, index = "kt")) %>%
        bind_rows(., tibble(error = error.kb[select, j, k], mixture = mixture.tmp, metric = metric.tmp, RCC = RCC.tmp, index = "kb")) %>%
        bind_rows(., tibble(error = error.kd[select, j, k], mixture = mixture.tmp, metric = metric.tmp, RCC = RCC.tmp, index = "kd"))
    } # end k
  } # end j
} # end l

data.plot.box$index <- factor(data.plot.box$index, levels = c("kappa", "kt", "kb", "kd"), labels = c("Clear~sky~index*\', \'*~kappa", "Clearness~index*\', \'*~italic(k)[italic(t)]", "Beam~transmittance*\', \'*~italic(k)[italic(b)]", "Diffuse~transmittance*\', \'*~italic(k)[italic(d)]"))
data.plot.box$metric <- factor(data.plot.box$metric, levels = c("AIC", "BIC", "WD", "KSD"))

data.plot.box1 <- data.plot.box %>%
  filter(mixture %in% c("3-component N", "3-component TSN")) %>%
  filter(error < 0.049)

p2 <- ggplot(data.plot.box1) +
  geom_hline(yintercept = 0, color = "black", linetype = "dotted") +
  geom_boxplot(aes(x = RCC, y = error, fill = mixture), size = line.size, outlier.size = point.size, varwidth = TRUE) +
  scale_fill_manual(values = c("gray80", "gray50")) +
  facet_wrap(~ metric + index, labeller = label_parsed, ncol = 4, scales = "free_x") +
  coord_flip() +
  theme_minimal() +
  theme(plot.margin = unit(c(0.5, 0.2, 0, 0), "lines"), panel.spacing = unit(0.5, "lines"), text = element_text(family = "Times", size = plot.size), axis.title = element_blank(), strip.text.x = element_text(margin = margin(0, 0, 0, 0, "lines")), plot.title = element_text(family = "Times", size = plot.size, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.box.margin = unit(c(-0.7, 0, 0, 0.0), "lines"), legend.direction = "horizontal") # c(0.92,0.15)

p2

setwd(file.path(dir0, "tex"))
ggsave(filename = "PerfRCC.pdf", plot = p2, width = 160, height = 120, unit = "mm")


#################################################################################
# Kernel minus mixture
#################################################################################
data.plot.line <- NULL
pb <- txtProgressBar(max = length(files), style = 3) # progress bar
for (i in seq_along(files))
{
  load(files[i])

  for (j in seq_along(rvs))
  {
    obs <- data.hist$x[which(data.hist$quantity == rvs[j])]
    kde <- density(obs) # kernel density fits
    for (k in seq_along(mixtures))
    {
      tmp <- data.plot %>%
        filter(quantity == rvs[j]) %>%
        filter(group == mixtures[k])
      threshold <- pull(tmp, "x")
      density.mix <- pull(tmp, "y")
      density.kde <- approx(kde$x, kde$y, xout = threshold)$y
      density.diff <- density.kde - density.mix

      data.plot.line <- data.plot.line %>%
        bind_rows(., tibble(x = threshold, y = density.diff, station = substr(files[i], 1, 3), mixture = mixtures[k], index = rvs[j]))
    }
  }

  setTxtProgressBar(pb, i)
} # end i
close(pb)

data.plot.line$mixture <- factor(data.plot.line$mixture, levels = c("2-component N", "3-component N", "3-component TSN"), labels = c("2-comp~N", "3-comp~N", "3-comp~TSN"))
data.plot.line$index <- factor(data.plot.line$index, levels = c("kappa", "kt", "kb", "kd"), labels = c("Clear~sky~index*\', \'*~kappa", "Clearness~index*\', \'*~italic(k)[italic(t)]", "Beam~transmittance*\', \'*~italic(k)[italic(b)]", "Diffuse~transmittance*\', \'*~italic(k)[italic(d)]"))

p3 <- ggplot() +
  geom_line(aes(x = x, y = y, group = station), data = data.plot.line, linewidth = line.size, alpha = 0.5) +
  facet_wrap(~ mixture + index, labeller = label_parsed, ncol = 4, scales = "free_x") +
  scale_x_continuous(name = expression(paste("Clearness index, ", italic(k)[italic(t)], " [dimensionless]")), limits = c(0, 1.2), breaks = c(0, 0.5, 1)) +
  scale_y_continuous(name = expression(paste("KDE [dimensionless]")), limits = c(-10, 10)) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2, 0.2, 0, 0.4), "lines"), text = element_text(family = "Times", size = plot.size), axis.text = element_text(size = plot.size), panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "gray80", linewidth = line.size / 2), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), strip.text.x = element_text(size = plot.size, margin = unit(c(0.1, 0, 0.1, 0), "lines")), strip.text.y = element_text(size = plot.size, margin = unit(c(0, 0.1, 0, 0.1), "lines")), panel.spacing = unit(0, "lines"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = plot.size), legend.background = element_rect(fill = "transparent"), legend.key = element_rect(fill = "transparent"), legend.box.margin = unit(c(-0.7, 0, 0, 0.0), "lines"), legend.direction = "horizontal") # legend.position = "inside", legend.position.inside = c(0.3, 0.7),

p3
