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
# dir0 <- "/Users/dyang/Dropbox/Working papers/distribution"
dir0 <- "/Users/seryangd/Library/CloudStorage/Dropbox/Working papers/distribution"
dir.data <- "/Volumes/Macintosh Research/Data/Separation/tibbles"

# source the necessary functions for fitting mixture models
source(file.path(dir0, "Code/0.Functions.R"))

# get the variables of interest
rvs <- c("kappa", "kt", "kb", "kd")
################################################################################

# read the metadata of sites
setwd(file.path(dir0, "Data"))
locations <- read.csv("location.csv", header = TRUE)

# access the folder contained arranged data from 126 stations
setwd(dir.data)
files <- dir()

# the loop to fit the mixture distribution and save results
set.seed(123) # for reproducibility
for (i in seq_along(files))
{
  # load data from a site
  load(files[i])

  # manipulate the data and get kappa, kt, kb, and kd
  data <- data %>%
    rename(kappa = kc) %>%
    mutate(Bh = pmax(Gh - Dh, 0)) %>%
    filter(Bh > 2) %>%
    mutate(kb = Bh / E0, kd = Dh / E0)

  # construct data for saving
  data.pit <- data.perf <- data.hist <- data.plot <- NULL
  mixn2.all <- mixn3.all <- mixtsn3.all <- list()
  for (j in seq_along(rvs))
  {
    rv <- rvs[j]
    dat <- pull(data, rv)
    dat <- sample(dat, size = min(30000, nrow(data)), replace = FALSE)

    # two-component Gaussian mixture
    mixn2 <- normalmixEM(dat, mu = c(0.2, 0.7), sigma = c(0.1, 0.1), k = 2)
    mixn2.all[[j]] <- mixn2
    # three-component Gaussian mixture
    mixn3 <- normalmixEM(dat, mu = c(0.2, 0.4, 0.7), sigma = c(0.1, 0.2, 0.1), k = 3)
    mixn3.all[[j]] <- mixn3

    # two-component truncated skew-normal mixture
    # get the bound of random variable of interest
    lower_lim <- 0
    upper_lim <- switch(rv,
      "kappa" = 1.6,
      "kt" = 1.2,
      "kb" = 0.9,
      "kd" = 0.7
    )

    params <- get_tsn3_initial_params(dat)
    # Stage 1: Fast coarse optimization
    fit_coarse <- mle2(tsn3_loglik,
      data = list(x = dat, a = lower_lim, b = upper_lim), # Provide fixed values here
      start = params,
      method = "Nelder-Mead",
      control = list(maxit = 400, trace = 1)
    )
    # Stage 2: Refine with better method
    if (fit_coarse@details$convergence == 1) {
      fit_final <- mle2(tsn3_loglik,
        data = list(x = dat, a = lower_lim, b = upper_lim), # Provide fixed values here
        start = as.list(coef(fit_coarse)),
        method = "BFGS", # Switch to faster method
        control = list(maxit = 200, trace = 1)
      )
    } else {
      fit_final <- fit_coarse
    }

    # combine the parameters into a list similar to "mixEM" class
    mixtsn3 <- list(
      p = c(
        plogis(as.numeric(coef(fit_final)["logit_p1"])),
        plogis(as.numeric(coef(fit_final)["logit_p2"])),
        1 - plogis(as.numeric(coef(fit_final)["logit_p1"])) - plogis(as.numeric(coef(fit_final)["logit_p2"]))
      ),
      xi = c(
        as.numeric(coef(fit_final)["xi1"]),
        as.numeric(coef(fit_final)["xi2"]),
        as.numeric(coef(fit_final)["xi3"])
      ),
      omega = c(
        exp(as.numeric(coef(fit_final)["log_omega1"])),
        exp(as.numeric(coef(fit_final)["log_omega2"])),
        exp(as.numeric(coef(fit_final)["log_omega3"]))
      ),
      alpha = c(
        as.numeric(coef(fit_final)["alpha1"]),
        as.numeric(coef(fit_final)["alpha2"]),
        as.numeric(coef(fit_final)["alpha3"])
      )
    )

    # rank the components from smallest to largest location
    ranking <- order(mixtsn3$xi)
    mixtsn3 <- lapply(mixtsn3, function(x) {
      x[ranking]
    })

    # save the estimated parameters
    mixtsn3.all[[j]] <- mixtsn3

    threshold <- seq(lower_lim - 0.1, upper_lim + 0.1, 0.005)
    data.plot <- data.plot %>%
      bind_rows(., tibble(x = threshold, y = sapply(threshold, function(x) {
        sum(mixn2$lambda * dnorm(x, mixn2$mu, mixn2$sigma))
      }), group = "2-component N", quantity = rv)) %>%
      bind_rows(., tibble(x = threshold, y = sapply(threshold, function(x) {
        sum(mixn3$lambda * dnorm(x, mixn3$mu, mixn3$sigma))
      }), group = "3-component N", quantity = rv)) %>%
      bind_rows(., tibble(x = threshold, y = dtsn3(threshold, a = lower_lim, b = upper_lim, fit = mixtsn3), group = "3-component TSN", quantity = rv))

    data.hist <- data.hist %>%
      bind_rows(., tibble(x = dat, quantity = rv))

    # compute pit and save into tibble
    data.pit <- data.pit %>%
      bind_rows(., tibble(pit = pit_residuals_n(dat, mixn2), group = "2-component N", quantity = rv)) %>%
      bind_rows(., tibble(pit = pit_residuals_n(dat, mixn3), group = "3-component N", quantity = rv)) %>%
      bind_rows(., tibble(pit = pit_residuals_tsn3(dat, a = lower_lim, b = upper_lim, fit = mixtsn3), group = "3-component TSN", quantity = rv))

    # log-likelihood, information criteria, MSE, Wasserstein,
    # ICs
    perf_n2 <- as.character(format(round(calculate_information_criteria(loglik = mixn2$loglik, n_params = 5, n_obs = length(dat)), 2), nsmall = 0))
    perf_n3 <- as.character(format(round(calculate_information_criteria(loglik = mixn3$loglik, n_params = 8, n_obs = length(dat)), 2), nsmall = 0))
    perf_tsn3 <- as.character(format(round(calculate_information_criteria(loglik = logLik(fit_final), n_params = 11, n_obs = length(dat)), 2), nsmall = 0))

    # Wasserstein distance
    mixn2_sample <- generate_mixture_samples(n = length(dat), fit = mixn2)
    mixn3_sample <- generate_mixture_samples(n = length(dat), fit = mixn3)
    mixtsn3_sample <- generate_tsn_mixture(n = length(dat), a = lower_lim, b = upper_lim, fit = mixtsn3)
    wasserstein <- c(
      transport::wasserstein1d(dat, mixn2_sample),
      transport::wasserstein1d(dat, mixn3_sample),
      transport::wasserstein1d(dat, mixtsn3_sample)
    )
    wasserstein <- latex_sci(wasserstein)

    # KS tests
    ks <- c(
      ks.test(dat, mixn2_sample)$statistic,
      ks.test(dat, mixn3_sample)$statistic,
      ks.test(dat, mixtsn3_sample)$statistic
    )

    var_latex <- switch(rv,
      "kappa" = "$\\kappa$",
      "kt" = "$k_t$",
      "kb" = "$k_b$",
      "kd" = "$k_d$"
    )
    data.perf <- data.perf %>%
      bind_rows(., tibble(quantity = "", group = "2-component N", loglik = perf_n2[1], npar = perf_n2[2], aic = perf_n2[3], bic = perf_n2[4], was = wasserstein[1], ksd = ks[1])) %>%
      bind_rows(., tibble(quantity = var_latex, group = "3-component N", loglik = perf_n3[1], npar = perf_n3[2], aic = perf_n3[3], bic = perf_n3[4], was = wasserstein[2], ksd = ks[2])) %>%
      bind_rows(., tibble(quantity = "", group = "3-component TSN", loglik = perf_tsn3[1], npar = perf_tsn3[2], aic = perf_tsn3[3], bic = perf_tsn3[4], was = wasserstein[3], ksd = ks[3]))
  }

  # save data
  old.dir <- getwd()
  setwd(file.path(dir0, "Data/Results"))
  save(list = c("data.pit", "data.plot", "data.hist", "data.perf", "mixn2.all", "mixn3.all", "mixtsn3.all"), file = paste0(locations$stn[i], ".RData"))


  # plotting and tables
  data.plot$group <- factor(data.plot$group, levels = c("2-component N", "3-component N", "3-component TSN"))
  data.pit$group <- factor(data.pit$group, levels = c("2-component N", "3-component N", "3-component TSN"), labels = c("2-comp~N", "3-comp~N", "3-comp~TSN"))
  data.plot$quantity <- factor(data.plot$quantity, levels = c("kappa", "kt", "kb", "kd"), labels = c("Clear~sky~index*\', \'*~kappa", "Clearness~index*\', \'*~italic(k)[italic(t)]", "Beam~transmittance*\', \'*~italic(k)[italic(b)]", "Diffuse~transmittance*\', \'*~italic(k)[italic(d)]"))
  data.hist$quantity <- factor(data.hist$quantity, levels = c("kappa", "kt", "kb", "kd"), labels = c("Clear~sky~index*\', \'*~kappa", "Clearness~index*\', \'*~italic(k)[italic(t)]", "Beam~transmittance*\', \'*~italic(k)[italic(b)]", "Diffuse~transmittance*\', \'*~italic(k)[italic(d)]"))
  data.pit$quantity <- factor(data.pit$quantity, levels = c("kappa", "kt", "kb", "kd"), labels = c("Clear~sky~index*\', \'*~kappa", "Clearness~index*\', \'*~italic(k)[italic(t)]", "Beam~transmittance*\', \'*~italic(k)[italic(b)]", "Diffuse~transmittance*\', \'*~italic(k)[italic(d)]"))

  # distribution fit
  p1 <- ggplot() +
    geom_histogram(aes(x = x, y = ..density..), data = data.hist, bins = 60, fill = "grey80", color = "gray30", linewidth = line.size) +
    geom_line(aes(x = x, y = y, linetype = group), data = data.plot, linewidth = line.size * 2, alpha = 0.8) +
    scale_linetype_manual(values = c("dashed", "dotdash", "solid")) +
    # scale_color_manual(values = colorblind_pal()(8)[2:4]) +
    facet_wrap(~quantity, labeller = label_parsed, nrow = 2, scales = "free") +
    scale_x_continuous(name = expression(paste(italic(k), "-index [dimensionless]"))) +
    scale_y_continuous(name = expression(paste("Prob. density [dimensionless]"))) +
    theme_bw() +
    theme(plot.margin = unit(c(0.2, 0.2, 0, 0.4), "lines"), text = element_text(family = "Times", size = plot.size), axis.text = element_text(size = plot.size), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), strip.text.x = element_text(size = plot.size, margin = unit(c(0.1, 0, 0.1, 0), "lines")), strip.text.y = element_text(size = plot.size, margin = unit(c(0, 0.1, 0, 0.1), "lines")), panel.spacing = unit(0, "lines"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = plot.size), legend.background = element_rect(fill = "transparent"), legend.key = element_rect(fill = "transparent"), legend.box.margin = unit(c(-0.7, 0, 0, 0.0), "lines"), legend.direction = "horizontal") # legend.position = "inside", legend.position.inside = c(0.3, 0.7),

  setwd(file.path(dir0, "Data/FigDist"))
  ggsave(filename = paste0(locations$stn[i], "dist.pdf"), plot = p1, width = 160, height = 90, unit = "mm")

  # pit histogram
  p2 <- ggplot() +
    geom_histogram(aes(x = pit, y = ..density..), data = data.pit, bins = 40, fill = "grey80", color = "gray30", linewidth = line.size) +
    facet_grid(group ~ quantity, labeller = label_parsed) +
    scale_x_continuous(name = expression(paste(italic(k), "-index [dimensionless]")), breaks = seq(0, 1, by = 0.25), labels = c("0", "0.25", "0.5", "0.75", "1")) +
    scale_y_continuous(name = expression(paste("Prob. density [dimensionless]"))) +
    theme_bw() +
    theme(plot.margin = unit(c(0.2, 0.2, 0, 0.4), "lines"), text = element_text(family = "Times", size = plot.size), axis.text = element_text(size = plot.size), panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "gray80", linewidth = line.size / 2), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), strip.text.x = element_text(size = plot.size, margin = unit(c(0.1, 0, 0.1, 0), "lines")), strip.text.y = element_text(size = plot.size, margin = unit(c(0, 0.1, 0, 0.1), "lines")), panel.spacing = unit(0, "lines"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = plot.size), legend.background = element_rect(fill = "transparent"), legend.key = element_rect(fill = "transparent"), legend.box.margin = unit(c(-0.7, 0, 0, 0.0), "lines"), legend.direction = "horizontal") # legend.position = "inside", legend.position.inside = c(0.3, 0.7),

  setwd(file.path(dir0, "Data/FigPIT"))
  ggsave(filename = paste0(locations$stn[i], "pit.pdf"), plot = p2, width = 160, height = 90, unit = "mm")

  setwd(old.dir)
}
