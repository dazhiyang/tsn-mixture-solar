################################################################################
# This code is written by Dazhi Yang with Antigravity
# School of Electrical Engineering and Automation
# Harbin Institute of Technology
# emails: yangdazhi.nus@gmail.com
################################################################################

rm(list = ls(all = TRUE))
libs <- c("dplyr", "lubridate", "mixtools", "sn", "bbmle", "ggplot2", "xtable", "transport")
invisible(lapply(libs, library, character.only = TRUE))

# Helper for LaTeX-style negative numbers
latex_neg <- function(x, digits = 2) {
    sapply(x, function(val) {
        if (is.na(val)) {
            return("NA")
        }
        formatted <- sprintf("%.*f", digits, val)
        if (val < 0) {
            return(paste0("$", formatted, "$"))
        } else {
            return(formatted)
        }
    })
}

################################################################################
# global input
################################################################################
dir0 <- "/Users/seryangd/Library/CloudStorage/Dropbox/Working papers/distribution"
dir.data <- "/Volumes/Macintosh Research/Data/Separation/tibbles"
dir.sens <- file.path(dir0, "Data/Sensitivity")
if (!dir.exists(dir.sens)) dir.create(dir.sens, recursive = TRUE)

# source the necessary functions for fitting mixture models
source(file.path(dir0, "Code/0.Functions.R"))

# get the variables of interest
rvs <- c("kappa", "kt", "kb", "kd")
n_samples_test <- c(10000, 30000, 50000)
################################################################################

# read the metadata of sites
setwd(file.path(dir0, "Data"))
locations <- read.csv("location.csv", header = TRUE)

# Define the stations for sensitivity analysis
# Feel free to change these names to test other locations
stn_ids <- c("BON", "TAM")

cat("Selected locations for sensitivity analysis:", paste(stn_ids, collapse = " and "), "\n")

results_all <- NULL

for (i in seq_along(stn_ids)) {
    stn_id <- stn_ids[i]
    file_path <- file.path(dir.data, paste0(stn_id, ".RData"))

    if (!file.exists(file_path)) {
        cat("File not found for", stn_id, "- skipping.\n")
        next
    }

    # Load data
    load(file_path) # This loads 'data'

    # Data preprocessing (from MainLoop.R)
    data <- data %>%
        rename(kappa = kc) %>%
        mutate(Bh = pmax(Gh - Dh, 0)) %>%
        filter(Bh > 2) %>%
        mutate(kb = Bh / E0, kd = Dh / E0)

    for (n_size in n_samples_test) {
        cat("Processing", stn_id, "with sample size", n_size, "...\n")

        set.seed(123) # consistent sampling for each size comparison
        mixtsn3.all.sens <- list()

        for (j in seq_along(rvs)) {
            rv <- rvs[j]
            dat_full <- pull(data, rv)

            # Ensure we don't ask for more than available
            actual_n <- min(n_size, length(dat_full))
            dat <- sample(dat_full, size = actual_n, replace = FALSE)

            # Bounds
            lower_lim <- 0
            upper_lim <- switch(rv,
                "kappa" = 1.6,
                "kt" = 1.2,
                "kb" = 0.9,
                "kd" = 0.7
            )

            # Fit 3-component TSN mixture (only testing the proposed model for sensitivity)
            params <- get_tsn3_initial_params(dat)

            # Stage 1: Fast coarse optimization
            fit_coarse <- tryCatch(
                {
                    mle2(tsn3_loglik,
                        data = list(x = dat, a = lower_lim, b = upper_lim),
                        start = params,
                        method = "Nelder-Mead",
                        control = list(maxit = 400, trace = 0)
                    )
                },
                error = function(e) {
                    return(NULL)
                }
            )

            if (is.null(fit_coarse)) next

            # Stage 2: Refine
            fit_final <- tryCatch(
                {
                    if (fit_coarse@details$convergence == 1) {
                        mle2(tsn3_loglik,
                            data = list(x = dat, a = lower_lim, b = upper_lim),
                            start = as.list(coef(fit_coarse)),
                            method = "BFGS",
                            control = list(maxit = 200, trace = 0)
                        )
                    } else {
                        fit_final <- fit_coarse
                    }
                },
                error = function(e) {
                    return(fit_coarse)
                }
            )

            # Parameters
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
            ranking <- order(mixtsn3$xi)
            mixtsn3 <- lapply(mixtsn3, function(x) {
                x[ranking]
            })
            mixtsn3.all.sens[[j]] <- mixtsn3

            # Performance metrics
            loglik_val <- as.numeric(logLik(fit_final))
            perf <- calculate_information_criteria(loglik = loglik_val, n_params = 11, n_obs = length(dat))

            mixtsn3_sample <- generate_tsn_mixture(n = length(dat), a = lower_lim, b = upper_lim, fit = mixtsn3)
            was_dist <- transport::wasserstein1d(dat, mixtsn3_sample)
            ks_stat <- ks.test(dat, mixtsn3_sample)$statistic

            results_all <- results_all %>%
                bind_rows(tibble(
                    Station = stn_id,
                    Size = n_size,
                    Var = rv,
                    p1 = mixtsn3$p[1], xi1 = mixtsn3$xi[1], omega1 = mixtsn3$omega[1], alpha1 = mixtsn3$alpha[1],
                    p2 = mixtsn3$p[2], xi2 = mixtsn3$xi[2], omega2 = mixtsn3$omega[2], alpha2 = mixtsn3$alpha[2],
                    p3 = mixtsn3$p[3], xi3 = mixtsn3$xi[3], omega3 = mixtsn3$omega[3], alpha3 = mixtsn3$alpha[3]
                ))
        }
        # Save results for this station and sample size
        save(mixtsn3.all.sens, file = file.path(dir.sens, paste0(stn_id, "_", n_size, ".RData")))
    }
}

# Format results for LaTeX
results_formatted <- results_all %>%
    mutate(
        p1 = latex_neg(p1), xi1 = latex_neg(xi1), omega1 = latex_neg(omega1), alpha1 = latex_neg(alpha1),
        p2 = latex_neg(p2), xi2 = latex_neg(xi2), omega2 = latex_neg(omega2), alpha2 = latex_neg(alpha2),
        p3 = latex_neg(p3), xi3 = latex_neg(xi3), omega3 = latex_neg(omega3), alpha3 = latex_neg(alpha3),
        size = as.integer(Size),
        var_order = match(Var, rvs),
        variable = case_when(
            Var == "kappa" ~ "$\\kappa$",
            Var == "kt" ~ "$k_t$",
            Var == "kb" ~ "$k_b$",
            Var == "kd" ~ "$k_d$",
            TRUE ~ Var
        )
    ) %>%
    arrange(Station, var_order, size) %>%
    select(Station, size, variable, p1, xi1, omega1, alpha1, p2, xi2, omega2, alpha2, p3, xi3, omega3, alpha3)

# Create LaTeX table
print(
    xtable(results_formatted,
        caption = "Fitted parameters of 3-component TSN mixture model across different sample sizes.",
        label = "tab:sensitivity_params"
    ),
    include.rownames = FALSE,
    sanitize.text.function = identity,
    file = file.path(dir0, "tex/sensitivity_table.tex")
)

cat("Sensitivity analysis complete. Table saved to tex/sensitivity_table.tex\n")
