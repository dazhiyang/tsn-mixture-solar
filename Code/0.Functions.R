################################################################################
# This code is written by Dazhi Yang
# School of Electrical Engineering and Automation
# Harbin Institute of Technology
# emails: yangdazhi.nus@gmail.com
################################################################################

# Density function for truncated skew-normal
dtsn <- function(x, xi, omega, alpha, a, b, log = FALSE) {
  # Untruncated skew-normal density
  sn_dens <- sn::dsn(x, xi = xi, omega = omega, alpha = alpha)

  # Truncation constant
  cdf_a <- sn::psn(a, xi = xi, omega = omega, alpha = alpha)
  cdf_b <- sn::psn(b, xi = xi, omega = omega, alpha = alpha)
  trunc_const <- cdf_b - cdf_a

  # Truncated density
  dens <- sn_dens / trunc_const
  dens[x < a | x > b] <- 0

  if (log) {
    return(log(dens))
  } else {
    return(dens)
  }
}

# distribution function for truncated skew-normal
ptsn <- function(q, xi, omega, alpha, a, b, lower.tail = TRUE, log.p = FALSE) {
  # Untruncated skew-normal CDF
  sn_cdf <- sn::psn(q, xi = xi, omega = omega, alpha = alpha)

  # Truncation constants
  cdf_a <- sn::psn(a, xi = xi, omega = omega, alpha = alpha)
  cdf_b <- sn::psn(b, xi = xi, omega = omega, alpha = alpha)
  trunc_constant <- cdf_b - cdf_a

  # Truncated CDF
  cdf <- (sn_cdf - cdf_a) / trunc_constant

  # Apply bounds
  cdf[q < a] <- 0
  cdf[q > b] <- 1

  # Handle lower.tail option
  if (!lower.tail) {
    cdf <- 1 - cdf
  }

  # Handle log.p option
  if (log.p) {
    cdf <- log(cdf)
  }

  return(cdf)
}

# random number generator for truncated skew-normal (use the inverse transform method)
rtsn <- function(n, xi, omega, alpha, a, b) {
  F_a <- sn::psn(a, xi = xi, omega = omega, alpha = alpha)
  F_b <- sn::psn(b, xi = xi, omega = omega, alpha = alpha)
  u <- stats::runif(n, F_a, F_b)
  result <- tryCatch(
    {
      # Code that may produce an error (the "try" block)
      sn::qsn(u, xi = xi, omega = omega, alpha = alpha)
    },
    error = function(e) {
      # Code to execute if an error occurs (the "catch" block)
      cat("An error occurred:", conditionMessage(e), "\n")
      # You can return a default value or switch to another part of the code here
      sn::qsn(u, xi = xi, omega = omega, alpha = alpha, solver = "RFB")
    }
  )
  return(result)
}

# Density function for three-component truncated skew-normal mixture
dtsn3 <- function(x, a, b, fit) {
  # weighted component distributions
  cp1 <- fit$p[1] * dtsn(x, xi = fit$xi[1], omega = fit$omega[1], alpha = fit$alpha[1], a = a, b = b)
  cp2 <- fit$p[2] * dtsn(x, xi = fit$xi[2], omega = fit$omega[2], alpha = fit$alpha[2], a = a, b = b)
  cp3 <- fit$p[3] * dtsn(x, xi = fit$xi[3], omega = fit$omega[3], alpha = fit$alpha[3], a = a, b = b)

  dens <- cp1 + cp2 + cp3
  return(dens)
}

# Three-component truncated skew-normal mixture
tsn3_loglik <- function(logit_p1, logit_p2,
                        xi1, log_omega1, alpha1,
                        xi2, log_omega2, alpha2,
                        xi3, log_omega3, alpha3,
                        x, a, b) {
  # Transform parameters
  p1 <- stats::plogis(logit_p1)
  p2 <- stats::plogis(logit_p2)
  omega1 <- exp(log_omega1)
  omega2 <- exp(log_omega2)
  omega3 <- exp(log_omega3)

  # Ensure probabilities sum to 1
  p3 <- 1 - p1 - p2
  if (p3 < 0.001 || p3 > 0.999) {
    return(1e10) # Penalize invalid probabilities
  }

  # Add heavy penalty if distribution is outside data bounds
  penalty <- 0
  if (xi1 + 3 * omega1 < a || xi1 - 3 * omega2 > b) {
    penalty <- penalty + 1000 # Large penalty for unreasonable parameters
  } else if (xi2 + 3 * omega2 < a || xi2 - 3 * omega2 > b) {
    penalty <- penalty + 1000 # Large penalty for unreasonable parameters
  } else if (xi3 + 3 * omega3 < a || xi3 - 3 * omega3 > b) {
    penalty <- penalty + 1000 # Large penalty for unreasonable parameters
  }

  # Calculate component densities
  tsn1_dens <- dtsn(x, xi1, omega1, alpha1, a, b)
  tsn2_dens <- dtsn(x, xi2, omega2, alpha2, a, b)
  tsn3_dens <- dtsn(x, xi3, omega3, alpha3, a, b)

  # Mixture density
  mixture_dens <- p1 * tsn1_dens + p2 * tsn2_dens + p3 * tsn3_dens

  # Handle very small values
  mixture_dens <- pmax(mixture_dens, .Machine$double.eps)

  # Return negative log-likelihood
  -sum(log(mixture_dens)) + penalty
}

# Improved initial parameter estimation for three-component truncated skew-normal
get_tsn3_initial_params <- function(dat) {
  # Estimate mixing proportions using k-means clustering
  km <- stats::kmeans(dat, centers = 3, nstart = 25)
  props <- table(km$cluster) / length(dat)

  # Sort clusters by their means
  cluster_means <- tapply(dat, km$cluster, mean)
  sorted_clusters <- order(cluster_means)

  # Assign components based on cluster characteristics
  cluster1 <- sorted_clusters[1]
  cluster2 <- sorted_clusters[2]
  cluster3 <- sorted_clusters[3]

  cluster1_dat <- dat[km$cluster == cluster1]
  cluster2_dat <- dat[km$cluster == cluster2]
  cluster3_dat <- dat[km$cluster == cluster3]

  # Estimate skew-normal parameters using method of moments
  estimate_params <- function(x) {
    # Calculate sample moments
    m <- mean(x)
    s <- stats::sd(x)
    skew <- if (s > 0) mean((x - m)^3) / s^3 else 0

    # Solve for delta from skewness (numerically)
    delta <- find_delta_from_skewness(skew)

    # Now back-calculate omega and xi
    omega <- s / sqrt(1 - (2 * delta^2) / pi)
    xi <- m - omega * delta * sqrt(2 / pi)

    # Convert delta to alpha
    alpha <- delta / sqrt(1 - delta^2)

    # Regularize extreme values
    alpha <- sign(alpha) * pmin(10, abs(alpha))

    list(xi = xi, omega = omega, alpha = alpha)
  }
  # estimate_params <- function(x) {
  #   m <- mean(x)
  #   s <- sd(x)
  #   skew <- if (s > 0) mean((x - m)^3) / s^3 else 0
  #   alpha <- sign(skew) * min(10, abs(skew) * 2)  # Regularized estimate
  #   list(xi = m, omega = s, alpha = alpha)
  # }

  # Estimate parameters for all three components
  tsn1_params <- estimate_params(cluster1_dat)
  tsn2_params <- estimate_params(cluster2_dat)
  tsn3_params <- estimate_params(cluster3_dat)

  return(list(
    logit_p1 = stats::qlogis(props[cluster1]),
    logit_p2 = stats::qlogis(props[cluster2]),
    xi1 = tsn1_params$xi,
    log_omega1 = log(pmax(0.05, tsn1_params$omega)),
    alpha1 = tsn1_params$alpha,
    xi2 = tsn2_params$xi,
    log_omega2 = log(pmax(0.05, tsn2_params$omega)),
    alpha2 = tsn2_params$alpha,
    xi3 = tsn3_params$xi,
    log_omega3 = log(pmax(0.05, tsn3_params$omega)),
    alpha3 = tsn3_params$alpha
  ))
}

# Helper function to solve for delta from skewness
find_delta_from_skewness <- function(skew, tol = 1e-6) {
  if (abs(skew) < tol) {
    return(0)
  }

  # Numerical solution for delta
  delta_vals <- seq(-0.99, 0.99, by = 0.001)
  theoretical_skew <- sapply(delta_vals, function(delta) {
    (4 - pi) / 2 * (delta * sqrt(2 / pi))^3 / (1 - 2 * delta^2 / pi)^(3 / 2)
  })

  # Find delta that minimizes difference
  delta_vals[which.min(abs(theoretical_skew - skew))]
}

# Calculate PIT residuals for normal mixtures
pit_residuals_n <- function(data, fit) {
  cdf_values <- numeric(length(data))
  for (i in seq_along(data)) {
    cdf_values[i] <- sum(fit$lambda * pnorm(data[i], fit$mu, fit$sigma))
  }
  return(cdf_values)
}

# Calculate PIT residuals for tsn mixtures
pit_residuals_tsn3 <- function(data, a, b, fit) {
  cdf_values <- numeric(length(data))
  for (i in seq_along(data)) {
    cdf_values[i] <- sum(fit$p * ptsn(data[i], fit$xi, fit$omega, fit$alpha, a = a, b = b))
  }

  # PIT residuals should be ~Uniform(0,1) if model is correct
  return(cdf_values)
}

# Calculate AIC and BIC
calculate_information_criteria <- function(loglik, n_params, n_obs) {
  aic <- 2 * n_params - 2 * loglik
  bic <- n_params * log(n_obs) - 2 * loglik
  # aicc <- aic + (2 * n_params * (n_params + 1)) / (n_obs - n_params - 1)

  return(data.frame(loglik = loglik, n_params = n_params, AIC = aic, BIC = bic))
}

# Vectorized way to generate normal mixture samples
generate_mixture_samples <- function(n, fit) {
  # Extract parameters
  lambda <- fit$lambda
  mu <- fit$mu
  sigma <- fit$sigma

  # Sample components for each observation
  comps <- sample(seq_along(lambda), n, replace = TRUE, prob = lambda)
  # Sample from the corresponding normal distributions
  rnorm(n, mu[comps], sigma[comps])
}

# Vectorized function to generate from truncated skew-normal mixture
generate_tsn_mixture <- function(n, a, b, fit) {
  # Extract parameters
  p <- fit$p
  xi <- fit$xi
  omega <- fit$omega
  alpha <- fit$alpha

  # Sample components for each observation
  comps <- sample(seq_along(p), n, replace = TRUE, prob = p)

  # Generate samples using vectorized approach
  samples <- numeric(n)
  for (j in 1:3) {
    comp_indices <- which(comps == j)
    if (length(comp_indices) > 0) {
      samples[comp_indices] <- rtsn(length(comp_indices), xi[j], omega[j], alpha[j], a, b)
    }
  }

  return(samples)
}

# a function to convert a numeric vector to latex-ready format to be used by xtable
latex_sci <- function(x, digits = 3) {
  sapply(x, function(num) {
    if (num == 0) {
      return("$0$")
    }
    if (abs(num) >= 1) {
      # For numbers >= 1, use regular formatting
      return(sprintf("$%.*f$", digits, num))
    }

    exponent <- floor(log10(abs(num)))
    coefficient <- num / (10^exponent)

    sprintf("$%.*f \\times 10^{%d}$", digits, coefficient, exponent)
  })
}

# a function to convert latex-ready format to be used by xtable to a numeric vector
latex_sci_to_numeric <- function(latex_str_vector) {
  # Vectorized conversion
  sapply(latex_str_vector, function(latex_str) {
    # Remove LaTeX formatting
    clean_str <- gsub("\\$|\\\\times|\\s", "", latex_str)

    # Split into coefficient and exponent parts
    parts <- strsplit(clean_str, "10\\^")[[1]]

    if (length(parts) == 1) {
      # No exponent part (regular number)
      return(as.numeric(parts[1]))
    }

    # Extract coefficient and exponent
    coefficient <- as.numeric(parts[1])

    # Handle exponent (may have braces)
    exponent_str <- parts[2]
    exponent_str <- gsub("\\{|\\(|\\]|\\)|\\}", "", exponent_str) # Remove brackets
    exponent <- as.numeric(exponent_str)

    # Calculate final value
    coefficient * 10^exponent
  })
}

# convert station coordinates to the continent it locates
coords2continent <- function(points) {
  countriesSP <- rworldmap::getMap(resolution = "low")
  # countriesSP <- getMap(resolution='high') #you could use high res map from rworldxtra if you were concerned about detail

  # converting points to a SpatialPoints object, setting CRS directly to that from rworldmap
  pointsSP <- sp::SpatialPoints(points, proj4string = sp::CRS(sp::proj4string(countriesSP)))

  # use 'over' to get indices of the Polygons object containing each point
  indices <- sp::over(pointsSP, countriesSP)

  # indices$continent   # returns the continent (6 continent model)
  indices$REGION # returns the continent (7 continent model)
  # indices$ADMIN  #returns country name
  # indices$ISO3 # returns the ISO3 code
}

# Function to calculate weighted Wasserstein distance between two stations
calculate_station_distance <- function(station1, station2, weights = NULL) {
  if (is.null(weights)) {
    weights <- c(1, 1, 1, 1) # Equal weights by default
  }

  distances <- numeric(4)
  trans_types <- c("kappa", "kt", "kb", "kd")

  for (i in seq_along(trans_types)) {
    trans_type <- trans_types[i]

    pdf1 <- station1[[trans_type]]
    pdf2 <- station2[[trans_type]]

    # Normalize PDFs to ensure they sum to 1
    pdf1_norm <- pdf1$y / sum(pdf1$y)
    pdf2_norm <- pdf2$y / sum(pdf2$y)

    # Calculate cumulative distribution functions
    cdf1 <- cumsum(pdf1_norm)
    cdf2 <- cumsum(pdf2_norm)

    # Wasserstein-1 distance is the integral of |CDF1^{-1} - CDF2^{-1}|
    # For 1D, it's equivalent to the integral of |CDF1(x) - CDF2(x)| dx
    wasserstein_dist <- sum(abs(cdf1 - cdf2)) * diff(pdf1$x[1:2])

    distances[i] <- wasserstein_dist * weights[i]
  }

  # Combine distances (you can use sum, weighted sum, or other combination)
  total_distance <- sum(distances) # Simple sum
  # total_distance <- sqrt(sum(distances^2))  # Euclidean combination
  # total_distance <- max(distances)  # Worst-case distance

  return(total_distance)
}
