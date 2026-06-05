#' HBomics algorithm: calling clonal allele imbalance events from multiomics data with a hierarchical Bayesian framework
#' Author: Zhen Y
#' Utils function
#' Calculate the allele imbalance ratio (the summit) from the AI posterior distribution
#' @param x A numeric vector of MCMC samples
#' @return The inferred allele imbalance ratio
#' @importFrom stats density
#' @export
Mode2 <- function(x) {
  if (!is.numeric(x) || length(x) < 2 || anyNA(x)) {
    abort_hbomics("Error: 'x' must be a numeric vector with at least two non-missing values.")
  }
  ux <- density(x)
  mode_val <- ux$x[which.max(ux$y)]
  return(mode_val)
}

#' Get Confidence/Credible Interval boundaries of the AI posterior distribution
#' @param x A numeric vector of MCMC samples
#' @param bound The interval boundary (default is 0.95)
#' @return A vector with the lower and upper quantiles
#' @importFrom stats quantile
#' @export
CI_quant <- function(x, bound = 0.95) {
  if (!is.numeric(x) || length(x) < 1 || anyNA(x)) {
    abort_hbomics("Error: 'x' must be a numeric vector with no missing values.")
  }
  validate_bound(bound)
  b1 <- 1 - bound
  b2 <- bound
  return(quantile(x, c(b1, b2)))
}

#' Compute KL Divergence between two MCMC posterior samples (sample-wise comparison)
#' @param p_samples Numeric vector of MCMC samples from distribution P
#' @param q_samples Numeric vector of MCMC samples from distribution Q
#' @param n_grid Number of points to evaluate the density over (default 512)
#' @return Numeric KL divergence value
#' @importFrom stats density
#' @export
calc_kl_divergence <- function(c_samples, t_samples, n_grid = 512) {
  if (!is.numeric(c_samples) || length(c_samples) < 2 || anyNA(c_samples)) {
    abort_hbomics("Error: 'c_samples' must be a numeric vector with at least two non-missing values.")
  }
  if (!is.numeric(t_samples) || length(t_samples) < 2 || anyNA(t_samples)) {
    abort_hbomics("Error: 't_samples' must be a numeric vector with at least two non-missing values.")
  }
  validate_integerish_scalar(n_grid, "n_grid", lower = 2)

  c_dens <- density(c_samples, from = 0, to = 1, n = n_grid)$y
  t_dens <- density(t_samples, from = 0, to = 1, n = n_grid)$y

  # Add a tiny epsilon to prevent log(0) or division by zero errors
  eps <- 1e-10
  c_dens <- c_dens + eps
  t_dens <- t_dens + eps

  # Normalize so the area under the curve equals 1
  c_dens <- c_dens / sum(c_dens)
  t_dens <- t_dens / sum(t_dens)

  # Calculate standard KL divergence: sum(P * log(P/Q))
  kl_div <- sum(c_dens * log(c_dens / t_dens))

  return(kl_div)
}