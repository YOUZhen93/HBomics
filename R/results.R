#' HBomics algorithm: calling clonal allele imbalance events from multiomics data with a hierarchical Bayesian framework
#' Author: Zhen Y
#'
#' @description Result functions for the HBomics package
#' @keywords internal
NULL
#
new_hbomics_result <- function(hb_df, diagnostics, post_a = NULL) {
  structure(
    list(
      hb_df = hb_df,
      diagnostics = diagnostics,
      post_a = post_a
    ),
    class = "HBomicsResult"
  )
}

# adding allele imbalance to feature levels
new_hbomics_allele_result <- function(gene_allele_imb, HBomics) {
  structure(
    list(
      gene_allele_imb = gene_allele_imb,
      HBomics = HBomics
    ),
    class = "HBomicsAlleleResult"
  )
}

new_hbomics_kl_result <- function(x) {
  class(x) <- c("HBomicsKLResult", class(x))
  x
}

extract_convergence_diagnostics <- function(fit, mutation_id) {
  message(date(), "     INFO: Extracting convergence diagnostics (R-hat and ESS)...")

  # Extract global max R-hat and min ESS across ALL parameters
  # including a, u, v, mu, and theta
  fit_summary <- rstan::summary(fit)$summary
  params <- rownames(fit_summary)
  valid_params <- params[params != "lp__"]

  max_rhat <- max(fit_summary[valid_params, "Rhat"], na.rm = TRUE)
  min_ess <- min(fit_summary[valid_params, "n_eff"], na.rm = TRUE)
  
  # Extract metrics for selection strength on allele imbalance a
  a_idx <- grep("^a\\[", params)
  # tidy things up
  diag_df <- data.frame(
    mutation_id = mutation_id,
    Rhat = fit_summary[a_idx, "Rhat"],
    ESS = fit_summary[a_idx, "n_eff"]
  )

  message(sprintf("INFO: Convergence === Global Max R-hat: %.4f", max_rhat))
  message(sprintf("INFO: Convergence === Global Min ESS:   %.0f", min_ess))

  if (max_rhat > 1.05) {
    warning(
      "Some parameters have R-hat > 1.05. Chains may not have fully converged. ",
      "Consider increasing 'iter' or 'warmup'.",
      call. = FALSE
    )
  }

  list(
    global_max_rhat = max_rhat,
    global_min_ess = min_ess,
    a_diagnostics = diag_df
  )
}

extract_hbomics_summary <- function(fit, mutation_id, meta, bound = 0.95) {
  validate_bound(bound)
  # Process the allele imbalance posterior samples
  posterior <- rstan::extract(fit)
  posterior_a <- posterior$a
  denominator <- nrow(posterior_a)

  allele_imbalance <- apply(posterior_a, 2, Mode2)
  ci <- t(apply(posterior_a, 2, CI_quant, bound = bound))
  colnames(ci) <- c("low", "high")

  meta <- meta[match(mutation_id, meta$mutation_id), , drop = FALSE]
  extra_meta_cols <- setdiff(names(meta), "mutation_id")

  hb_df <- data.frame(
    mutation_id = mutation_id,
    low = ci[, "low"],
    high = ci[, "high"],
    mode = allele_imbalance,
    stringsAsFactors = FALSE
  )

  if (length(extra_meta_cols) > 0) {
    hb_df[extra_meta_cols] <- meta[extra_meta_cols]
  }
  # one-sided empirical p values
  p_val <- vapply(
    seq_along(allele_imbalance),
    function(i) {
      if (allele_imbalance[i] > 0.5) {
        sum(posterior_a[, i] < 0.5) / denominator
      } else {
        sum(posterior_a[, i] > 0.5) / denominator
      }
    },
    numeric(1)
  )

  hb_df$p_val <- p_val

  list(
    post_a = posterior_a,
    hb_df = hb_df
  )
}
