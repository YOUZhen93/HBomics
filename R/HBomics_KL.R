#' HBomics algorithm: calling clonal allele imbalance events from multiomics data with a hierarchical Bayesian framework
#' Author: Zhen Y
#' @title KL divergence of sample-wise comparison 
#' @description This function uses posterior samples from HBomics object for site-by-site comparison of different samples
#' It also provides feature-by-feature comparison with allele depth weighted KL scores after
#' running allele_integ function
#' @param obj_C First fitted HBomics result object (e.g., control or primary)
#' @param obj_T Second fitted HBomics result object (e.g., treatment)
#' @param symmetric Logical; if TRUE, calculates Jeffrey's divergence (default TRUE)
#' @return A data frame with class `HBomicsKLResult`, sorted by descending KL
#'   divergence.
#' @export
compare_hbomics_kl <- function(obj_C, obj_T, symmetric = TRUE) {
  obj_C <- coerce_hbomics_result(obj_C, require_posterior = TRUE, arg_name = "obj_C")
  obj_T <- coerce_hbomics_result(obj_T, require_posterior = TRUE, arg_name = "obj_T")
  validate_flag(symmetric, "symmetric")

  mut_ids_C <- obj_C$hb_df$mutation_id
  mut_ids_T <- obj_T$hb_df$mutation_id

  # Find overlapping mutations
  common_muts <- intersect(mut_ids_C, mut_ids_T)
  intersect_muts_num <- length(common_muts)
  if (intersect_muts_num == 0) {
    abort_hbomics("Error: No overlapping mutation IDs found between the two objects.")
  }
  if (intersect_muts_num < max(length(mut_ids_C), length(mut_ids_T))) {
    warning(
      "WARNING: Not all mutations overlap. ",
      intersect_muts_num,
      " mutations will be used for the KL analysis.",
      call. = FALSE
    )
  }

  # getting the posterior samples
  index_C <- match(common_muts, mut_ids_C)
  index_T <- match(common_muts, mut_ids_T)
  post_C <- obj_C$post_a[, index_C, drop = FALSE]
  post_T <- obj_T$post_a[, index_T, drop = FALSE]

  kl_scores <- numeric(intersect_muts_num)

  # compute KL divergence
  message(sprintf("%s     INFO: calculating KL Divergence for %d shared het. SNVs ...", date(), intersect_muts_num))

  for (i in seq_len(intersect_muts_num)) {
    p_samples <- post_C[, i]
    q_samples <- post_T[, i]

    kl_pq <- calc_kl_divergence(p_samples, q_samples)

    if (symmetric) {
      kl_qp <- calc_kl_divergence(q_samples, p_samples)
      kl_scores[i] <- kl_pq + kl_qp
    } else {
      kl_scores[i] <- kl_pq
    }
  }

  # format and return results
  res_df <- data.frame(
    mutation_id = common_muts,
    KL_divergence = kl_scores,
    stringsAsFactors = FALSE
  )

  # Sort descending to bring the largest subclonal shifts to the top
  res_df <- res_df[order(-res_df$KL_divergence), ]
  rownames(res_df) <- NULL

  new_hbomics_kl_result(res_df)
}

