library(HBomics)

assert_true <- function(x, message) {
  if (!isTRUE(x)) {
    stop(message, call. = FALSE)
  }
}

fake_result <- structure(
  list(
    hb_df = data.frame(
      mutation_id = c("m1", "m2"),
      mode = c(0.60, 0.40),
      total_count = c(10, 20),
      alt_count = c(6, 8),
      gene_symbol = c("G1", "G1"),
      stringsAsFactors = FALSE
    ),
    diagnostics = list(
      global_max_rhat = 1.00,
      global_min_ess = 100,
      a_diagnostics = data.frame(
        mutation_id = c("m1", "m2"),
        Rhat = c(1.00, 1.00),
        ESS = c(100, 120)
      )
    ),
    post_a = matrix(
      c(0.60, 0.62,
        0.40, 0.38),
      nrow = 2
    )
  ),
  class = "HBomicsResult"
)

kl_res <- compare_hbomics_kl(fake_result, fake_result, symmetric = TRUE)
assert_true(inherits(kl_res, "HBomicsKLResult"), "compare_hbomics_kl() should return an HBomicsKLResult.")
assert_true(all(c("mutation_id", "KL_divergence") %in% names(kl_res)), "KL result should expose mutation_id and KL_divergence columns.")
assert_true(all(is.finite(kl_res$KL_divergence)), "KL scores should be finite.")

allele_res <- allele_integ(
  fake_result,
  phased = data.frame(
    ID = c("m1", "m2"),
    GT = c("0|1", "0|1"),
    stringsAsFactors = FALSE
  ),
  feature = "gene_symbol"
)

assert_true(inherits(allele_res, "HBomicsAlleleResult"), "allele_integ() should return an HBomicsAlleleResult.")
assert_true(all(c("gene_allele_imb", "HBomics") %in% names(allele_res)), "Allele result should expose gene_allele_imb and HBomics elements.")
assert_true("GT" %in% names(allele_res$HBomics), "Integrated HBomics output should contain the GT column.")
