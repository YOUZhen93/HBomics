library(HBomics)

assert_error_contains <- function(expr, pattern) {
  msg <- tryCatch(
    {
      force(expr)
      NULL
    },
    error = function(e) conditionMessage(e)
  )

  if (is.null(msg)) {
    stop("Expected an error but the expression succeeded.", call. = FALSE)
  }
  if (!grepl(pattern, msg, fixed = TRUE)) {
    stop(sprintf("Expected error containing '%s' but got '%s'.", pattern, msg), call. = FALSE)
  }
}

valid_data_list <- list(
  T = 2,
  C = 1,
  S = 1,
  segment_id = c(1, 1),
  subclone_id = 1,
  n = 2,
  A = c(1, 1),
  N = c(2, 2),
  b = c(0.50, 0.60),
  v_t = 0.10,
  v_j = 0.20,
  lambda = 0.20,
  mu0 = 0.50,
  sigma0 = 0.10,
  c = 1,
  s = 0.50
)

valid_meta <- data.frame(
  mutation_id = c("m1", "m2"),
  total_count = c(10, 12),
  alt_count = c(6, 7),
  gene_symbol = c("G1", "G2"),
  stringsAsFactors = FALSE
)

missing_b <- valid_data_list
missing_b$b <- NULL

assert_error_contains(
  run_HBomics_pal(
    data_list = missing_b,
    mutation_id = c("m1", "m2"),
    meta = valid_meta
  ),
  "'data_list' is missing required fields"
)

assert_error_contains(
  run_HBomics_pal(
    data_list = valid_data_list,
    mutation_id = "m1",
    meta = valid_meta
  ),
  "'mutation_id' must have the same length"
)

assert_error_contains(
  run_HBomics_pal(
    data_list = valid_data_list,
    mutation_id = c("m1", "m2"),
    meta = valid_meta,
    bound = 1
  ),
  "'bound' must satisfy"
)
