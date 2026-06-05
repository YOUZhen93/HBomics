#' HBomics algorithm: calling clonal allele imbalance events from multiomics data with a hierarchical Bayesian framework
#' Author: Zhen Y
#'
#' @description Validation functions for the HBomics package
#' @keywords internal
NULL

required_hbomics_data_list_fields <- function() {
  c(
    "T", "C", "S", "segment_id", "subclone_id", "n",
    "A", "N", "b", "v_t", "v_j", "lambda", "mu0",
    "sigma0", "c", "s"
  )
}

abort_hbomics <- function(..., call. = FALSE) {
  stop(..., call. = call.)
}

validate_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    abort_hbomics(sprintf("Error: '%s' must be TRUE or FALSE.", name))
  }
}

validate_string <- function(x, name) {
  if (!is.character(x) || length(x) != 1 || is.na(x) || !nzchar(x)) {
    abort_hbomics(sprintf("Error: '%s' must be a non-empty string.", name))
  }
}

validate_numeric_scalar <- function(x, name, lower = -Inf, upper = Inf,
                                    inclusive_lower = TRUE, inclusive_upper = TRUE) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || !is.finite(x)) {
    abort_hbomics(sprintf("Error: '%s' must be a finite numeric scalar.", name))
  }

  lower_ok <- if (inclusive_lower) x >= lower else x > lower
  upper_ok <- if (inclusive_upper) x <= upper else x < upper

  if (!lower_ok || !upper_ok) {
    comparator_lower <- if (inclusive_lower) ">=" else ">"
    comparator_upper <- if (inclusive_upper) "<=" else "<"
    abort_hbomics(
      sprintf("Error: '%s' must satisfy %s %s and %s %s.", name, comparator_lower, lower, comparator_upper, upper)
    )
  }
}

validate_integerish_scalar <- function(x, name, lower = -Inf, upper = Inf,
                                       inclusive_lower = TRUE, inclusive_upper = TRUE) {
  validate_numeric_scalar(
    x = x,
    name = name,
    lower = lower,
    upper = upper,
    inclusive_lower = inclusive_lower,
    inclusive_upper = inclusive_upper
  )

  if (!isTRUE(all.equal(x, as.integer(x)))) {
    abort_hbomics(sprintf("Error: '%s' must be an integer-like scalar.", name))
  }
}

validate_required_columns <- function(df, cols, object_name) {
  missing_cols <- setdiff(cols, names(df))
  if (length(missing_cols) > 0) {
    abort_hbomics(
      sprintf(
        "Error: %s is missing required columns: %s.",
        object_name,
        paste(missing_cols, collapse = ", ")
      )
    )
  }
}

validate_bound <- function(bound) {
  validate_numeric_scalar(bound, "bound", lower = 0, upper = 1, inclusive_lower = FALSE, inclusive_upper = FALSE)
}

validate_sampling_args <- function(chains, warmup, iter) {
  validate_integerish_scalar(chains, "chains", lower = 1)
  validate_integerish_scalar(warmup, "warmup", lower = 0)
  validate_integerish_scalar(iter, "iter", lower = 1)

  if (warmup >= iter) {
    abort_hbomics("Error: 'warmup' must be smaller than 'iter'.")
  }
}

validate_hbomics_inputs <- function(data_list, mutation_id, meta) {
  if (!is.list(data_list)) {
    abort_hbomics("'data_list' must be a named list.")
  }

  missing_fields <- setdiff(required_hbomics_data_list_fields(), names(data_list))
  if (length(missing_fields) > 0) {
    abort_hbomics(
      sprintf(
        "Error: 'data_list' is missing required fields: %s.",
        paste(missing_fields, collapse = ", ")
      )
    )
  }

  validate_integerish_scalar(data_list$T, "data_list$T", lower = 1)
  validate_integerish_scalar(data_list$C, "data_list$C", lower = 1)
  validate_integerish_scalar(data_list$S, "data_list$S", lower = 1)
  validate_numeric_scalar(data_list$lambda, "data_list$lambda", lower = 0, inclusive_lower = FALSE)
  validate_numeric_scalar(data_list$mu0, "data_list$mu0", lower = 0, upper = 1)
  validate_numeric_scalar(data_list$sigma0, "data_list$sigma0", lower = 0, inclusive_lower = FALSE)
  validate_integerish_scalar(data_list$c, "data_list$c", lower = 0)
  validate_numeric_scalar(data_list$s, "data_list$s", lower = 0, inclusive_lower = FALSE)

  T <- as.integer(data_list$T)
  C <- as.integer(data_list$C)
  S <- as.integer(data_list$S)

  if (length(mutation_id) != T) {
    abort_hbomics("Error: 'mutation_id' must have the same length as 'data_list$T'.")
  }
  if (anyNA(mutation_id) || anyDuplicated(mutation_id)) {
    abort_hbomics("Error: 'mutation_id' must contain unique, non-missing values.")
  }

  vector_requirements <- list(
    "segment_id" = T,
    "A" = T,
    "N" = T,
    "b" = T,
    "subclone_id" = S,
    "n" = S,
    "v_t" = C,
    "v_j" = S
  )

  for (field_name in names(vector_requirements)) {
    if (length(data_list[[field_name]]) != vector_requirements[[field_name]]) {
      abort_hbomics(
        sprintf(
          "Error: 'data_list$%s' must have length %d.",
          field_name,
          vector_requirements[[field_name]]
        )
      )
    }
  }

  if (anyNA(data_list$segment_id) || any(data_list$segment_id < 1 | data_list$segment_id > S, na.rm = TRUE)) {
    abort_hbomics("Error: 'data_list$segment_id' values must be between 1 and 'data_list$S'.")
  }
  if (anyNA(data_list$subclone_id) || any(data_list$subclone_id < 1 | data_list$subclone_id > C, na.rm = TRUE)) {
    abort_hbomics("Error: 'data_list$subclone_id' values must be between 1 and 'data_list$C'.")
  }
  if (any(data_list$A < 0 | data_list$N < 0 | data_list$A > data_list$N, na.rm = TRUE)) {
    abort_hbomics("Error: 'data_list$A' and 'data_list$N' must be non-negative and satisfy A (alt count) <= N (total count).")
  }
  if (any(data_list$b < 0 | data_list$b > 1, na.rm = TRUE)) {
    abort_hbomics("Error: 'data_list$b' values must be between 0 and 1.")
  }
  if (any(data_list$v_t < 0, na.rm = TRUE) || any(data_list$v_j < 0, na.rm = TRUE)) {
    abort_hbomics("Error: 'data_list$v_t' and 'data_list$v_j' must be non-negative.")
  }

  if (!is.data.frame(meta)) {
    abort_hbomics("Error: 'meta' must be a data.frame.")
  }
  validate_required_columns(meta, "mutation_id", "meta")
  if (anyDuplicated(meta$mutation_id) || anyNA(meta$mutation_id)) {
    abort_hbomics("Error: 'meta$mutation_id' must contain unique, non-missing values.")
  }
  if (!setequal(mutation_id, meta$mutation_id)) {
    abort_hbomics("Error: 'meta$mutation_id' must match 'mutation_id' exactly.")
  }

  meta <- meta[match(mutation_id, meta$mutation_id), , drop = FALSE]

  list(
    data_list = data_list,
    mutation_id = as.character(mutation_id),
    meta = meta
  )
}

coerce_hbomics_result <- function(x, require_posterior = FALSE, arg_name = "hbomics") {
  if (!is.list(x) || is.null(x$hb_df) || !is.data.frame(x$hb_df)) {
    abort_hbomics(
      sprintf("Error: '%s' must be an HBomicsResult object containing a data-frame 'hb_df'.", arg_name)
    )
  }

  if (require_posterior && is.null(x$post_a)) {
    abort_hbomics(
      sprintf("Error: '%s' must contain posterior samples in 'post_a'. Re-run HBomics with flag: KL = TRUE.", arg_name)
    )
  }

  x
}

validate_phased_input <- function(phased) {
  if (!is.data.frame(phased)) {
    abort_hbomics("Error: 'phased' must be a data.frame.")
  }
  validate_required_columns(phased, c("ID", "GT"), "phased")
}
