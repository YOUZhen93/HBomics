#' HBomics algorithm: calling clonal allele imbalance events from multiomics data with a hierarchical Bayesian framework
#' Author: Zhen Y
#' @title HBomics Class
#' @description R6 class for fitting the HBomics model and retrieving validated result objects.
#' @import R6
#' @import rstan
#' @export
HBomics <- R6::R6Class(
  "HBomics",
  public = list(
    # Initialize with input validation
    # requires meta (feature info) and data_list (main input)
    initialize = function(data_list, mutation_id, meta) {
      validated <- validate_hbomics_inputs(data_list, mutation_id, meta)
      private$data_list <- validated$data_list
      private$mutation_id <- validated$mutation_id
      private$meta <- validated$meta
      invisible(self)
    },

    fit_model = function(chains = 4, warmup = 1000, iter = 2000, cores = 1, seed = 123, ...) {
      validate_sampling_args(
        chains = chains,
        warmup = warmup,
        iter = iter
      )

      message(date(), "     INFO: Fitting pre-compiled HBomics model...")
      private$fit <- rstan::sampling(
        stanmodels$HBomics_pal,
        data = private$data_list,
        chains = chains,
        iter = iter,
        warmup = warmup,
        seed = seed,
        cores = cores,
        refresh = 100,
        ...
      )

      invisible(self)
    },
    # checking convergence on posterior samples
    convergence = function() {
      if (is.null(private$fit)) {
        abort_hbomics("Error: HBomics model not fitted. Run $fit_model() first.")
      }

      private$diagnostics <- extract_convergence_diagnostics(
        fit = private$fit,
        mutation_id = private$mutation_id
      )
      private$diagnostics
    },
    # Process the allele imbalance posterior samples
    result = function(KL = FALSE, bound = 0.95) {
      validate_flag(KL, "KL")
      validate_bound(bound)

      if (is.null(private$fit)) {
        abort_hbomics("Error: HBomics model not fitted. Run $fit_model() first.")
      }
      if (is.null(private$diagnostics)) {
        self$convergence()
      }

      summary_data <- extract_hbomics_summary(
        fit = private$fit,
        mutation_id = private$mutation_id,
        meta = private$meta,
        bound = bound
      )

      private$posterior_a <- summary_data$post_a
      private$hb_df <- summary_data$hb_df

      new_hbomics_result(
        hb_df = private$hb_df,
        diagnostics = private$diagnostics,
        post_a = if (KL) private$posterior_a else NULL
      )
    },

    get_fit = function() {
      private$fit
    },

    get_diagnostics = function() {
      private$diagnostics
    },

    get_hb_df = function() {
      private$hb_df
    },

    get_posterior_samples = function() {
      private$posterior_a
    }
  ),
  private = list(
    data_list = NULL,
    mutation_id = NULL,
    meta = NULL,
    fit = NULL,
    posterior_a = NULL,
    hb_df = NULL,
    diagnostics = NULL
  )
)

#' Main wrapper to run the HBomics pipeline
#' @param data_list A named list containing the het. SNV information required by
#'   the Stan model.
#' T: total SNP number (int); 
#' C: total subclone number (int);
#' S: total segment number (int);
#' segment_id: segment id, should be unique (array);
#' subclone_id: subclone id, should be unique (array);  
#' n: SNV number of each segment (array);
#' A: alternative allele count of each SNV (array);
#' N: total allele count of each SNV (array);
#' b: genomic allele ratios of each SNV (array);
#' v_t: subclonal genomic variance (array);
#' v_j: segmental genomic variance (array);
#' lambda: Gamma scale paramter (numeric), default is 0.2;
#' mu0: mapping bias from enrichment data (numeric); 1 - mean(ref allele ratio);
#' sigma0: standard deviation of allele ratios from enrichment data (numeric) 
#' c: smallest median SNV density for segments (int)
#' s: scale parameter in sigmoid function, default is 0.5
#' @param mutation_id Character vector of unique mutation IDs in the same order
#'   as the site-level inputs in `data_list`.
#' @param meta User-defined metadata with a required `mutation_id` column and one
#'   row per mutation.
#' this data frame can contains features (genes/regulatory elements per SNV resides) 
#' and other metrics
#' be sure it includes all mutations from the data_list and mutation_id
#' @param chains Number of Markov chains.
#' @param warmup Number of warmup iterations.
#' @param iter Total number of iterations.
#' @param cores Number of CPU cores to use.
#' @param seed Random seed.
#' @param KL Logical flag indicating whether posterior samples should be retained
#'   for downstream KL-divergence computation.
#' @param bound Credible interval boundary for posterior summaries.
#' @param ... Additional arguments passed to [rstan::sampling()].
#' @return An `HBomicsResult` object containing `hb_df`, `diagnostics`, and
#'   optionally `post_a` when `KL = TRUE`.
#' @export
run_HBomics_pal <- function(data_list, mutation_id, meta, bound = 0.95,
                            chains = 4, warmup = 1000, iter = 2000,
                            cores = 4, seed = 123, KL = FALSE, ...) {
  model_obj <- HBomics$new(
    data_list = data_list,
    mutation_id = mutation_id,
    meta = meta
  )

  model_obj$fit_model(
    chains = chains,
    iter = iter,
    warmup = warmup,
    cores = cores,
    seed = seed,
    ...
  )
  model_obj$convergence()
  model_obj$result(KL = KL, bound = bound)
}
