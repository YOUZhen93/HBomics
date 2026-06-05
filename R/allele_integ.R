#' HBomics algorithm: calling clonal allele imbalance events from multiomics data with a hierarchical Bayesian framework
#' Author: Zhen Y
#' @name allele_integ
#' @rdname allele_integ
#' 
#' @title feature-level allele imbalance integration
#'
#' @description Function for allele imbalance ratio integration to feature-level (gene/peak)
#'   with user provided phase information
#'   for feature with multiple SNVs, allele imbalance ratios are weighted by the allele depth from the enrichment data, 
#'   i.e., the most expressed allele weights more on the integrative allele imbalance ratios  
#' @param hbomics The hbomics object from HBomics function,
#' @param phased phased SNV file with two columns: ID contain SNV id in the HBomics object;
#' GT column contain phase information such as "1|0" or "0|1"
#' @param feature indicate the feature column in the HBomics object hb_df data frame;
#' @return An `HBomicsAlleleResult` object with feature-level summaries in
#'   `gene_allele_imb` and the integrated site-level data in `HBomics`.
#' 
#' @export
#' @import dplyr
allele_integ <- function(hbomics, phased, feature, ...) {
  hbomics <- coerce_hbomics_result(hbomics, arg_name = "hbomics")
  validate_phased_input(phased)
  validate_string(feature, "feature")

  hb_df <- hbomics$hb_df
  validate_required_columns(
    hb_df,
    c("mutation_id", "mode", "total_count", "alt_count", feature),
    "hbomics$hb_df"
  )

  if (length(intersect(hb_df$mutation_id, phased$ID)) < 1) {
    abort_hbomics("No mutation found in the phased SNV set. Confirm that IDs match the HBomics result object.")
  }

  # a dict for genotype mapping (to directional numeric values)
  c_dict <- c(-1, 1)
  names(c_dict) <- c("1|0", "0|1")

  # remove SNVs without feature information
  hb_df <- hb_df %>%
    filter(!is.na(.data[[feature]]), .data[[feature]] != "") %>%
    # merge by mutation id
    left_join(phased, by = c("mutation_id" = "ID")) %>%
    # Keep only informative phased genotypes
    filter(GT %in% c("1|0", "0|1")) %>%
    # Apply numeric scoring
    mutate(GT = c_dict[GT])
  
gene_allele_imb <- hb_df %>%
    group_by(.data[[feature]]) %>%
    summarize(
      total_count_agg = sum(total_count, na.rm = TRUE),
      alt_count_agg = sum(alt_count, na.rm = TRUE),
      total_ar = alt_count_agg / total_count_agg,
      
      # Calculate the weighted, haplotype-aware imbalance score
      raw_agg_ar = sum((total_count / total_count_agg) * (mode - 0.5) * GT, na.rm = TRUE),
      
      # Apply directionality imbalance based on the total allele ratio
      allele_imb = if_else(total_ar >= 0.5, abs(raw_agg_ar), -abs(raw_agg_ar)),
      .groups = "drop"
    ) %>%
    # Sort by absolute imbalance (largest shifts at the top)
    arrange(desc(abs(allele_imb))) %>%
    # unify output column names
    rename(
      id = all_of(feature), 
      total_count = total_count_agg, 
      alt_count = alt_count_agg
    )

  new_hbomics_allele_result(
    gene_allele_imb = as.data.frame(gene_allele_imb),
    HBomics = hb_df
  )
}
