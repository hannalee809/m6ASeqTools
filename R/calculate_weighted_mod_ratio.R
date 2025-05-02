#' Calculate Weighted Modification Ratio
#'
#' Calculates the sum of a given modification ratio column for each transcript,
#' then normalizes it by transcript length.
#'
#' @param data_frame A data frame with m6A site-level data with transcript regions (from `get_transcript_region_lengths()`)
#' @param mod_ratio_column Character string, name of the column with modification ratios (default: "mod_ratio")
#'
#' @return A data frame with weighted modification ratios per transcript
#'
#' @importFrom dplyr group_by summarise
#' @importFrom magrittr %>%
#' @export
calculate_weighted_mod_ratio <- function(data_frame, mod_ratio_column = "mod_ratio") {
  # Group by transcript-related info
  grouped <- data_frame %>%
    dplyr::group_by(
      ensembl_transcript_name,
      ensembl_transcript_id,
      ensembl_gene_id,
      ensembl_gene_name,
      gene_biotype,
      transcript_length
    )

  # Sum the modification ratio column dynamically
  weighted_df <- grouped %>%
    dplyr::summarise(
      sum_mod = sum(.data[[mod_ratio_column]], na.rm = TRUE),
      n_sites = dplyr::n(),
      sites_in_5utr = sum(isin_5utr),
      sites_in_cds = sum(isin_cds),
      sites_in_3utr = sum(isin_3utr),
      .groups = "drop"
    )

  # Calculate weighted modification ratio
  weighted_df$weighted_mod_ratio <- weighted_df$sum_mod / weighted_df$transcript_length

  return(weighted_df)
}
