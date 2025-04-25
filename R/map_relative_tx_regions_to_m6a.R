#' Map Relative Transcript Regions to m6A Sites
#'
#' This function calculates the relative transcript regions of the 5' UTR, CDS,
#' and 3' UTR for each transcript using the region length df from
#' `get_transcript_region_lengths()`
#'
#' @param data.site_proba_df The dataframe from m6Anet output containing columns
#'   `ensembl_transcript_id` and `transcript_position`
#' @param df_transcript_lengths The output dataframe from
#'   `get_transcript_region_lengths()` with transcript IDs and lengths/positions
#'   of UTRs and CDS regions.
#' @return A data frame with annotated transcript positions per m6A site,
#'   relative to the region of the transcript
#'
#' @importFrom GenomicFeatures makeTxDbFromGFF cdsBy fiveUTRsByTranscript
#'   threeUTRsByTranscript
#' @importFrom dplyr group_by summarise rename full_join mutate across
#'   everything
#' @importFrom tidyr replace_na
#' @importFrom magrittr %>%
#' @export
#' @examples
#' mapped_regions <- map_relative_tx_regions_to_m6A(data.site_proba_df, tx_regions)

map_relative_tx_regions_to_m6A <- function(data.site_proba_df, df_transcript_lengths) {
  df <- df %>%
    dplyr::inner_join(df_transcript_lengths, by = c("ensembl_transcript_id"="transcript")) %>%
    dplyr::mutate(
      isin_5utr = transcript_position <= end_utr5,
      isin_cds = (transcript_position >= start_cds) & (transcript_position <= end_cds),
      isin_3utr = (transcript_position >= start_utr3) & (transcript_position <= end_utr3),


      rel_5utr = ifelse(isin_5utr, transcript_position / utr5_length, 0),
      rel_cds = ifelse(isin_cds, (transcript_position - utr5_length) / cds_length_tx, 0),
      rel_3utr = ifelse(isin_3utr, (transcript_position - (utr5_length + cds_length_tx)) / utr3_length, 0),

    ) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ tidyr::replace_na(., 0)))
  return(df)
}
