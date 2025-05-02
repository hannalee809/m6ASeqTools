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
#'   everything case_when
#' @importFrom tidyr replace_na
#' @importFrom magrittr %>%
#' @export

map_relative_tx_regions_to_m6A <- function(data.site_proba_df, df_transcript_lengths) {
  df <- data.site_proba_df %>%
    dplyr::left_join(df_transcript_lengths, by = c("ensembl_transcript_id" = "transcript")) %>%
    dplyr::mutate(
      isin_5utr = dplyr::if_else(!is.na(end_utr5) & transcript_position <= end_utr5, 1L, 0L),
      isin_cds  = dplyr::if_else(!is.na(start_cds) & transcript_position >= start_cds & transcript_position <= end_cds, 1L, 0L),
      isin_3utr = dplyr::if_else(!is.na(start_utr3) & transcript_position >= start_utr3 & transcript_position <= end_utr3, 1L, 0L),

      rel_5utr = dplyr::if_else(isin_5utr == 1 & !is.na(utr5_length) & utr5_length != 0,
                                transcript_position / utr5_length, as.numeric(NA)),
      rel_cds = dplyr::if_else(isin_cds == 1 & !is.na(cds_length_tx) & cds_length_tx != 0,
                               (transcript_position - utr5_length) / cds_length_tx, as.numeric(NA)),
      rel_3utr = dplyr::if_else(isin_3utr == 1 & !is.na(utr3_length) & utr3_length != 0,
                                (transcript_position - (utr5_length + cds_length_tx)) / utr3_length, as.numeric(NA))
    )

  return(df)
}
