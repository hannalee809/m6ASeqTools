#' Get Transcript Region Lengths
#'
#' This function calculates the lengths of the 5' UTR, CDS, and 3' UTR for each transcript
#' from a GTF annotation file using GenomicFeatures and summarizes them in a data frame.
#'
#' @param gtf_path The path to the GTF file of Gencode annotations.
#' @return A data frame with transcript IDs and lengths/positions of UTRs and CDS regions.
#'
#' @importFrom GenomicFeatures makeTxDbFromGFF cdsBy fiveUTRsByTranscript threeUTRsByTranscript
#' @importFrom dplyr group_by summarise rename full_join mutate across
#' @importFrom tidyr replace_na
#' @importFrom stats setNames
#' @importFrom methods as
#' @importFrom S4Vectors DataFrame
#' @importFrom magrittr %>%
#' @import GenomicRanges
#' @import IRanges
#' @export
#' @examples
#' gtf_path <- 'gencode.v45.chr_patch_hapl_scaff.annotation.gtf'
#' tx_region <- get_transcript_region_lengths(gtf_path)

get_transcript_region_lengths <- function(gtf_path) {
  # Build TxDb from GTF
  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_path, format = "gtf")

  # Get CDS
  cds_gr <- GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)
  cds_df <- as.data.frame(cds_gr) %>%
    dplyr::group_by(group_name) %>%
    dplyr::summarise(cds_length_tx = sum(width), .groups = "drop") %>%
    dplyr::rename(transcript = group_name)

  # Get 5' UTR
  utr5_gr <- GenomicFeatures::fiveUTRsByTranscript(txdb, use.names = TRUE)
  utr5_df <- as.data.frame(utr5_gr) %>%
    dplyr::group_by(group_name) %>%
    dplyr::summarise(utr5_length = sum(width), .groups = "drop") %>%
    dplyr::rename(transcript = group_name)

  # Get 3' UTR
  utr3_gr <- GenomicFeatures::threeUTRsByTranscript(txdb, use.names = TRUE)
  utr3_df <- as.data.frame(utr3_gr) %>%
    dplyr::group_by(group_name) %>%
    dplyr::summarise(utr3_length = sum(width), .groups = "drop") %>%
    dplyr::rename(transcript = group_name)

  # Merge all and fill NAs with 0
  tx_regions <- dplyr::full_join(utr5_df, utr3_df, by = "transcript") %>%
    dplyr::full_join(cds_df, by = "transcript") %>%
    dplyr::mutate(across(c(utr5_length, utr3_length, cds_length_tx), ~tidyr::replace_na(., 0)))

  # Add region positions (start/end)
  tx_regions <- tx_regions %>%
    dplyr::mutate(
      start_utr5 = ifelse(is.na(utr5_length), NA, 1),
      end_utr5 = ifelse(is.na(utr5_length), NA, utr5_length),
      start_cds = ifelse(is.na(cds_length_tx), NA, ifelse(is.na(utr5_length), 1, utr5_length + 1)),
      end_cds = ifelse(is.na(cds_length_tx), NA, ifelse(is.na(utr5_length), cds_length_tx, utr5_length + cds_length_tx - 1)),
      start_utr3 = ifelse(is.na(utr3_length), NA, ifelse(is.na(end_cds), 1, end_cds + 1)),
      end_utr3 = ifelse(is.na(utr3_length), NA, ifelse(is.na(end_cds), utr3_length, end_cds + utr3_length))
    )

  return(tx_regions)
}

