#' Create BED and BedGraph Files from m6Anet Output
#'
#' This function takes transcript-level m6A positions, maps them to genomic
#' coordinates using a GTF file, creates standard BED6 and BedGraph data frames,
#' and optionally writes them to files.
#'
#' @param data.site_proba_df A data frame containing m6Anet output, must include
#'   'ensembl_transcript_id', 'transcript_position', 'ensembl_gene_name',
#'   'mod_ratio'.
#' @param gtf_path Path to the GENCODE/Ensembl GTF or GFF3 annotation file.
#' @param output_bed Optional file path to write the BED6 file (for IGV feature
#'   track).
#' @param output_bedgraph Optional file path to write the BedGraph file (for
#'   coverage/score track).
#' @return A named list containing two data frames: -
#'   m6a_df_with_genomic_coordinates: The input data merged with genomic
#'   coordinates (chr, genomic_pos, strand). The `name` column is constructed by
#'   joining the transcript ID and the transcript position with an underscore
#'   (i.e., `transcript_id_transcript_position`). - bed_df: The final, formatted
#'   BED6 data frame.
#' @importFrom dplyr select mutate transmute filter left_join
#' @importFrom GenomicFeatures exonsBy mapFromTranscripts
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom AnnotationDbi select
#' @importFrom readr write_delim
#' @importFrom txdbmaker makeTxDbFromGFF
#' @export

create_bed_file <- function(data.site_proba_df,
                            gtf_path,
                            output_bed = NULL,
                            output_bedgraph = NULL) {

  # Create the Transcript Database (TxDb) from gtf_path
  txdb <- txdbmaker::makeTxDbFromGFF(gtf_path)

  # Create row id in m6a dataframe for ID matching
  data.site_proba_df$row_id <- 1:nrow(data.site_proba_df)

  # Convert m6a df to GRanges object
  m6a_gr <- GenomicRanges::GRanges(
    seqnames = data.site_proba_df$ensembl_transcript_id,
    ranges = IRanges::IRanges(
      start = as.numeric(data.site_proba_df$transcript_position),
      end = as.numeric(data.site_proba_df$transcript_position)
    )
  )

  names(m6a_gr) <- data.site_proba_df$ensembl_transcript_id

  # Extract exon-based transcript structure from the TxDb
  transcripts_granges_list <- GenomicFeatures::exonsBy(txdb, by = "tx")


  # Get ALL transcript-level annotations from the TxDb
  all_transcripts_df <- AnnotationDbi::select(
    txdb,
    keys = names(transcripts_granges_list), # Use the internal IDs as keys
    columns = "TXNAME",                    # Request the actual transcript names
    keytype = "TXID"
  )

  # Create a mapping vector: [Internal ID] -> [Actual TXNAME]
  id_map <- all_transcripts_df$TXNAME
  names(id_map) <- all_transcripts_df$TXID

  # Assign the actual transcript names to the GRangesList
  names(transcripts_granges_list) <- id_map[names(transcripts_granges_list)]

  # Run the mapping function
  m6a_gr_genomic_coords <- GenomicFeatures::mapFromTranscripts(m6a_gr, transcripts_granges_list)

  # Convert to Data Frame
  m6a_gr_genomic_coords <- as.data.frame(m6a_gr_genomic_coords, row.names = NULL)


  # Merge the original m6A data to keep all columns
  # This requires a unique ID for each m6A site, which the 'queryHits' column provides.
  m6a_gr_genomic_coords$query_hit_id <- 1:nrow(m6a_gr_genomic_coords)


  # Prepare the Merged m6A Data with Genomic Coordinates
  # The 'xHits' column in final_genomic_df corresponds to the row index
  # in the original m6a_gr_genomic_coords (which is the same order as data.site_proba_df).

  # We will use the 'xHits' column to merge the data.
  m6a_gr_genomic_coords <- m6a_gr_genomic_coords %>%
    # Select and rename columns for clarity
    dplyr::select(
      chr = seqnames,
      genomic_pos = start, # Position is start
      strand,
      original_m6a_row_index = xHits # The index of the row in the original df
    ) %>%
    # Join back to the original m6A data frame to get the transcript ID
    dplyr::mutate(
      ensembl_transcript_id = data.site_proba_df$ensembl_transcript_id[original_m6a_row_index]
    )


  # merge with original m6A df
  merged_genomic_coords <- dplyr::left_join(
    data.site_proba_df,
    m6a_gr_genomic_coords %>% select(-ensembl_transcript_id), # remove duplicate
    by = c("row_id" = "original_m6a_row_index")
  )

  # create name column in m6A df with the transcript name and position
  # this will correspond to the name column in the bed file
  merged_genomic_coords <- merged_genomic_coords %>%
    dplyr::mutate(
      name = paste0(ensembl_transcript_name, "_", transcript_position)
    )

  # BED6 format (chr, start, end, name, score, strand)
  bed_df <- merged_genomic_coords %>%
    dplyr::mutate(
      start = genomic_pos - 1, # 0-based
      end = genomic_pos, # 1-based end (for single base),
      name = paste0(ensembl_transcript_name, "_", transcript_position),
      score = mod_ratio * 1000 # Should be between 0 and 1000 for IGV color
    ) %>%
    dplyr::select(
      chr,
      start,
      end,
      name,
      score,
      strand
    )

  bedgraph_df <- merged_genomic_coords %>%
    dplyr::transmute(
      chr,
      start = genomic_pos - 1, # 0-based
      end = genomic_pos,
      score = mod_ratio # make the mod ratio the score
    )

  # write files
  if (!is.null(output_bed)) {
    write.table(bed_df,
                file = output_bed,
                sep = "\t",
                quote = FALSE,
                row.names = FALSE,
                col.names = FALSE)
  }

  if (!is.null(output_bedgraph)) {
    write.table(bedgraph_df,
                file = output_bedgraph,
                sep = "\t",
                quote = FALSE,
                row.names = FALSE,
                col.names = FALSE)
  }


  # return merged m6a data with genomic coordinates and
  # bed file into a named list

  results_list <- list(
    m6a_df_with_genomic_coordinates = merged_genomic_coords,
    bed_df = bed_df
  )

  return(results_list)
}
