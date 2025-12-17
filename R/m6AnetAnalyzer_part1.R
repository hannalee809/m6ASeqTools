#' Automate core steps of m6ASeqTools pipeline (Part 1)
#'
#' This function performs a series of common m6ASeqTools analyses including QC summary, filtering, region annotation,
#' chromosome location annotation, and weighted modification ratio calculation.
#'
#' @param m6Aout Data frame from m6Anet output (e.g., `data.site_proba.csv`)
#' @param gtf_path File path to GTF annotation file (e.g., GENCODE)
#' @param output_directory Path where outputs will be saved
#' @param name A string identifier used to name output files and the resulting object
#' @param probability Numeric threshold for high-confidence modification sites (default = 0.9)
#'
#' @return Writes multiple output files (CSV, PDF, HTML) to `output_directory` and assigns
#'         the weighted modification ratio dataframe to the global environment
#' @export
m6AnetAnalyzer_part1 <- function(m6Aout,
                    gtf_path,
                    output_directory,
                    name,
                    probability = 0.9) {

  # Create output directory if it doesn't exist
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
  }

  # Set working directory temporarily
  old_dir <- getwd()
  setwd(output_directory)
  on.exit(setwd(old_dir), add = TRUE)

  # 1. Summarize m6Anet output (generates HTML QC report)
  print("Summarizing m6Anet output...")
  summarize_m6anet_output(
    data.site_proba_df = m6Aout,
    output_file = paste0("QC_", name, ".html"),
    output_dir = output_directory
  )
  print("Done")

  # 2. Filter by modification probability
  print("Filtering sites based on probability cutoff...")
  filtered_df <- filter_m6a_sites(m6Aout, prob_modified = probability)
  write.csv(filtered_df, paste0("filtered_m6a_", name, ".csv"), row.names = FALSE)
  print("Done")

  # 3. Get transcript region lengths from GTF
  print("Extracting transcript regions from annotation...")
  tx_df <- get_transcript_region_lengths(gtf_path)
  write.csv(tx_df, paste0("transcript_region_lengths_", name, ".csv"), row.names = FALSE)
  print("Done")

  # 4. Map filtered sites to transcript regions
  print("Mapping m6A sites along the region...")
  map_positions <- map_relative_tx_regions_to_m6A(filtered_df, tx_df)
  write.csv(map_positions, paste0("m6a_relative_tx_regions_", name, ".csv"), row.names = FALSE)
  print("Done")

  # 5. Plot relative transcript position density
  print("Generating region plot...")
  plot_density <- plot_relative_positions(map_positions, "Transcriptomic m6A Distribution Across Regions")
  ggplot2::ggsave(paste0("Transcript_Region_Density_", name, ".pdf"), plot = plot_density)
  print("Done")

  # 6. Chromosome location annotation
  print("Calculating chromosom location...")
  calculate_chromosome_location(
    filtered_df,
    output_csv = paste0("chromosome_location_", name, ".csv"),
    output_plot = paste0("chromosome_location_", name, ".pdf")
  )
  print("Done")

  # 7. Weighted modification ratio per transcript
  print("Calculating weighted mod ratio...")
  wmr <- calculate_weighted_mod_ratio(map_positions)
  write.csv(wmr, paste0("weighted_mod_ratio_", name, ".csv"), row.names = FALSE)

  # Save result to global environment for convenience
  assign(paste0("weighted_mod_ratio_", name), wmr, envir = .GlobalEnv)
  print("Pipeline finished.")
}
