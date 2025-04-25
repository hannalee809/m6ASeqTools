#' Compare weighted mod ratio for two groups
#'
#' Merges two data frames on commonly modified transcripts and calculates the
#' log2 fold change of the weighted modification ratio between groups.
#'
#' @param group1_df Data frame for Group 1 (experimental), must include
#'   `weighted_mod_ratio`
#' @param group2_df Data frame for Group 2 (control), must include
#'   `weighted_mod_ratio`
#' @param group1_name Name of Group 1 (used as suffix in output columns)
#' @param group2_name Name of Group 2 (used as suffix in output columns)
#'
#' @return A data frame with merged columns and `log2fc_weighted_mod_ratio`
#' @export
#' @examples
#' A_vs_B_log2fc_weighted_mod_ratio <- log2fc_weighted_mod_ratio(A_weighted_mod_df, B_weighted_mod_df, group1_name = "A", group2_name = "B")
log2fc_weighted_mod_ratio <- function(group1_df, group2_df,
                                      group1_name = "Group1",
                                      group2_name = "Group2") {
  # Merge on common identifiers
  common_df <- merge(
    group1_df, group2_df,
    by = c("ensembl_transcript_name", "ensembl_transcript_id", "ensembl_gene_id",
           "ensembl_gene_name", "gene_biotype", "transcript_length"),
    suffixes = c(paste0(".", group1_name), paste0(".", group2_name))
  )

  # Dynamically compute log2FC using the column names
  col1 <- paste0("weighted_mod_ratio.", group1_name)
  col2 <- paste0("weighted_mod_ratio.", group2_name)

  # Avoid divide-by-zero or NaN by adding a small pseudocount if needed
  common_df$log2fc_weighted_mod_ratio <- log2((common_df[[col1]] + 1e-6) / (common_df[[col2]] + 1e-6))

  return(common_df)
}
