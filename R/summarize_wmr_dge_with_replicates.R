#' Summarize Weighted Modification Ratio and DGE and Provide annotations for clusters
#'
#' This function groups transcripts into functional categories based on gene expression and m6A methylation changes,
#' and summarizes transcript biotype and regional site distributions.
#'
#' @param df Output from weighted_mod_ratio_and_DGE
#' @param group1_name Name for group 1
#' @param group2_name Name for group 2
#' @param log2FC_dge Column name for log2 fold change (DGE)
#' @param diff_wmr Column name for the calculated difference in median weighted mod ratios
#' @param sig_col Column name for significance boolean (e.g., adjusted p < 0.05)
#'
#' @importFrom dplyr filter group_by summarize mutate n
#' @importFrom rlang expr sym
#'
#' @return A list containing filtered data frames and summaries
#' @export

summarize_wmr_dge_with_replicates <- function(df,
                                              group1_name = "Group1",
                                              group2_name = "Group2",
                                              log2FC_dge,
                                              log2fc_wmr = "median_diff",
                                              sig_col = "significant") {

  # Use tidy evaluation to extract columns dynamically
  log2FC_dge <- rlang::sym(log2FC_dge)
  log2fc_wmr <- rlang::sym(log2fc_wmr)
  sig_col <- rlang::sym(sig_col)

  group_list <- list()

  # Define filtering rules as expressions
  filters <- list(
    sig_upreg_sig_hyperm6A_group1        = rlang::expr((!!sig_col) == TRUE  & (!!log2FC_dge) > 0 & (!!log2fc_wmr) > 0),
    sig_upreg_group1_sig_hyperm6a_group2 = rlang::expr((!!sig_col) == TRUE  & (!!log2FC_dge) > 0 & (!!log2fc_wmr) < 0),
    sig_upreg_sig_hyperm6A_group2        = rlang::expr((!!sig_col) == TRUE  & (!!log2FC_dge) < 0 & (!!log2fc_wmr) < 0),
    sig_upreg_group2_sig_hyperm6a_group1 = rlang::expr((!!sig_col) == TRUE  & (!!log2FC_dge) < 0 & (!!log2fc_wmr) > 0),
    noDGE_not_sig_hyperm6A_group1        = rlang::expr((!!sig_col) == FALSE & (!!log2fc_wmr) > 0),
    noDGE_not_sig_hyperm6A_group2        = rlang::expr((!!sig_col) == FALSE & (!!log2fc_wmr) < 0)
  )

  for (name in names(filters)) {
    group_df <- dplyr::filter(df, !!filters[[name]])

    # Biotype summary
    biotype_summary <- group_df %>%
      dplyr::group_by(gene_biotype) %>%
      dplyr::summarize(count = dplyr::n(), .groups = "drop") %>%
      dplyr::mutate(frequency = count / sum(count), group = name)

    group_list[[name]] <- list(
      data = group_df,
      biotype_summary = biotype_summary
    )
  }

  return(group_list)
}
