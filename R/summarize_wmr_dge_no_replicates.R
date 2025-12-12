#' Summarize Weighted Modification Ratio and DGE and Provide annotations for clusters
#'
#' This function groups transcripts into functional categories based on gene expression and m6A methylation changes,
#' and summarizes transcript biotype and regional site distributions.
#'
#' @param df Output from weighted_mod_ratio_and_DGE
#' @param group1_name Name for group 1
#' @param group2_name Name for group 2
#' @param log2FC_dge Column name for log2 fold change (DGE)
#' @param log2fc_wmr Column name for log2 fold change (weighted mod ratio)
#' @param sig_col Column name for significance boolean (e.g., adjusted p < 0.05)
#' @param combine Logical; whether to return combined summaries
#'
#' @return A list containing filtered data frames and summaries
#' @export

summarize_wmr_dge_no_replicates <- function(df,
                                                 group1_name = "Group1",
                                                 group2_name = "Group2",
                                                 log2FC_dge,
                                                 log2fc_wmr,
                                                 sig_col = "significant",
                                                 combine = TRUE) {
  group_list <- list()
  all_biotype <- list()
  all_region <- list()

  # Use tidy evaluation to extract columns dynamically
  log2FC_dge <- rlang::sym(log2FC_dge)
  log2fc_wmr <- rlang::sym(log2fc_wmr)
  sig_col <- rlang::sym(sig_col)

  # Define filtering rules as expressions
  filters <- list(
    upreg_hyperm6A_group1        = rlang::expr((!!sig_col) == TRUE  & (!!log2FC_dge) > 0 & (!!log2fc_wmr) > 0),
    upreg_group1_hyperm6a_group2 = rlang::expr((!!sig_col) == TRUE  & (!!log2FC_dge) > 0 & (!!log2fc_wmr) < 0),
    upreg_hyperm6A_group2        = rlang::expr((!!sig_col) == TRUE  & (!!log2FC_dge) < 0 & (!!log2fc_wmr) < 0),
    upreg_group2_hyperm6a_group1 = rlang::expr((!!sig_col) == TRUE  & (!!log2FC_dge) < 0 & (!!log2fc_wmr) > 0),
    noDGE_hyperm6A_group1        = rlang::expr((!!sig_col) == FALSE & (!!log2fc_wmr) > 0),
    noDGE_hyperm6A_group2        = rlang::expr((!!sig_col) == FALSE & (!!log2fc_wmr) < 0)
  )

  for (name in names(filters)) {
    group_df <- dplyr::filter(df, !!filters[[name]])

    # Biotype summary
    biotype_summary <- group_df %>%
      dplyr::group_by(gene_biotype) %>%
      dplyr::summarize(count = dplyr::n(), .groups = "drop") %>%
      dplyr::mutate(frequency = count / sum(count), group = name)

    # Region summary
    region_summary <- tibble::tibble(
      Region = c("5'UTR", "CDS", "3'UTR", "Total Sites"),
      !!group1_name := c(
        sum(group_df[[9]], na.rm = TRUE),
        sum(group_df[[10]], na.rm = TRUE),
        sum(group_df[[11]], na.rm = TRUE),
        sum(group_df[[8]], na.rm = TRUE)
      ),
      !!group2_name := c(
        sum(group_df[[15]], na.rm = TRUE),
        sum(group_df[[16]], na.rm = TRUE),
        sum(group_df[[17]], na.rm = TRUE),
        sum(group_df[[14]], na.rm = TRUE)
      ),
      group = name
    )

    group_list[[name]] <- list(
      data = group_df,
      biotype_summary = biotype_summary,
      region_summary = region_summary
    )

    if (combine) {
      all_biotype[[name]] <- biotype_summary
      all_region[[name]] <- region_summary
    }
  }

  if (combine) {
    return(list(
      groups = group_list,
      biotype_summary_all = dplyr::bind_rows(all_biotype),
      region_summary_all = dplyr::bind_rows(all_region)
    ))
  } else {
    return(group_list)
  }
}
