#' Summarize Weighted Modification Ratio and DGE and Provide annotations for clusters
#'
#' This function filters m6A sites from m6anet output based on probability_modified criteria.
#'
#' @param df Output from weighted_mod_ratio_and_DGE
#' @param group_1_name Name for group 1
#'  @param group_2_name Name for group 2
#' @return filtered dataframe of m6A sites
#' @importFrom dplyr group_by summarize mutate n bind_rows
#' @importFrom magrittr %>%
#' @importFrom rlang sym !!
#' @export

summarize_weighted_mod_ratio_and_DGE <- function(df, group1_name = "Group1", group2_name = "Group2", combine = TRUE) {
  group_list <- list()
  all_biotype <- list()
  all_region <- list()

  # Define filtering rules for each group
  filters <- list(
    upreg_hyperm6A_group1        = df[, 22] == TRUE  & df[, 20] > 0 & df[, 19] > 0,
    upreg_group1_hyperm6a_group2 = df[, 22] == TRUE  & df[, 20] > 0 & df[, 19] < 0,
    upreg_hyperm6A_group2        = df[, 22] == TRUE  & df[, 20] < 0 & df[, 19] < 0,
    upreg_group2_hyperm6a_group1 = df[, 22] == TRUE  & df[, 20] < 0 & df[, 19] > 0,
    noDGE_hyperm6A_group1        = df[, 22] == FALSE & df[, 20] > 0,
    noDGE_hyperm6A_group2        = df[, 22] == FALSE & df[, 20] < 0
  )

  for (name in names(filters)) {
    group_df <- df[filters[[name]], ]

    # Biotype summary
    biotype_summary <- group_df %>%
      dplyr::group_by(gene_biotype) %>%
      dplyr::summarize(count = dplyr::n(), .groups = "drop") %>%
      dplyr::mutate(frequency = count / sum(count),
                    group = name)

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

  # If combine = TRUE, return merged summaries too
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

