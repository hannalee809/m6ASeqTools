#' Compute Wilcoxon Test (Mann-Whitney U test) for Differential Methylation
#'
#' This function computes a statistical test for differential methylation for
#' conditions with replicates. Specifically, it performs a Wilcoxon rank-sum
#' test per transcript to compare weighted modification ratios (WMR) between the
#' two conditions.
#'
#' @param condition1_list A list of data frames for condition 1. Each must contain
#'   `ensembl_transcript_id` and `weighted_mod_ratio`.
#' @param condition2_list A list of data frames for condition 2. Same required
#'   columns as `condition1_list`.
#' @param condition1_name A string naming condition 1.
#' @param condition2_name A string naming condition 2.
#'
#' @return A tibble with one row per transcript containing:
#'   \itemize{
#'     \item `ensembl_transcript_id`
#'     \item `median_diff` — median WMR( condition2 – condition1 )
#'     \item `p_value` — Wilcoxon rank-sum test p-value
#'   }
#'
#' @details Transcripts with all-zero WMR values across replicates are removed.
#'   Transcripts without data in both conditions return `NA` for all outputs.
#'
#' @importFrom dplyr bind_rows mutate group_by filter group_modify
#' @importFrom tibble tibble
#' @importFrom stats wilcox.test
#' @export

run_wmr_differential_test <- function(
  condition1_list,
  condition2_list,
  condition1_name,
  condition2_name
) {

  # first define that combines replicates per condition and labels them
  combine_reps <- function(df_list, condition_name) {
    dplyr::bind_rows(lapply(seq_along(df_list), function(i) {
      df_list[[i]] %>% dplyr::mutate(condition = condition_name, replicate = i)
    }))
  }

  # combine all replicates
  all_wmr <- dplyr::bind_rows(
    combine_reps(condition1_list, condition1_name),
    combine_reps(condition2_list, condition2_name)
  )

  # filter all-zero transcripts
  filtered <- all_wmr %>%
    dplyr::group_by(ensembl_transcript_id) %>%
    dplyr::filter(any(weighted_mod_ratio > 0)) %>%
    ungroup()

  # run Wilcoxon per transcript
  results <- filtered %>%
    dplyr::group_by(ensembl_transcript_id) %>%
    dplyr::group_modify(~ {

      # Skip transcripts that do not have measurements in both conditions
      if (length(unique(.x$condition)) < 2)
        return(tibble::tibble(median_diff = NA, p_value = NA))

      # Perform Wilcoxon rank-sum test comparing WMR between conditions
      # exact = FALSE avoids errors when there are ties or small sample sizes
      wt <- stats::wilcox.test(weighted_mod_ratio ~ condition, data = .x, exact = FALSE)

      # Compute effect size: median difference (condition2 − condition1)
      tibble::tibble(
        median_diff =
          median(.x$weighted_mod_ratio[.x$condition == condition2_name]) -
          median(.x$weighted_mod_ratio[.x$condition == condition1_name]),
        p_value = wt$p.value
      )
    }) %>%
    ungroup()

  return(results)
}

