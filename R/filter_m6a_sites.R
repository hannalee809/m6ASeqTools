#' Filter m6A Sites
#'
#' This function filters m6A sites from m6anet output based on probability_modified criteria.
#'
#' @param data.site_proba_df The datadata frame from m6Anet output containing column 'probability_modified'
#' @param prob_modified The cut-off for probability modified. The default is 0.9 and above.
#' @return filtered dataframe of m6A sites
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @export

filter_m6a_sites <- function(data.site_proba_df, prob_modified=0.9) {
  filt_data <- data.site_proba_df %>%
    dplyr::filter(probability_modified >= prob_modified)

  return(filt_data)
}
