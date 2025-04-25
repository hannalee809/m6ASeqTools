#' Summarize m6anet output
#'
#' This function summarizes the output of m6anet. This results in a html file with plots and
#' metrics. You can also put pooled data from m6anet as input.
#'
#' @param data.site_proba_df The 'datadata frame from m6Anet output.
#' @param output_file Name of the output HTML file.
#' @return Renders an HTML summary QC report.
#' @examples
#' summarize_m6anet_output(data.site_proba_df, "qc_report.html")


summarize_m6anet_output <- function(data.site_proba_df, output_file="m6anet_qc_report.html") {
  rmarkdown::render(
    input = system.file("rmd", "m6anet_qc_template.Rmd", package = "m6ASeqTools"),
    output_file = output_file,
    params = list(data.site_proba_df = data.site_proba_df),
    envir = new.env(parent = globalenv())
  )
}

