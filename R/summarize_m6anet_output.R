#' Summarize m6anet output
#'
#' This function summarizes the output of m6anet. This results in a html file with plots and
#' metrics. You can also put pooled data from m6anet as input.
#'
#' @param data.site_proba_df The 'datadata frame from m6Anet output.
#' @param output_file Name of the output HTML file.
#' @return Renders an HTML summary QC report.
#' @export

summarize_m6anet_output <- function(data.site_proba_df, output_file="m6anet_qc_report.html", output_dir = getwd()) {
  rmarkdown::render(
    input = system.file("rmd", "m6anet_qc_template.Rmd", package = "m6AnetAnalyzer"),
    output_file = output_file,
    output_dir = output_dir,
    params = list(data.site_proba_df = data.site_proba_df),
    envir = new.env(parent = globalenv())
  )
}
