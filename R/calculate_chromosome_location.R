#' Calculate Chromosome Location
#'
#' Adds chromosome annotation to a data frame of high-confidence m6A sites using `org.Hs.eg.db`.
#'
#' @param data.site_proba_df m6anet output dataframe containing at least a column 'ensembl_gene_name'.
#' @param output_csv Optional path to save chromosome summary CSV.
#' @param output_plot Optional path to save barplot of chromosome distribution.
#'
#' @return chr_locations Merged data frame with chromosome info.
#'
#' @importFrom dplyr select distinct left_join filter group_by summarize mutate
#' @importFrom ggplot2 ggplot aes geom_bar labs theme_minimal ggsave
#' @importFrom magrittr %>%
#' @import org.Hs.eg.db
#' @examples
#' # Load example data
#' df_m6a <- read.csv(system.file("extdata", "sampleA_site_proba.csv", package = "m6ASeqTools"))
#'
#' # Run chromosome location summary
#' df_with_chr <- calculate_chromosome_location(
#'   df_m6a,
#'   output_csv = tempfile(fileext = ".csv"),
#'   output_plot = tempfile(fileext = ".pdf")
#' )
#' head(df_with_chr)
#' @export

calculate_chromosome_location <- function(data.site_proba_df,
                                          output_csv = NULL,
                                          output_plot = NULL) {
  # Extract uniquely modified genes
  gene_df <- dplyr::select(data.site_proba_df, ensembl_gene_name) %>%
    dplyr::distinct()

  # Retrieve chromosome info
  gene_locations <- AnnotationDbi::select(org.Hs.eg.db,
                                          keys = gene_df$ensembl_gene_name,
                                          columns = c("SYMBOL", "CHR"),
                                          keytype = "SYMBOL")
  # Merge chromosome info with m6a data
  chr_locations <- dplyr::left_join(data.site_proba_df, gene_locations,
                                    by = c("ensembl_gene_name" = "SYMBOL"))

  # Merge chromosome info with uniquely modified gene data
  gene_df <- dplyr::left_join(gene_df, gene_locations,
                              by = c("ensembl_gene_name" = "SYMBOL"))

  # Summarize by chromosome
  chr_summary <- gene_df %>%
    dplyr::filter(!is.na(CHR)) %>%
    dplyr::group_by(CHR) %>%
    dplyr::summarize(n = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(percentage = (n / sum(n)) * 100)

  chr_summary$CHR <- factor(chr_summary$CHR,
                            levels = c(as.character(1:22), "X", "Y"))

  p <- ggplot2::ggplot(chr_summary, ggplot2::aes(x = CHR, y = percentage)) +
    ggplot2::geom_bar(stat = "identity", fill = "cornflowerblue") +
    ggplot2::labs(title = "Chromosome Location Distribution (Percentage)",
                  x = "Chromosome", y = "Percentage of Genes (%)") +
    ggplot2::theme_light()

  if (!is.null(output_plot)) {
    ggplot2::ggsave(output_plot, plot = p, width = 8, height = 6)
  }

  if (!is.null(output_csv)) {
    utils::write.csv(chr_summary, output_csv, row.names = FALSE)
  }

  return(chr_locations)
}
