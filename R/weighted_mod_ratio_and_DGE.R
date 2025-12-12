#' Weighted Modification Ratio vs Differential Gene Expression Plot
#'
#' This function merges weighted modification ratio (WMR) results with DESeq2
#' differential gene expression results and generates a scatter plot comparing
#' WMR changes to gene-level log2 fold changes. Users can either input the
#' results from `log2fc_weighted_mod_ratio` or `run_wmr_differential_test`,
#' which they specify in the function argument. Users may specify significance
#' thresholds for WMR if they used the differential test. (p-value and/or
#' median/log2FC difference), which highlight points in blue Genes significant
#' in both WMR and DGE are highlighted in purple.
#'
#' @param wmr_df Data frame containing weighted modification ratio metrics.
#'               Must include either \code{log2fc_weighted_mod_ratio} or
#'               \code{median_diff} depending on \code{input_type}.
#' @param deseq DESeq2 results data frame containing gene-level differential
#'              expression statistics.
#' @param input_type Character. Choose either:
#'        \itemize{
#'          \item \code{"log2fc_wmr"} – use log2FC of weighted mod ratio.
#'          \item \code{"wmr_diff"} – use median WMR difference.
#'        }
#' @param wmr_gene_col Column name in \code{wmr_df} storing gene IDs.
#' @param wmr_sig List specifying optional WMR significance thresholds for
#'   `input_type = "wmr_diff`:
#'        \itemize{
#'          \item \code{p_val}: p-value cutoff for WMR significance.
#'          \item \code{diff}: minimum absolute log2FC or median difference.
#'        }
#'        Use \code{NULL} to disable a criterion.
#' @param deseq_gene_col Column name in \code{deseq} storing gene IDs.
#' @param log2fc_col Column name in \code{deseq} for gene log2 fold change.
#' @param padj_col Column name in \code{deseq} for adjusted p-values.
#' @param group1_name Name of condition/group 1 (used in labels).
#' @param group2_name Name of condition/group 2 (used in labels).
#'
#' @return A list containing:
#' \describe{
#'   \item{plot}{A ggplot2 scatter plot.}
#'   \item{data}{Merged data frame containing WMR + DGE values.}
#' }
#'
#' @details
#' Significance coloring:
#' \itemize{
#'   \item Purple = significant in both WMR and DGE.
#'   \item Red = DGE only.
#'   \item Blue = WMR only.
#'   \item Black = not significant.
#' }
#'
#' DGE significance is defined as:
#' \code{|log2FC| > 0.58} and \code{padj < 0.05}.
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_identity geom_hline theme_minimal labs geom_vline
#' @importFrom dplyr case_when
#'
#' @examples
#' weighted_mod_ratio_and_DGE_NEW(
#'   wmr_df,
#'   deseq_results,
#'   input_type = "wmr_diff",
#'   wmr_sig = list(p_val = 0.05, diff = 0.2)
#' )
#'
#' @export
weighted_mod_ratio_and_DGE <- function(
  wmr_df,
  deseq,
  input_type = c("log2fc_wmr", "wmr_diff"),
  wmr_gene_col = "ensembl_gene_id",
  wmr_sig = list(p_val = NULL, diff = NULL),
  deseq_gene_col = "ensembl_gene_id",
  log2fc_col = "log2FoldChange",
  padj_col = "padj",
  group1_name = "Group1",
  group2_name = "Group2"
) {

  input_type <- match.arg(input_type)

  # ---------- 1. Merge ----------
  merged <- merge(wmr_df, deseq, by.x = wmr_gene_col, by.y = deseq_gene_col, all.x = TRUE)

  # ---------- 2. X-axis ----------
  if (input_type == "log2fc_wmr") {
    if (!"log2fc_weighted_mod_ratio" %in% colnames(merged))
      stop("log2fc_weighted_mod_ratio missing from wmr_df")

    merged$x_value <- merged$log2fc_weighted_mod_ratio
    x_label <- paste("Log2FC (Weighted Mod Ratio):", group1_name, "vs", group2_name)

  } else if (input_type == "wmr_diff") {
    if (!"median_diff" %in% colnames(merged))
      stop("median_diff missing from wmr_df")

    merged$x_value <- merged$median_diff
    x_label <- paste("Median Difference (WMR):", group2_name, "-", group1_name)
  }

  # ---------- 3. DGE ----------
  merged$gene_log2fc <- merged[[log2fc_col]]
  merged$gene_padj <- merged[[padj_col]]

  merged$dge_sig <- abs(merged$gene_log2fc) > 0.58 & merged$gene_padj < 0.05

  # ---------- 4. WMR significance ----------
  wmr_p <- wmr_sig$p_val
  wmr_diff <- wmr_sig$diff

  merged$wmr_sig <- TRUE

  if (!is.null(wmr_p) && "p_value" %in% colnames(merged))
    merged$wmr_sig <- merged$wmr_sig & merged$p_value < wmr_p

  if (!is.null(wmr_diff))
    merged$wmr_sig <- merged$wmr_sig & abs(merged$x_value) > wmr_diff

  if (is.null(wmr_p) && is.null(wmr_diff))
    merged$wmr_sig <- FALSE

  # ---------- 5. Colors ----------
  merged$point_color <- dplyr::case_when(
    merged$dge_sig & merged$wmr_sig ~ "purple",
    merged$dge_sig & !merged$wmr_sig ~ "red",
    !merged$dge_sig & merged$wmr_sig ~ "blue",
    TRUE ~ "black"
  )

  # ---------- 6. Title ----------
  title_text <- paste0("WMR (", group1_name, " vs ", group2_name, ") vs Gene Expression")

  # ---------- 7. Plot ----------
  p <- ggplot2::ggplot(
    merged,
    ggplot2::aes(x = x_value, y = gene_log2fc)
  ) +
    ggplot2::geom_point(ggplot2::aes(color = point_color), alpha = 0.7) +
    ggplot2::scale_color_identity(
      guide = "legend",
      labels = c(
        "purple" = "DGE + WMR Significant",
        "red"    = "DGE Significant",
        "blue"  = "WMR Significant",
        "black"  = "Not Significant"
      ),
      breaks = c("purple", "red", "blue", "black")
    ) +
    ggplot2::geom_vline(
      xintercept = c(0,0),
      colour = "gray",
      linetype = "dashed",
      linewidth = 0.5
    ) +
    ggplot2::geom_hline(
      yintercept = c(-0.58, 0.58),
      colour = "gray",
      linetype = "dashed",
      linewidth = 0.5
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = title_text,
      x = x_label,
      y = paste("Log2FC (", log2fc_col, ")"),
      color = "Significance"
    )

  return(list(plot = p, data = merged))
}



