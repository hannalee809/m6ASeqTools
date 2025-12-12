#' Plot Weighted Modification Ratio vs. Differential Gene Expression
#'
#' This function merges transcript-level weighted modification ratio data with differential gene expression results
#' and generates a scatter plot comparing the two, highlighting significantly differentially expressed genes.
#'
#' @param log2df A data frame containing weighted modification ratios with a column named `log2fc_weighted_mod_ratio`
#' @param deseq A DESeq2 results data frame with log2 fold changes and adjusted p-values
#' @param gene_col Character. Column name in `deseq` containing gene IDs (default: "ensembl_gene_id")
#' @param log2fc_col Character. Column name in `deseq` for log2 fold change values (default: "log2FoldChange")
#' @param padj_col Character. Column name in `deseq` for adjusted p-values (default: "padj")
#' @param group1_name Character. Name of the CONTROL group
#' @param group2_name Character. Name of the EXPERIMENTAL group
#'
#' @return A list containing:
#' \describe{
#'   \item{plot}{A \code{ggplot2} object comparing weighted mod ratio and gene expression}
#'   \item{data}{The merged data frame used to generate the plot}
#' }
#'
#' @import ggplot2
#' @export
#' @examples
#' m6a_DGE_A_vs_B <- weighted_mod_ratio_and_DGE(log2df, deseq, group1_name="A", group2_name="B")
#' df <- m6a_DGE_A_vs_B$data
#' df <- m6a_DGE_A_vs_B$plot
weighted_mod_ratio_and_DGE <- function(log2df, deseq,
                                       gene_col = "ensembl_gene_id",
                                       log2fc_col = "log2FoldChange",
                                       padj_col = "padj",
                                       group1_name = "Group1",
                                       group2_name = "Group2") {
  # Merge by gene ID
  merged <- merge(log2df, deseq, by.x = "ensembl_gene_id", by.y = gene_col)

  # Extract and evaluate columns
  merged$log2fc_value <- merged[[log2fc_col]]
  merged$padj_value <- merged[[padj_col]]
  merged$significant <- abs(merged$log2fc_value) > 0.58 & merged$padj_value < 0.05

  # Create title dynamically
  title_text <- paste("Weighted Mod Ratio (", group1_name, " vs ", group2_name, ") vs. Gene Expression", sep = "")

  # Plot
  p <- ggplot2::ggplot(merged,
                       ggplot2::aes(x = log2fc_weighted_mod_ratio,
                                    y = log2fc_value)) +
    ggplot2::geom_point(ggplot2::aes(color = significant), alpha = 0.7) +
    ggplot2::scale_color_manual(
      values = c("TRUE" = "red", "FALSE" = "black"),
      labels = c("TRUE" = "Significant", "FALSE" = "Not Significant")
    ) +
    ggplot2::geom_hline(yintercept = c(-0.58, 0.58), colour = "gray", linetype = "dashed", linewidth = 0.5) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = title_text,
      x = paste("Log2FC (Weighted Mod Ratio):", group1_name, "vs", group2_name),
      y = paste("Log2FC (", log2fc_col, ")"),
      color = "Significant DEG"
    )
  # Clean up temp helper columns if they exist
  merged$log2fc_value <- NULL
  merged$padj_value <- NULL

  return(list(plot = p, data = merged))
}


