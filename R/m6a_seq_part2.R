#' Weighted Modification Ratio and Differential Gene Expression Analysis
#'
#' This function compares m6A methylation patterns across two groups, calculates
#' the log2 fold change of weighted modification ratios, integrates differential
#' gene expression (DGE) results, and summarizes findings including transcript
#' biotype and region composition.
#'
#' @param vec_list A named list of vectors containing m6A sites, transcripts, or
#'   genes for comparison (used in `compare_m6a_distribution`). This step is
#'   optional.
#' @param g1_df Data frame containing weighted modification ratios for group 1
#' @param g2_df Data frame containing weighted modification ratios for group 2
#' @param group1_name Name for group 1 (default: "Group1")
#' @param group2_name Name for group 2 (default: "Group2")
#' @param deseq DESeq2 results data frame (must include gene ID, log2FC, and
#'   padj columns)
#' @param gene_col Column name in DESeq2 results for gene ID (default:
#'   "ensembl_gene_id")
#' @param log2fc_col Column name in DESeq2 results for log2 fold change
#'   (default: "log2FoldChange")
#' @param padj_col Column name in DESeq2 results for adjusted p-value (default:
#'   "padj")
#' @param output_dir Directory where output files will be saved
#'
#' @return Writes several CSV and PDF files to the output directory and returns
#'   invisibly
#' @export
m6a_seq_part2 <- function(vec_list,
                        g1_df,
                        g2_df,
                        group1_name = "Group1",
                        group2_name = "Group2",
                        deseq,
                        gene_col = "ensembl_gene_id",
                        log2fc_col = "log2FoldChange",
                        padj_col = "padj",
                        output_dir) {

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Set working directory temporarily
  old_dir <- getwd()
  setwd(output_dir)
  on.exit(setwd(old_dir), add = TRUE)

  # Optional: Compare m6A distributions
  if (compare_sites && !is.null(vec_list)) {
    comparison <- compare_m6a_distribution(vec_list)
    new_compare <- lapply(comparison, function(x) {
      length(x) <- max(lengths(comparison))
      return(x)
    })
    write.csv(do.call(cbind, new_compare), "compare_m6A_distribution.csv", row.names = FALSE)
  }

  # 2. Calculate log2FC of weighted mod ratio
  log2fc_df <- log2fc_weighted_mod_ratio(g1_df, g2_df, group1_name, group2_name)
  write.csv(log2fc_df, "log2fc_weighted_mod_ratio.csv", row.names = FALSE)

  # 3. Merge with DGE and plot
  wmr <- weighted_mod_ratio_and_DGE(
    log2df = log2fc_df,
    deseq = deseq,
    gene_col = gene_col,
    log2fc_col = log2fc_col,
    padj_col = padj_col,
    group1_name = group1_name,
    group2_name = group2_name
  )
  ggsave("weighted_mod_ratio_and_DGE_plot.pdf", plot = wmr$plot)
  write.csv(wmr$data, "weighted_mod_ratio_and_DGE.csv", row.names = FALSE)

  # 4. Summarize biotype + region by expression clusters
  sum_wmr_dge <- summarize_weighted_mod_ratio_and_DGE(
    df = wmr$data,
    group1_name = group1_name,
    group2_name = group2_name,
    log2FC_dge = paste0(group1_name, "_vs_", group2_name, "_log2FC"),
    log2fc_wmr = "log2fc_weighted_mod_ratio",
    sig_col = "significant",
    combine = TRUE
  )

  # Save summaries
  write.csv(sum_wmr_dge$biotype_summary_all, "summarized_weighted_mod_ratio_and_DGE_biotype_summary.csv", row.names = FALSE)
  write.csv(sum_wmr_dge$region_summary_all, "summarized_weighted_mod_ratio_and_DGE_region_summary.csv", row.names = FALSE)

  # Save group-level cluster breakdown
  capture.output(str(sum_wmr_dge$groups), file = "summarized_weighted_mod_ratio_and_DGE_groups.txt")

  invisible(NULL)
}
