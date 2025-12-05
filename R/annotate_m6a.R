#' Annotate m6A Sites with External Feature Overlaps
#'
#' This function identifies overlaps between m6A sites and a user-provided
#' genomic annotation dataset (e.g., RBP binding sites, enhancers, promoters,
#' etc.). The user specifies the required genomic coordinate columns
#' (`chr`, `start`, `end`, `strand`) and the column containing the
#' annotation name. All remaining columns in the
#' annotation dataset are retained and merged into the overlap table.
#'
#' @param m6a_df A data frame containing m6A site coordinates. Must include
#'   columns for chromosome, start, end, and strand.
#' @param annot_df A data frame containing annotation genomic intervals.
#'   Must include columns for chromosome, start, end, and strand, plus an
#'   annotation name column specified in `annot_name_col`.
#' @param m6a_chr_col Name of the m6A chromosome column.
#' @param m6a_start_col Name of the m6A start coordinate column.
#' @param m6a_end_col Name of the m6A end coordinate column.
#' @param m6a_strand_col Name of the m6A strand column.
#' @param annot_chr_col Name of the annotation chromosome column.
#' @param annot_start_col Name of the annotation start coordinate column.
#' @param annot_end_col Name of the annotation end coordinate column.
#' @param annot_strand_col Name of the annotation strand column.
#' @param feature_col Name of feature column in the annotation. (i.e. RBP_name, miRNA_ID, SNP_ID, etc.)
#' @param zero_based Whether input files contain 0-based genomic coordinates. Default is TRUE.
#' @param ignore_strand Logical; if FALSE (default), overlaps must occur on the
#'   same strand. Set to TRUE to ignore strand information, useful for
#'   annotations such as SNPs or regions lacking strand data.
#'
#'#' @return A list containing two data frames:
#' \describe{
#'   \item{\code{overlap_df}}{
#'     A row-level table where each row corresponds to a single overlap event.
#'     Includes:
#'     \itemize{
#'       \item \code{m6A_index}: index of the m6A site.
#'       \item \code{feature_col}: annotation name provided by \code{annot_name_col}.
#'       \item All additional annotation metadata columns.
#'     }
#'   }
#'
#'   \item{\code{annotated_m6A_df}}{
#'     A per–m6A-site summary table created by merging the m6A input with a
#'     summarised annotation table. Contains:
#'     \itemize{
#'       \item Original m6A coordinates and metadata.
#'       \item \code{Overlapping_Features}: semicolon-separated annotation names.
#'       \item \code{Feature_Count}: number of unique annotations per m6A site.
#'     }
#'   }
#' }
#'
#' Both inputs are internally converted to \code{GenomicRanges} objects.
#' Overlaps are detected using \code{findOverlaps()}. All annotation metadata
#' columns are preserved and attached to the overlap table.
#' A warning may be issued if the chromosome naming or reference genome differs
#' between the two inputs.
#'
#' @examples
#' \dontrun{
#' res <- annotate_m6a(
#'   m6a_df,
#'   annot_df,
#'   m6a_chr_col = "chr",
#'   m6a_start_col = "start",
#'   m6a_end_col = "end",
#'   m6a_strand_col = "strand",
#'   annot_chr_col = "chrom",
#'   annot_start_col = "start",
#'   annot_end_col = "end",
#'   annot_strand_col = "strand",
#'   annot_name_col = "RBP",
#'   zero_based = TRUE,
#'   ignore_strand = FALSE
#' )
#'
#' res$overlap_df
#' res$annotated_m6A_df
#' }
#'
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom IRanges IRanges
#' @importFrom dplyr group_by summarise
#' @importFrom magrittr %>%
#' @export

annotate_m6a <- function(
  m6a_df,
  annot_df,
  m6a_chr_col = "chr",
  m6a_start_col = "start",
  m6a_end_col = "end",
  m6a_strand_col = "strand",
  annot_chr_col = "chr",
  annot_start_col = "start",
  annot_end_col = "end",
  annot_strand_col = "strand",
  feature_col = NULL,
  zero_based = TRUE,
  ignore_strand = FALSE
) {


  # 1. Validate columns
  if (!all(m6a_cols %in% names(m6a_df)))
    stop("Missing required m6A coordinate columns in m6a_df")

  if (!all(ann_cols %in% names(annotation_df)))
    stop("Missing required annotation coordinate columns in annotation_df")

  # Feature column MUST be provided
  if (is.null(feature_col) || !(feature_col %in% names(annotation_df)))
    stop("feature_col must be provided and must exist in annotation_df")


  # 2. Convert m6A to GRanges
  m6a_start <- if (zero_based) m6a_df[[m6a_cols["start"]]] + 1 else m6a_df[[m6a_cols["start"]]]

  m6a_gr <- GenomicRanges::GRanges(
    seqnames = m6a_df[[m6a_cols["chr"]]],
    ranges = IRanges::IRanges(start = m6a_start,
                     end = m6a_df[[m6a_cols["end"]]]),
    strand = if ("strand" %in% names(m6a_df)) m6a_df[[m6a_cols["strand"]]] else "*"
  )

  GenomicRanges::mcols(m6a_gr) <- m6a_df[, setdiff(names(m6a_df), m6a_cols), drop = FALSE]


  # 3. Convert annotation to GRanges
  ann_start <- if (zero_based) annotation_df[[ann_cols["start"]]] + 1 else annotation_df[[ann_cols["start"]]]

  ann_gr <- GenomicRanges::GRanges(
    seqnames = annotation_df[[ann_cols["chr"]]],
    ranges = IRanges::IRanges(start = ann_start,
                     end = annotation_df[[ann_cols["end"]]]),
    strand = if ("strand" %in% names(annotation_df)) annotation_df[[ann_cols["strand"]]] else "*"
  )

  GenomicRanges::mcols(ann_gr) <- annotation_df[, setdiff(names(annotation_df), ann_cols), drop = FALSE]


  # 4. Overlaps
  overlaps <- GenomicRanges::findOverlaps(m6a_gr, ann_gr, ignore.strand = ignore_strand)

  if (length(overlaps) == 0) {
    message("No overlaps detected.")
    return(list(overlap_df = NULL, summary_df = NULL))
  }

  m6a_idx <- S4Vectors::queryHits(overlaps)
  ann_idx <- S4Vectors::subjectHits(overlaps)


  # 5. Overlap table (full metadata)
  ann_meta <- as.data.frame(GenomicRanges::mcols(ann_gr))

  overlap_df <- data.frame(
    m6A_index = m6a_idx,
    ann_meta[ann_idx, , drop = FALSE],
    stringsAsFactors = FALSE
  )


  # 6. Summary table (per m6A site)
  summary_df <- overlap_df %>%
    dplyr::group_by(m6A_index) %>%
    dplyr::summarise(
      Overlapping_Features = paste(unique(.data[[feature_col]]), collapse = "; "),
      Feature_Count = dplyr::n_distinct(.data[[feature_col]]),
      .groups = "drop"
    )

  # 7. Merge Summary table (per m6A site) with original m6A bed file
  m6a_df$row_id <- 1:nrow(m6a_df)
  m6a_df_with_annotations <- merge(m6a_df,
                                   summary_df,
                                   by.x="row_id",
                                   by.y="m6A_index",
                                   all=TRUE)


  # 7. Return both outputs
  return(list(
    overlap_df = overlap_df,
    annotated_m6A_df = m6a_df_with_annotations
  ))
}
