#' Example m6A sites dataset: Gli36 IGF2BP2 knockdown
#'
#' A small example dataset of predicted m6A sites for the Gli36 glioma cell line
#' after IGF2BP2 siRNA knockdown. This dataset is used to illustrate the
#' workflow in vignettes.
#'
#' @format A data frame with 100 rows and 13 variables:
#' \describe{
#'   \item{ensembl_transcript_id}{Ensembl transcript ID}
#'   \item{ensembl_gene_id}{Ensembl gene ID}
#'   \item{otthumg_gene_id}{otthumg gene ID}
#'   \item{otthumg_transcript_id}{otthumg transcript ID}
#'   \item{ensembl_transcript_name}{Ensembl transcript name}
#'   \item{ensembl_gene_name}{Ensembl gene name}
#'   \item{transcript_length}{Transcript length}
#'   \item{gene_biotype}{Transcript biotype}
#'   \item{transcript_position}{Transcript position}
#'   \item{n_reads}{Number of reads for this site}
#'   \item{probability_modified}{Probability of modification}
#'   \item{kmer}{Modified kmer}
#'   \item{mod_ratio}{Modification ratio}
#' }
#' @source Generated for package examples
"example_m6A_igf_gli36_cells"
