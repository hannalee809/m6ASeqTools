#' Compare m6A Distribution between Samples
#'
#' Compares common and unique m6A sites, transcripts, or genes for 2 or more groups
#'
#' @param named_list list of vectors (gene names, transcript ids, or site ids per group)
#'
#' @return a list of unique elements per group, pairwise intersections, common elements
#' to all, pairwise intersections excluding the third (if exactly 3 groups are input)
#'
#' @importFrom utils combn
#' @export
compare_m6a_distribution <- function(named_list) {
  group_names <- names(named_list)
  if (is.null(group_names)) stop("Input list must be named.")

  # Unique to each group (not found in any of the others)
  unique_per_group <- lapply(seq_along(named_list), function(i) {
    others <- unlist(named_list[-i])
    setdiff(named_list[[i]], others)
  })
  names(unique_per_group) <- paste0("unique_", group_names)

  # Pairwise intersections
  pairwise <- utils::combn(group_names, 2, simplify = FALSE)
  pairwise_inters <- lapply(pairwise, function(groups) {
    intersect(named_list[[groups[1]]], named_list[[groups[2]]])
  })
  names(pairwise_inters) <- sapply(pairwise, function(x) paste(x, collapse = "_"))

  # Common across all
  common_all <- Reduce(intersect, named_list)

  result <- c(unique_per_group, pairwise_inters, list(common_all = common_all))

  # If exactly 3 groups, add "only in pair but not the third"
  if (length(group_names) == 3) {
    exclusive_pairs <- list()
    for (combo in pairwise) {
      third <- setdiff(group_names, combo)
      label <- paste0(paste(combo, collapse = "_"), "_only")
      exclusive_pairs[[label]] <- setdiff(
        intersect(named_list[[combo[1]]], named_list[[combo[2]]]),
        named_list[[third]]
      )
    }
    result <- c(result, exclusive_pairs)
  }

  # Also return lengths
  summary_counts <- lapply(result, length)
  attr(result, "summary") <- summary_counts

  return(result)
}


