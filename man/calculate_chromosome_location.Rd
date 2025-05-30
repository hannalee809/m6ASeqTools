% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_chromosome_location.R
\name{calculate_chromosome_location}
\alias{calculate_chromosome_location}
\title{Calculate Chromosome Location}
\usage{
calculate_chromosome_location(
  data.site_proba_df,
  output_csv = NULL,
  output_plot = NULL
)
}
\arguments{
\item{data.site_proba_df}{m6anet output dataframe containing at least a column 'ensembl_gene_name'.}

\item{output_csv}{Optional path to save chromosome summary CSV.}

\item{output_plot}{Optional path to save barplot of chromosome distribution.}
}
\value{
chr_locations Merged data frame with chromosome info.
}
\description{
Adds chromosome annotation to a data frame of high-confidence m6A sites using \code{org.Hs.eg.db}.
}
