# m6ASeqTools

An R package for downstream analysis of m6A RNA modifications detected by **m6Anet**.  
Includes transcript-level summaries, regional distributions, QC report generation, and comparative analyses across groups.

---

## 📦 Installation

```r
# Install directly from GitHub
devtools::install_github("hannalee809/m6ASeqTools")
```

---

## Features

### Descriptive Statistics
- Generate quality control reports from m6Anet outputs (summarize_m6anet_output)
- Map m6A to 5′UTR, CDS, and 3′UTR using transcript models (get_transcript_region_lengths, map_relative_tx_regions_to_m6A,
plot_relative_positions)
- chromosome location (calculate_chromosome_location)
- filter sites (filter_m6a_sites)
- overall transcript level methylation (calculate_weighted_mod_ratio)

### Comparative analysis across experimental groups
- Compare common/unique methylated genes, transcripts, and sites across groups (compare_m6a_distribution)
- Integrate **gene expression** (DESeq2) for combined analysis (weighted_mod_ratio_and_DGE)

---

## Usage

Input the data frame generated from m6A and the name for your output file and the function will output a HTML QC summary report.
```r
summarize_m6anet_output(m6Aout, output_file="qc_report.html")
```

Input the path to the GTF file of annotations and the function will output a data frame filled with transcript IDs and lengths of UTRs and CDS regions.
```r
get_transcript_region_lengths("/Downloads/"reference.gtf")
```

Input the data frame generated from "get_transcript_region_lengths" and the data frame from m6Anet and the function will output a data frame with the transcript positions for every m6A site.
```r
map_relative_tx_regions_to_m6A(m6Aout,tx_regions)
```

Input the data frame generated from "map_relative_tx_regions_to_m6A" output and a string for the label of the plot and the function will return a plot showing density of m6A sites across the transcript.
```r
plot_rel_positions(df,"label")
```

Input the m6a output dataframe and then optional paths to save the csv summary and barplot.
```r
calculate_chromosome_location(m6Aout,output_csv=NULL, output_plot=NULL)
```

Input the m6a data frame output and then the probability modified cut-off which by default is >= 0.9.
```r
filter_m6a_sites(m6Aout, prob_modified=0.9)
```

Input a m6a data frame containing site-level data and a string with what the user would like to name the modified ratio column.
```r
calculate_weighted_mod_ratio(m6a_df, mod_ratio_column = "mod_ratio")
```

Input a list of vectors containing gene names, transcript ids, or site ids for every group and the function will return a list of elements unique to each group, pairwise overlaps, elements shared by all groups, and each pair-wise overlap excluding the third group.
```r
gene_list <- list(
  sample1 = unique(sample1_m6a$ensembl_gene_name),
  sample2 = unique(sample2_m6a$ensembl_gene_name),
  sample3 = unique(sample3_m6a$ensembl_gene_name)
)

compare_m6a_distribution(gene_list)
```

Input two dataframes, one for group1 and group2, along with column names for group1 and group2 and the function will return a data frame with the merged columns and a column for log2 fold change.
```r
log2fc_weighted_mod_ratio(group1,group2,group1_name="Group 1",group2_name="Group 2")
```

Input a data frame which must have a column named "log2fc_weighted_mod_ratio", a DESEQ dataframe, multiple column names along with the name of the experimental and control group. The function will return a plot comparing weighted mod ratio with gene expression.
```r
weighted_mod_ratio_and_DGE(log2fc_df, deseq_df, gene_col = "gene_id", log2fc_col = "log2fc", padj_col = "padj", group1_name = "Group1", group2_name = "Group2")
```

---
