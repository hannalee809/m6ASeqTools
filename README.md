# m6ASeqTools

An R package for post-processing and interpreting site-level m6A predictions from **m6Anet**. It summarizes m6A distribution across genes, transcript biotypes, and regions, and compares conditions using a weighted modification ratio. Integration with gene expression data reveals biotype- and region-specific relationships between m6A localization and transcriptional regulation.

---

## Installation

```r
# Install directly from GitHub
devtools::install_github("hannalee809/m6ASeqTools")
```

---

## Features

### Descriptive Statistics
- Generate quality control reports from m6Anet outputs
  - Number of m6A sites/transcripts/genes
  - single- (1 m6A site) and multi- (>1 m6A site) per transcript distribution
  - modified kmer distribution
  - modified biotype distribution
  - modified transcript length distribution
- Map m6A sites to 5′UTR, CDS, and 3′UTR using GTF file
- Chromosomal distribution
- Filter sites based on probability_modified
- Calculate **weighted_modification_ratio** to summarize site-level methylation levels to the transcript

### Comparative analysis across experimental groups
- Compare common/unique methylated genes, transcripts, and sites across groups
- Integrate **gene expression** (DESeq2) for combined analysis between m6A and Gene expression
- Analyze biotype distributions and m6A site localizations for clusters of methylation and expression

---

## Usage

### Part 1: Descriptive Statistics

**m6A_seq_part1**: A function that does all descriptive statistics functions at once. Input the m6Anet output data frame, the path to the GTF file, the path to the output directory to store results, a name for output files, and the threshold for modification sites (default = 0.9). Specific functions can be run separately and usage for each are listed below. 
```r
m6A_seq_part1(m6Aout, gtf_path, output_directory, mod_ratio_df, probability = 0.9) 
```
**Summarized report of the m6anet output:** Input the *data.site_proba.csv* data frame generated from m6anet and the name for your output file and the function will output a HTML QC summary report.
```r
summarize_m6anet_output(m6Aout, output_file="qc_report.html")
```
**Filter m6A sites:** Input the m6a data frame output and then the probability modified cut-off which by default is >= 0.9. We recommend saving this dataframe to your working environment to use to subsequent analysis, since these only include high confidence m6A sites.
```r
filter_m6a_sites(m6Aout, prob_modified=0.9)
```
**Transcript region annotation and distribution:**
- First, input the path to the reference GTF file (i.e. Ensembl annotations) and the function will output a data frame with transcript IDs and associated UTR and CDS lengths.
```r
get_transcript_region_lengths("/Downloads/"reference.gtf")
```
- Input the m6A data frame (from "filter_m6a_sites") and the data frame (from "get_transcript_region_lengths") and the function will output a data frame with the information on whether the site is located in the 5'UTR, CDS, 3'UTR, or NA and the relative position of the site on that region. This output will be necessary for part 2 of m6ASeqTools to compare transcript regions between experimental conditions.
```r
map_relative_tx_regions_to_m6A(m6Aout,tx_regions)
```
- Optionally, you can input the data frame generated from "map_relative_tx_regions_to_m6A" output and a string for the plot label and the function will return a plot showing KDE density estimates of m6A sites across the transcript.
```r
plot_rel_positions(df,"label")
```
**Chromosomal annotation and distribution:**
Input the filtered m6a output and then optional filepaths to save the output CSV and plot. Function will return a CSV file containing percentage of each chromosome and a barplot showing chromosomal distribution.
```r
calculate_chromosome_location(m6Aout,output_csv=NULL, output_plot=NULL)
```
**weighted modification ratio calculation:**
Input m6a data frame (needs to contain transcript region analysis from **map_relative_tx_regions_to_m6A**) and a string with what the user would like to name the modified ratio column and the function will return a data frame containing a column with weighted mod ratio values.
```r
calculate_weighted_mod_ratio(m6a_df, mod_ratio_column = "mod_ratio")
```

### Part 2: Comparative Analysis between Conditions

**m6a_seq_part2** An integrated function for a comparative analysis between two groups
Input a list of vectors containing m6A sites, transcript, or genes if the user would like to use the compare_m6A_distrubition function (optional), the data frame generated from m6A_seq_part1 for one group, the resulting data frame from m6A_seq_part1 for another group, a string for group 1, and a string for group 2, DESEQ resulting data frame, the column name from DESEQ for gene ID, the DESEQ column name for log2 fold change, the DESEQ column name for the adjusted p-value, and the path to the directory to save resulting outputs.
```r
m6a_seq_part2(vec_list, g1_df, g2_df, group1_name = "Group1", group2_name = "Group2", deseq, gene_col = "ensembl_gene_id", log2fc_col = "log2FoldChange", padj_col = "padj", output_dir)
```
**Comparing unique and common sites/transcripts/genes across conditions:** Input a list of vectors containing gene names, transcript ids, or site ids for every group and the function will return a list of elements unique to each group, pairwise overlaps, elements shared by all groups, and each pair-wise overlap excluding the third group.
```r
gene_list <- list(
  sample1 = unique(sample1_m6a$ensembl_gene_name),
  sample2 = unique(sample2_m6a$ensembl_gene_name),
  sample3 = unique(sample3_m6a$ensembl_gene_name)
)
site_list <- list(
  sample1=unique(paste(sample1_m6a$ensembl_transcript_id, sample1_m6a$transcript_position, sep="_")),
  sample2=unique(paste(sample2_m6a$ensembl_transcript_id, sample2_m6a$transcript_position, sep="_")),
  sample3=unique(paste(sample3_m6a$ensembl_transcript_id, sample3_m6a$transcript_position, sep="_")),
)
compare_m6a_distribution(gene_list)
compare_m6a_distribution(site_list)
```
**Differential methylation levels between two groups using weighted modification ratios:**
Input two dataframes, one for group1 (experimental group) and group2 (control group), from the **calculate_weighted_mod_ratio** along with column names for group1 and group2. and the function will return a data frame with the merged columns and a column for log2 of ratio of the weighted mod ratio for group 1 over group 2.
```r
log2fc_weighted_mod_ratio(group1,group2,group1_name="Group 1",group2_name="Group 2")
```
**Combine differential methylation with differential gene expression**
Input output from **log2fc_weighted_mod_ratio**, DESEQ dataframe with log2 fold change and p-values. Specify the gene column name, log2fc column name, and pvalue name in the DESEQ dataframe. Specify group 1 (experimental group) and group 2 (control group) from **log2fc_weighted_mod_ratio.** **Make sure that the group order matches in the differential methylation and expression outputs, meaning that the experimental and control groups are identical!**
```r
weighted_mod_ratio_and_DGE(log2fc_df, deseq_df, gene_col = "gene_id", log2fc_col = "log2fc", padj_col = "padj", group1_name = "Group1", group2_name = "Group2")
```
**Summarize clusters of methylation and expression and summarize m6A localization and biotypes for each cluster**
Input the data frame generated from **weighted_mod_ratio_and_DGE**, the group names, the name of the log2 fold change column from DESEQ, the name of the log2 fold change column from weighted mod ratio, the name of the significance column containing boolean values, and whether you want a combined summary.  
```r
summarize_weighted_mod_ratio_and_DGE(df,group1_name = "Group1",group2_name = "Group2", log2FC_dge, log2fc_wmr, sig_col = "significant", combine = TRUE)
```
---
