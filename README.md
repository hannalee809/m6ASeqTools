# m6ASeqTools

An R package for downstream analysis of m6A RNA modifications detected by **m6Anet**.  
Includes transcript-level summaries, regional distributions, QC report generation, and comparative analyses across groups.

---

## Features

### Descriptive Statistics
- Generate quality control reports from m6Anet outputs (summarize_m6anet_output)
- Map m6A to 5′UTR, CDS, and 3′UTR using transcript models (get_transcript_region_legnths, map_relative_tx_regions_to_m6A,
plot_relative_positions)
- chromosome location (calculate_chromosome_location)
- filter sites (filter_m6a_sites)
- overall transcript level methylation (calculate_weighted_mod_ratio)

### Comparative analysis across experimental groups
- Compare common/unique methylated genes, transcripts, and sites across groups (compare_m6a_distribution)
- Integrate **gene expression** (DESeq2) for combined analysis (weighted_mod_ratio_and_DGE)

---

## 📦 Installation

```r
# Install directly from GitHub
devtools::install_github("YOUR_USERNAME/m6ASeqTools")
