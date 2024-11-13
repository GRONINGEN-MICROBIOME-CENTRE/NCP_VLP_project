
# Negativeome in Early-Life Virome Studies: Characterization and Decontamination

This repository contains the scripts used during the data analysis for the paper **“Negativeome in Early-Life Virome Studies: Characterization and Decontamination”**.

You can download the current version of the paper from [BioRxiv](https://www.biorxiv.org/content/10.1101/2024.10.14.618243v1).

## Short Summary

This study aims to evaluate contamination patterns in viral-like particle (VLP) sequencing data by analyzing publicly available samples and negative controls (NCs) from various mother-infant cohorts. It also proposes a decontamination strategy to effectively remove contaminants from the data. Four studies were selected for this analysis:

1. **Garmaeva, Sanzhima et al.**  
   *“Transmission and dynamics of mother-infant gut viruses during pregnancy and early life.”*  
   Nature Communications, vol. 15, no. 1, 1945, 2 Mar. 2024.  
   [doi:10.1038/s41467-024-45257-4](https://doi.org/10.1038/s41467-024-45257-4)
   
2. **Liang, Guanxiang et al.**  
   *“The stepwise assembly of the neonatal virome is modulated by breastfeeding.”*  
   Nature, vol. 581, no. 7809, 2020, pp. 470-474.  
   [doi:10.1038/s41586-020-2192-1](https://doi.org/10.1038/s41586-020-2192-1)
   
3. **Maqsood, Rabia et al.**  
   *“Discordant transmission of bacteria and viruses from mothers to babies at birth.”*  
   Microbiome, vol. 7, no. 1, 156, 10 Dec. 2019.  
   [doi:10.1186/s40168-019-0766-7](https://doi.org/10.1186/s40168-019-0766-7)
   
4. **Shah, Shiraz A et al.**  
   *“Expanding known viral diversity in the healthy infant gut.”*  
   Nature Microbiology, vol. 8, no. 5, 2023, pp. 986-998.  
   [doi:10.1038/s41564-023-01345-7](https://doi.org/10.1038/s41564-023-01345-7)

Details regarding metadata and the study setups can be found in the original publications cited above.

## Scripts Overview

The scripts in this repository are divided into three categories:

### 1. Upstream Analysis
- **Location:** [Upstream Analysis Folder](https://github.com/GRONINGEN-MICROBIOME-CENTRE/NCP_VLP_project/tree/master/Upstream_analysis)
- This section contains scripts used for raw data processing up to the creation of RPKM tables. The scripts are organized into four study-specific folders, named after the first authors of the publications used in this study. 
- Additionally, the folder named "NCP_all" includes scripts for pooled data analysis, including dereplication, decontamination, and RPKM table creation.
- Further information can be found in the README files provided within each folder.

### 2. Downstream Analysis
- **Location:** [Downstream Analysis Folder](https://github.com/GRONINGEN-MICROBIOME-CENTRE/NCP_VLP_project/tree/master/Downstream_analysis)
- This section contains scripts used for statistical analysis and plots creation.

### 3. Metadata processing
- **Location:** [Metadata Folder](https://github.com/GRONINGEN-MICROBIOME-CENTRE/NCP_VLP_project/tree/master/Metadata_processing)
- This section contains scripts used for metadata creation.

## Associated data

The datasets analyzed in this study are publicly accessible. Sequencing data from Maqsood et al., Liang et al., and Shah et al. can be found on the European Nucleotide Archive under project identifiers PRJEB33578, PRJNA524703, and PRJEB46943, respectively. Sequencing data for Garmaeva et al. is hosted in the European Genome-Phenome Archive (EGA) under study ID EGAS00001005969. Relevant source data, such as redundant viral sequences, vOTU representatives, associated metadata, and a virus sequence database from negative control samples (version 1.0.0), are available through the [FigShare repository](https://doi.org/10.6084/m9.figshare.27170739).

While we recommend using the complete dataset available at the [FigShare repository](https://doi.org/10.6084/m9.figshare.27170739) to test downstream analysis scripts and ensure reproducibility, any subset of samples from the original studies can be used for testing upstream analysis scripts. Rerunning the downstream analysis script should yield the same results and figures as in the original publication. For the upstream analysis, the expected output may vary depending on the specific script and sample set; however, the final results should include the same set of viruses identified in the selected samples, along with their corresponding metadata. The final set will be available in the "All_identified_virus_sequences" folder. More information can be found within the scripts and README files.

## System requirements

A complete list of tools and packages, along with their respective versions used in this study, is available in the methods section of the original publication currently available at the [BioRxiv](https://www.biorxiv.org/content/10.1101/2024.10.14.618243v1).
**Note:** All analyses were conducted on the Gearshift and Hábrók high-performance computing clusters, provided by the Genomics Coordination Center and the Center for Information Technology at the University of Groningen. While some portions of the analysis can be executed on a standard laptop, others are computationally intensive and may require additional time and memory, making high-performance computing resources preferable for efficient processing.

## Installation and run time
Before running the provided code, ensure that all code dependencies are installed following each tool or package's installation instructions. Installation times vary but generally complete within a day. Note that some tools require significant storage due to large database files, making high-performance computing clusters preferable.

Once dependencies are installed, all code from this study is ready for execution. Estimated runtime and memory requirements are specified in the headers of the .sh scripts used in the upstream analysis. The metadata creation script has an expected runtime of under 10 minutes, while the downstream analysis typically completes within an hour.
