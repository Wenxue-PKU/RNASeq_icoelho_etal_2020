# Repository for the RNASeq Analysis of the Paper "Trem-2 promotes emergence of restorative macrophages and endothelial cells during recovery from hepatic tissue damage" (Coelho et al., 2020)

This repository contains the raw data (counts files) and scripts for the analysis of the RNASeq for the referenced paper.

##### Quality Control

This quality control is done using the program **Fastqc** (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/); considering this program produces a report per file, we can also use a program to merge all the reports into one - **MultiQC** (https://multiqc.info/).

The reports for *the files corresponding to the 4 lanes separated* are saved in *fastqc_reports_perlane*. A summary of all reports is under the name *mQC_raw_perlane*.

The reports for *the files corresponding to the 4 lanes merged* are saved in *fastqc_reports_merged*. A summary of all reports is under the name *mQC_raw_merged*.

Looking at the overall reports, the samples look like to have a good overall quality, sequencing depth and profile (some of the plots will give out warnings, such as in %GC and Overrepresented sequences: these warnings are not uncommon in RNA-Seq datasets but, are not very worrying).
