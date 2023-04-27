---
layout: page
title: RNA-seq Tutorial Outline
subtitle: Temporary Overview
---


# Packages Needed
SRA Toolkit
Fastqc
multiqc 
trimmomatic 
Kallisto
R
o	Tximport
o	DESeq2
 
 
## Pipeline Steps
o	Examining data on SRA Database
o	https://datacarpentry.org/organization-genomics/03-ncbi-sra/index.html
o	Download from SRA using fastq-dump 
o	Examining Read quality 
- FASTQC for quality control metrics
- summarise FASTQC with multiqc
o	https://datacarpentry.org/wrangling-genomics/02-quality-control/index.html
o	Trimming and Filtering 
- Trim reads using trimmomatic 
o	https://datacarpentry.org/wrangling-genomics/03-trimming/index.html
o	Pseudalignment using kallisto
o	Pseudoalignment vs alignment
o	<<kallisto.pdf>>
o	<<kallisto and pseudoalignment.pdf>>
 
o	Import DESeq2 into RNA-seq 
o	https://ycl6.gitbook.io/guide-to-rna-seq-analysis/differential-expression-analysis/differential-gene-expression/dge-analysis-with-salmon-input
o	https://lashlock.github.io/compbio/R_presentation.html
 
 
o	Introduction to shell 
o	https://datacarpentry.org/shell-genomics/
o	https://linuxjourney.com/lesson/the-shell
o	Introduction to R
o	https://datacarpentry.org/genomics-r-intro/
 
 
