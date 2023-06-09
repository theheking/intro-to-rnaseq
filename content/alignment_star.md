---
layout: page
title: Extension) Understanding and Using STAR
subtitle: Adapted from Meeta Mistry, Bob Freeman and Mary Piper
---

> Overview
> --------
> **Questions**
> 
> *  How do I perform pseudo-alignment to map the transcriptome of my sample?
>     
> 
> **Objectives**
> 
> * Understanding the alignment method STAR utilizes to align sequence reads to the reference genome
>     
> * Identifying the intricacies of alignment tools used in NGS analysis (parameters, usage, etc)
>     
> *  Submit your job to the cluster



## Read Alignment
==================

<img src="../assets/img/RNAseqWorkflow.png" width="400">

Now that we have explored the quality of our raw reads, we can move on to read alignment. We perform read alignment or mapping to determine where in the genome the reads originated from. The alignment process consists of choosing an appropriate reference genome to map our reads against and performing the read alignment using one of several splice-aware alignment tools such as [STAR](http://bioinformatics.oxfordjournals.org/content/early/2012/10/25/bioinformatics.bts635) or [HISAT2](http://ccb.jhu.edu/software/hisat2/index.shtml). The choice of aligner is often a personal preference and also dependent on the computational resources that are available to you.

## STAR Aligner

To determine where on the human genome our reads originated from, we will align our reads to the reference genome using [STAR](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/) (Spliced Transcripts Alignment to a Reference). STAR is an aligner designed to specifically address many of the challenges of RNA-seq data mapping using a strategy to account for spliced alignments. 

### STAR Alignment Strategy

STAR is shown to have high accuracy and outperforms other aligners by more than a factor of 50 in mapping speed, but it is memory intensive. The algorithm achieves this highly efficient mapping by performing a two-step process:

1. Seed searching
2. Clustering, stitching, and scoring

#### Seed searching

For every read that STAR aligns, STAR will search for the longest sequence that exactly matches one or more locations on the reference genome. These longest matching sequences are called the Maximal Mappable Prefixes (MMPs):

<img src="../assets/img/alignment_STAR_step1.png" width="400">

	
The different parts of the read that are mapped separately are called 'seeds'. So the first MMP that is mapped to the genome is called *seed1*.

STAR will then search again for only the unmapped portion of the read to find the next longest sequence that exactly matches the reference genome, or the next MMP, which will be *seed2*. 

<img src="../assets/img/alignment_STAR_step2.png" width="400">

This sequential searching of only the unmapped portions of reads underlies the efficiency of the STAR algorithm. STAR uses an uncompressed suffix array (SA) to efficiently search for the MMPs, this allows for quick searching against even the largest reference genomes. Other slower aligners use algorithms that often search for the entire read sequence before splitting reads and performing iterative rounds of mapping.

**If STAR does not find an exact matching sequence** for each part of the read due to mismatches or indels, the previous MMPs will be extended.

<img src="../assets/img/alignment_STAR_step3.png" width="400">

**If extension does not give a good alignment**, then the poor quality or adapter sequence (or other contaminating sequence) will be soft clipped.

<img src="../assets/img/alignment_STAR_step4.png" width="400">


#### Clustering, stitching, and scoring

The separate seeds are stitched together to create a complete read by first clustering the seeds together based on proximity to a set of 'anchor' seeds, or seeds that are not multi-mapping.

Then the seeds are stitched together based on the best alignment for the read (scoring based on mismatches, indels, gaps, etc.). 

<img src="../assets/img/alignment_STAR_step5.png" width="400">


## Running STAR

### Set-up

To get started with this lesson, start an interactive session with 6 cores:

```
$ qrsh -pe smp 6
```

You should have a directory tree setup similar to that shown below. it is best practice to have all files you intend on using for your workflow present within the same directory. In our case, we have our original FASTQ files generated in the previous section. 

To use the STAR aligner, load the module: 

```
	$ export MODULEPATH=/share/ClusterShare/Modules/modulefiles/contrib/centos7.8:$MODULEPATH
	$ module load centos7.8/joaach/STAR/2.7.7a
```

Aligning reads using STAR is a two-step process:   

1. Create a genome index 
2. Map reads to the genome


### Creating a genome index

For this workshop we are using reads that originate from a small subsection of chromosome 1 (~300,000 reads) and so we are using only chr1 as the reference genome. 

To store our genome indices, we will use the `/n/scratch2/` space with large temporary storage capacity. We need to create a directory for the indices within this space:

```bash
	mkdir -p /share/ScratchGeneral/[your_ID]/rnaseq_tutorial/star_human_ref/
	cd /share/ScratchGeneral/[your_ID]/rnaseq_tutorial/star_human_ref/
	wget https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
	wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.chr.gtf.gz

```

The basic options to **generate genome indices** using STAR are as follows:

* `--runThreadN`: number of threads
* `--runMode`: genomeGenerate mode
* `--genomeDir`: /path/to/store/genome_indices
* `--genomeFastaFiles`: /path/to/FASTA_file 
* `--sjdbGTFfile`: /path/to/GTF_file
* `--sjdbOverhang`: readlength -1

> *NOTE:* In case of reads of varying length, the ideal value for `--sjdbOverhang` is max(ReadLength)-1. In most cases, the default value of 100 will work similarly to the ideal value.

Now let's create a job submission script to generate the genome index:

```bash
	$ mkdir /share/ScratchGeneral/[your_ID]/scripts/
	$ vim /share/ScratchGeneral/[your_ID]/scripts/STAR_index.run
```
Within `vim` we now add our shebang line, the SLURM directives, and our STAR command. 

```bash
#$ -S /bin/sh
#$ -q short.q
#$ -pe smp 22
#$ -j y
#$ -b y

#$ -cwd
#making sure bashprofile is loaded -this depends on whether this is in your /home/user/ folder
#. ~/.bash_profile
#loading module path for setting up environment within qsub job
export MODULEPATH=/share/ClusterShare/Modules/modulefiles/contrib/centos7.8:$MODULEPATH

RUN_MODE="genomeGenerate"
GENOME_DIR="/share/ScratchGeneral/[your_ID]/rnaseq_tutorial/STAR/"
GENOME_FASTA="/share/ClusterShare/biodata/contrib/[your_ID]/GRCh38_STAR_index_ENSEMBL/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF_FILE="/share/ClusterShare/biodata/contrib/[your_ID]/Hg38_STAR_index_ENSEMBL/Homo_sapiens.GRCh38.102.chr.gtf"
STAR --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir ${genome_Dir} \
  --genomeFastaFiles ${GEN_FASTA} \
  --sjdbGTFfile ${GEN_GTF}


```

```bash
	qsub ~/unix_lesson/rnaseq/scripts/STAR_index.run
```

### Aligning reads

Create an output directory for our alignment files:

```bash
	mkdir -p /share/ScratchGeneral/[your_ID]/STAR_output/
	cd /share/ScratchGeneral/[your_ID]/STAR_output/
```

### STAR command in interactive bash

Details on STAR and its functionality can be found in the [user manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf); we encourage you to peruse through to get familiar with all available options.

The basic options for aligning reads to the genome using STAR are:

* `--runThreadN`: number of threads / cores
* `--readFilesIn`: /path/to/FASTQ_file
* `--genomeDir`: /path/to/genome_indices_directory
* `--outFileNamePrefix`: prefix for all output files

Listed below are additional parameters that we will use in our command:

* `--outSAMtype`: output filetype (SAM default)
* `--outSAMunmapped`: what to do with unmapped reads

> **NOTE:** Default filtering is applied in which the maximum number of multiple alignments allowed for a read is set to 10. If a read exceeds this number there is no alignment output. To change the default you can use `--outFilterMultimapNmax`, but for this lesson we will leave it as default. Also, note that "**STAR’s default parameters are optimized for mammalian genomes.** Other species may require significant modifications of some alignment parameters; in particular, the maximum and minimum intron sizes have to be reduced for organisms with smaller introns" [[1](http://bioinformatics.oxfordjournals.org/content/early/2012/10/25/bioinformatics.bts635.full.pdf+html)].

We can access the software by simply using the STAR command followed by the basic parameters described above and any additional parameters. The full command is provided below for you to copy paste into your terminal. If you want to manually enter the command, it is advisable to first type out the full command in a text editor (i.e. [Sublime Text](http://www.sublimetext.com/) or [Notepad++](https://notepad-plus-plus.org/)) on your local machine and then copy paste into the terminal. This will make it easier to catch typos and make appropriate changes. 

```bash
RUN_MODE="alignReads"
SAMPLE="/share/ScratchGeneral/[your_ID]/rnaseq_tutorial/TRIMMED_FASTA/SRR306844chr1_chr3.trimmed.fastq.gz"
OUT_PREFIX="SRR306844chr1_chr3"
SAM_TYPE="BAM SortedByCoordinate"
GENOME_DIR="/share/ScratchGeneral/[your_ID]/rnaseq_tutorial/STAR/"

STAR --runMode "${RUN_MODE}" \
--genomeDir ${GENOME_DIR} \
--runThreadN 6 \
--readFilesIn ${SAMPLE} \
--outFileNamePrefix ${OUT_PREFIX} \
--outSAMtype ${SAM_TYPE}

```

> Exercise
> =========
> Automate the command above using a for loop across every sample. Hint: look at the trimmomatic for loop

---
*Adapted from a lesson  developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
