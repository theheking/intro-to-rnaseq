---
layout: page
title: 5) Using DESeq2 for Differential Expression Analysis
---
> Overview
> --------
> **Questions**
> 
> *   What are the top DETs in our experiment? 
>     
> 
> **Objectives**
> 
> *   Use DESeq2 to output top differentially expressed transcripts  
>  
> *   Understand the best diagrams to show differentially expressed transcripts
> 

---------------------------------------

The next step in the RNA-seq workflow is the differential expression analysis. Differential expression testing aims to determine which genes are expressed at different levels between conditions. These genes can offer biological insight into the processes affected by the condition(s) of interest.

The steps outlined in the gray box below we have already discussed, and we will now continue to describe the steps in an **end-to-end gene-level RNA-seq differential expression workflow**.

![](../assets/img/de_workflow.png)

So what does the count data actually represent? The count data used for differential expression analysis represents the number of sequence reads that originated from a particular gene. The higher the number of counts, the more reads associated with that gene, and the assumption that there was a higher gene expression level in the sample.

![](../assets/img/deseq_counts_overview.png)

***Note: We are using features that are transcripts, not genes***

Counts and CPM
---------------
Kallisto counts the number of reads that align to one transcript. This is the *raw count*; however, normalisation is needed to make accurate comparisons of gene expression between samples. Normalisation is used to account for variabilities between or within *raw counts* due to technical differences such as read depth. The default in DEGUST is *Counts per million (CPM)*. CPM accounts for sequencing depth. This is not the best normalisation method for differential expression analysis between samples. However, we are not going to learn R in this course, so we have to work with what we have. 

Using DESeq
-------------
1. Opening up a project
2. Transferring locally
3. Install all packages
4. Import Kallisto output and metadata for DESeq analysis
5. Running DESeq
6. Taking into account confounding effects **(Extension)**


Opening up a project
---------------------
- Open up RStudio.
- Under the File menu, click on New project, choose New directory, then Empty project
- Enter a name for this new folder, and choose a convenient location for it. This will be your working directory for the rest of the day.
- Confirm that the folder named in the Create project as a sub-directory of the box is where you want the working directory created. Use the Browse button to navigate folders if changes are needed.
- Click on “Create project”
- Under the Files tab on the right of the screen, click on New Folder and create a folder named data within your newly created working directory. (e.g., ~/data-carpentry/data)
- Create a new R script (`File > New File > R script`) and save it in your working directory (e.g. data-carpentry-script.R)s
- We can open it by clicking the New File button or using the `Ctrl-Shift-N` keyboard shortcut (`Cmd-Shift-N`) on Mac

**Please be thoughtful about where you are saving your directory e.g. Desktop**


Transferring to local computer
--------------------------------

Log onto the Wolfpack. Change the directory into the location that contains your aligned kallisto output `abundance.tsv`.

        $ ssh [userid]@dice02.garvan.unsw.edu.au
        $ cd /srv/scratch/zID/data/SRR306844chr1_chr3/
        $ ls abundance.tsv
      
This file contains the counts of one sample. For input into DESeq, you have to upload all the abundance.tsv files for every sample found in their respective folder.  

Logout off the cluster and stay on your laptop. 
You will now be transferring recursively downloading your files to your local computer. First move into a directory that you can access. 

        $  rsync -r [userid]@dice02.garvan.unsw.edu.au:"[your_scratch]/FASTQ_TRIMMED/*chr1_chr3/" .
        
 This should download a directory per sample containing the file called `transcript_counts.csv`. 


## Install and load packages

First, we'll need to install some add-on packages. Most generic R packages are hosted on the Comprehensive R Archive Network (CRAN, <http://cran.us.r-project.org/>). To install one of these packages, you would use `install.packages("packagename")`. You only need to install a package once, then load it each time using `library(packagename)`. Let's install the **gplots** and **calibrate** packages.

```{r install_packages, eval=FALSE}
install.packages("gplots")
install.packages("calibrate")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("tximport")
```

Bioconductor packages work a bit differently, and are not hosted on CRAN. Go to <http://bioconductor.org/> to learn more about the Bioconductor project. To use any Bioconductor package, you'll need a few "core" Bioconductor packages. Run the following commands to (1) download the installer script, and (2) install some core Bioconductor packages. You'll need internet connectivity to do this, and it'll take a few minutes, but it only needs to be done once.

```{r bioclite, eval=FALSE}
# Download the installer script
source("http://bioconductor.org/biocLite.R")

# biocLite() is the bioconductor installer function. Run it without any
# arguments to install the core packages or update any installed packages. This
# requires internet connectivity and will take some time!
biocLite()
```

To install specific packages, first download the installer script if you haven't done so, and use `biocLite("packagename")`. This only needs to be done once then you can load the package like any other package. Let's download the [DESeq2 package](http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html):

```{r load_deseq2, eval=FALSE}
# Do only once
source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
```

Now let's load the packages we'll use:

```{r load_pkgs, eval=TRUE}
library(DESeq2)
library(gplots)
library(calibrate)
```


Uploading metadata and counts table to DESeq 
-----------------------------------------------
Back to your RStudio...

Load the txImportData library to DESeq
```
        library(tximport)
```

kallisto abundance.tsv files can be imported as well, but this is typically slower than the approach above. Note that we add an additional argument in this code chunk, ignoreAfterBar=TRUE. This is because the Gencode transcripts have names like “ENST00000456328.2|ENSG00000223972.5|…”, though our tx2gene table only includes the first “ENST” identifier. We therefore want to split the incoming quantification matrix rownames at the first bar “|”, and only use this as an identifier. We didn’t use this option earlier with Salmon, because we used the argument --gencode when running Salmon, which itself does the splitting upstream of tximport. Note that ignoreTxVersion and ignoreAfterBar are only to facilitating the summarization to gene level.

```
        files <- file.path(dir, "kallisto", samples$run, "abundance.tsv.gz")
        names(files) <- paste0("sample", 1:6)
        txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
        head(txi.kallisto.tsv$counts)
```

The data.frame contains information about transcripts (one transcript per row) with the gene positions in the first five columns and then information about the number of reads aligning to the gene in each experimental sample. There are three replicates for control (column names starting with "ctl") and three for samples treated with ultraviolet-B light (starting "uvb"). We don't need the information on gene position for this analysis, just the counts for each gene and sample, so we can remove it from the data frame.


```{r remove_metadata_cols}
        # Remove first five columns (chr, start, end, strand, length)
        countdata <- countdata[ ,-(1:5)]
        head(countdata)
        colnames(countdata)
```

We can rename the columns to something shorter and a bit more readable.

```{r bad_renaming, eval=FALSE}
        # Manually
        c("ctl1", "ctl2", "ctl3", "uvb1", "uvb2", "uvb3")
```

We can do it manually, but what if we have 600 samples instead of 6? This would become cumbersome. Also, it's always a bad idea to hard-code sample phenotype information at the top of the file like this. A better way to do this is to use the `gsub` command to strip out the extra information. This more robust to introduced errors, for example if the column order changes at some point in the future or you add additional replicates.

```{r rename_cols}
        # Using gsub -- robust. Get help with ?gsub
        gsub(pattern="trimmed_|.fastq_tophat.accepted_hits.bam", replacement="", x=colnames(countdata))
        colnames(countdata) <- gsub(pattern="trimmed_|.fastq_tophat.accepted_hits.bam", replacement="", x=colnames(countdata))
        head(countdata)
```

---

**EXERCISE**

There's an R function called `rowSums()` that calculates the sum of each row in a numeric matrix, like the count matrix we have here, and it returns a vector. There's also a function called `which.max()` that determines the index of the maximum value in a vector.

0. Find the gene with the highest expression across all samples -- remember, each row is a gene.
0. Extract the expression data for this gene for all samples.
0. In which sample does it have the highest expression?
0. What is the function of the gene? Can you suggest why this is the top expressed gene?

```{r, echo=FALSE, include=FALSE}
        topGene <- which.max(rowSums(countdata))
        topGene
        countdata[topGene, ]
        # this is a pseudogene - maybe an artefact of only aligning reads to a single chromosome?
```


## DESeq2 analysis

DESeq2 is an R package for analyzing count-based NGS data like RNA-seq. It is available from [Bioconductor](http://www.bioconductor.org/). Bioconductor is a project to provide tools for analysing high-throughput genomic data including RNA-seq, ChIP-seq and arrays. You can explore Bioconductor packages [here](http://www.bioconductor.org/packages/release/BiocViews.html#___Software).

Just like R packages from CRAN, you only need to install Bioconductor packages once, then load them every time you start a new R session.

```{r install_deseq2, eval=FALSE}
        source("http://bioconductor.org/biocLite.R")
        biocLite("DESeq2")
```

```{r load_deseq22}
        library("DESeq2")
        citation("DESeq2")
```

It requires the count data to be in matrix form, and an additional dataframe describing sample metadata. Notice that the **colnames of the countdata** match the **rownames of the metadata*.

```{r read_coldata}
        mycoldata <- read.csv("data/coldata.csv", row.names=1)
        mycoldata
```

DESeq works on a particular type of object called a DESeqDataSet. The DESeqDataSet is a single object that contains input values, intermediate calculations like how things are normalized, and all results of a differential expression analysis. You can construct a DESeqDataSet from a count matrix, a metadata file, and a formula indicating the design of the experiment.

```{r make_deseqdataset}
        dds <- DESeqDataSetFromMatrix(countData=countdata, colData=mycoldata, design=~condition)
        dds
```

Next, let's run the DESeq pipeline on the dataset, and reassign the resulting object back to the same variable. Before we start, `dds` is a bare-bones DESeqDataSet. The `DESeq()` function takes a DESeqDataSet and returns a DESeqDataSet, but with lots of other information filled in (normalization, results, etc). Here, we're running the DESeq pipeline on the `dds` object, and reassigning the whole thing back to `dds`, which will now be a DESeqDataSet populated with results.

```{r run_deseq}
        dds <- DESeq(dds)
```

Now, let's use the `results()` function to pull out the results from the `dds` object. Let's re-order by the adjusted p-value.

```{r}
        # Get differential expression results
        res <- results(dds)
        head(res)
        
        # Order by adjusted p-value
        res <- res[order(res$padj), ]
        head(res)
```

Combine DEseq results with the original counts data. Write significant results to a file.

```{r write_results}
        sig <- subset(res, padj<0.05)
        dir.create("results")
        write.csv(sig, file="results/sig.csv") # tab delim data
```


## Data Visualization

We can also do some exploratory plotting of the data.

The differential expression analysis above operates on the raw (normalized) count data. But for visualizing or clustering data as you would with a microarray experiment, you ned to work with transformed versions of the data. First, use a *regularlized log* transofmration while re-estimating the dispersion ignoring any information you have about the samples (`blind=TRUE`). Perform a principal components analysis and hierarchical clustering.

```{r}
        # Transform
        rld <- rlogTransformation(dds)
        
        # Principal components analysis
        plotPCA(rld, intgroup="condition")
        
        # Hierarchical clustering analysis
        ## let's get the actual values for the first few genes
        head(assay(rld))
        ## now transpose those
        t(head(assay(rld)))
        ## now get the sample distances from the transpose of the whole thing
        dist(t(assay(rld)))
        sampledist <- dist(t(assay(rld)))
        plot(hclust(sampledist))
```

Let's plot a heatmap.

```{r plot_heatmaps}
        # ?heatmap for help
        sampledist
        as.matrix(sampledist)
        sampledistmat <- as.matrix(sampledist)
        heatmap(sampledistmat)
```

That's a horribly ugly default. You can change the built-in heatmap function, but others are better.

```{r gplots_heatmap}
        # better heatmap with gplots
        library("gplots")
        heatmap.2(sampledistmat)
        heatmap.2(sampledistmat, key=FALSE, trace="none")
        colorpanel(10, "black", "white")
        heatmap.2(sampledistmat, col=colorpanel(64, "black", "white"), key=FALSE, trace="none")
        heatmap.2(sampledistmat, col=colorpanel(64, "steelblue", "white"), key=FALSE, trace="none")
        heatmap.2(sampledistmat, col=colorpanel(64, "red", "white", "blue"), key=FALSE, trace="none")
```

What about a histogram of the p-values?

```{r plot_pval_hist}
        # Examine plot of p-values
        hist(res$pvalue, breaks=50, col="grey")
```

Let's plot an MA-plot. This shows the fold change versus the overall expression values.

```{r MA_plot}
        with(res, plot(baseMean, log2FoldChange, pch=16, cex=.5, log="x"))
        with(subset(res, padj<.05), points(baseMean, log2FoldChange, col="red", pch=16))
        
        # optional: label the points with the calibrate package. see ?textxy for help
        library("calibrate")
        res$Gene <- rownames(res)
        head(res)
        with(subset(res, padj<.05), textxy(baseMean, log2FoldChange, labs=Gene, cex=1, col="red"))
```

Let's create a volcano plot.

```{r volcano_plot}
        par(pch=16)
        with(res, plot(log2FoldChange, -log10(pvalue), main="Volcano plot"))
        with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), col="red"))
        with(subset(res, abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), col="orange"))
        with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), col="green"))
        # Add legend
        legend("topleft", legend=c("FDR<0.05", "|LFC|>2", "both"), pch=16, col=c("red","orange","green"))
        # Label points
        with(subset(res, padj<.05 & abs(log2FoldChange)>2), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=1))
```



## Going further
* Read about multifactor designs in the [DESeq2 vignette](http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf) for cases where you have multiple variables of interest (e.g. irradiated vs controls in multiple tissue types).





Beginning section Edited from [Training-modules](https://github.com/hbctraining/Training-modules) 
