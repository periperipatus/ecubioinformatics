# Annotations & Gene Ontology

This tutorial was originally conceived by Chris Balakrishan and has been modified by Peri Bolton.

Here, we will use our Trinity assembled *Manacus vitellinus* transcriptome and annotate it based on the Zebra finch (*Taeniopygia guttata*) genome. This species has the best genome of the closely related birds. 

This is a similar process to what would happen if you did an RNA-seq experiment *de novo*. 


# Table of Contents

* Annotation with BLAST (#annotation-with-blast)
* Ensembl BioMart (#ensembl-biomart)
* Functional Enrichment (#functional-enrichment)

# Annotation with BLAST

Make sure you are in your directory.

We are going to download Ensembl Zebra finch cDNA file, that contains known and predicted transcripts from the genome assembly. Then we need to unzip it.

```bash
wget http://ftp.ensembl.org/pub/release-102/fasta/taeniopygia_guttata/cdna/Taeniopygia_guttata.bTaeGut1_v1.p.cdna.all.fa.gz
gunzip Taeniopygia_guttata.bTaeGut1_v1.p.cdna.all.fa.gz

```

**Question 1:** Which Ensembl release is this? When was it released?

NCBI has a way to run a BLAST search using your local computer. This is called [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/).

First, you need to tell it what species you want to search against. Make a BLAST database of the Zebra finch cDNA using ```makeblastdb```. 

```bash
makeblastdb -in Taeniopygia_guttata.bTaeGut1_v1.p.cdna.all.fa -out ZFdb -dbtype 'nucl' -parse_seqids -hash_index &
```

**Question 2:** What do the flags mean?

Now, we BLAST our transcriptome against the Zebra finch BLAST database.


```bash
blastn -db ZFdb -query MAVI_Trinity.fasta -num_threads 1 -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > MAVIvsZEFI.blast.txt
```

**Question 3:** What is the e-value?

**Question 4:** What do the flags ```-max_target_seqs``` and ```-outfmt 6``` do?

Now we have the output from our differential transcript expression analysis, and this file. We need to associate the BLAST results with the differential expression results.

Let's have a look at the BLAST results first. 

```bash

head MAVIvsZEFI.blast.txt

```

**Question 5:** What do each of the columns mean? Use the NCBI documentation. Which part of the header in Zebra finch cDNA file has it used in the results?



Now, how can we use these results to extract functional information from our gene expression results?
	1. We can create a key file that links our Trinity transcript name to a gene name for gene-level expression analysis (import using `tximport` and convert transcript counts to gene counts with `tx2gene` (see [documentation](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html))
	2. Directly annotate the results from differential transcript analysis then perform functional analysis.


For the purposes of this exercise we will go through option 2, using [Ensembl](https://uswest.ensembl.org/index.html) BioMart to convert transcript IDs to geneIDs for downstream GO analysis. This is in order to familiarise yourself with using BioMart and not repeat Yesterday's differential expression exercise. 


Now, let's download our Blast results file `MAVIvsZEFI.blast.txt` to our desktop. We are going to merge our differential expression results with the BLAST results using R.


Now, let's load all our data into R-Studio.

```r
library(dplyr)
dte<- read.csv("MAVI_results_DE_transcripts.csv")
dte<- dte %>% rename("transcript_name"="X")

blastres<- read.table("MAVIvsZEFI.blast.txt", header=FALSE, sep="\t")
blastres<- blastres  %>% rename("transcript_name"="V1","ZEFI_transcript"="V2")
blastres<- blastres[,-c(3:12)]
blastres$ZEFI_transcript<- gsub("\\.1","", blastres$ZEFI_transcript) #removing the version number
blastres<- blastres[!duplicated(blastres$transcript_name),]  #Keep only the first ('best') hit per query 
```

**Question 6:** Check on the number of row names in each. Why do you think there might be different number of rows in each? How does this relate to your BLAST parameters above?

Now, let's combine the results, based on the shared column `"transcript_name"`, and remove any p-values with `NA` values.

```r
dte_annotated<- merge(blastres, dte, by="transcript_name")
dte_annotated<- dte_annotated %>% filter(!is.na(padj))

```

**Question 7:** How many rows did we remove when we took out the `NA` values? Why were those p-values `NA`? Think back to the differential expression lesson yesterday.


Now, let's export these data and use a GUI for once!

```r
write.csv(dte_annotated, "MAVI_DTE_annotated.csv")
```

<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>

# Ensembl BioMart

Annotation with human Gene IDs.
Ensembl is handy in that you can associate stable GeneIDs across species. Here, for each 
zebra finch ID (and therefore manakin transcript), we will extract a human gene name and ID.

1. Open the file in excel (at this point, the file should be small enough)
2.  Go to ensemble.org in your web browser.
3.	At the top, click on “BioMart”
4.	Choose database (Ensembl Genes 102, Ensembl updates its gene annotations/databases periodically)
5.	Click the “dataset” pulldown and select “Zebra Finch” (Note: zebra finch ≠ zebra fish)
6.	Click on “Filters” Section and the “gene” subsection
7.	In the “input external references” pulldown select “Transcript ID”
8.	Copy/paste your zebra finch ID list (only one column) into the space below the pulldown.
9.	Click the "Attributes" section.
10.	Check "Gene Name" and "Description" in the "gene" subsection
11.	Click “Results” at the top, and download the file as a tsv or csv by clicking "Go"

Note, we can also use BioMart to find gene orthologues across species using the "Homologs" function.

## back to R-Studio...

Now, we need to realign these gene names with our original data and rank them in order of p-values. We need to rank them because we will use a tool that assesses GO enrichment based on pvalue rank-order.

Each step use the R data viewer to see what the transformations have done. 

```r
biomart<- read.table("mart_export.txt", sep="\t", header=TRUE)
biomart<- biomart %>% rename("ZEFI_name"="Transcript.stable.ID")
biomart<- biomart %>% select(ZEFI_name, Gene.stable.ID, Gene.name) #keeping only the columns of interest

dte_hgnc<- merge(biomart, dte_annotated, by="ZEFI_transcript") 
dte_hgnc<- dte_hgnc[dte_hgnc$Gene.name!="",] #remove gene names with no universal gene name
dte_hgnc<- dte_hgnc[order(dte_hgnc$padj, decreasing = FALSE),] #smallest p-values first

write.csv(dte_hgnc, "MAVI_DTE_hcnc_names.csv")
```

**Question 8:** Why do some genes with an ENSTGUG ID not have a Gene name?

**Question 9:** Now, have a look at the fasta headers in the Zebra finch transcriptome. How could we circumvent the BioMart process?

<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>

# Functional enrichment

We can use these gene names to test for functional enrichment. Actually we could do this with zebra finch gene IDs too, just using a different software. 
Two commonly used approaches are Gene Ontology (GO) and KEGG (Kyoto Encyclopedia of Genes and Genomes) Pathway Analysis. Furthermore, there are many different tools for testing functional enrichment. 
One popular R tool is `clusterProfiler`.

We are going to use web-based GOrilla for testing for GO enrichment based on rank order of p-values. 

You can figure out how to use them yourself. But basically, just paste the list of differentially
expressed genes where it is requested Using the "Gene.name" column. We will use the gene names to compare against the human database of GO terms.  Enrich amongst all GO categories

http://cbl-gorilla.cs.technion.ac.il
 
**Question 10:**	What GO categories are significantly enriched among differentially expressed genes according to Biological Process? Does this make sense based on our predictions?

**Question 11:**	What GO categories are significantly enriched among differentially expressed genes according to Function? Does this make sense based on our predictions?

** On your own:** Try this using the results from the gene-level differential expression analysis (the second exercise yesterday).

<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>

# Assignment

Try this on your own data.
