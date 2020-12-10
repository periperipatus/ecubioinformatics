# Read Mapping and Quantification

This tutorial was conceived by Chris Balakrishnan and further modified by Peri Bolton.

# Table of Contents

* [Transcript quantification using Kallisto](#transcript-quantification-using-kallisto)
	* [Indexing the transcriptome](#indexing-the-transcriptome)
	* [Quantify Reads](#quantify-reads)
* [Splice-aware mapping with STAR](#splice-aware-mapping-with-star)

# Transcript quantification using Kallisto

## Index the transcriptome

Today we are going to quantify our reads with [Kallisto](https://pachterlab.github.io/kallisto/about).
It is a pseudo-alignment program for rapid quantification of transcripts. It is paired with R package [Sleuth](https://pachterlab.github.io/sleuth/about) for analysis of differential expression, but its output can be used with anything.
```Salmon``` is another commonly used transcript quantification tool. These are both integrated into R package ```tximport``` for use in downstream differential expression analysis. 


To be able to quantify the reads, we need a transcriptome. This is a fasta file that contains all the transcript isoforms known for an organism and/or tissue.

We have a *de novo* transcriptome for the golden-collared manakin *Manacus vitellinus*, assembled using [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)

Make sure that you are in ```~/<YOUR DIR>/``` on the server before we start. 

```bash
cp ../Bioinformatics_Workshop/MAVI_Trinity.fasta ./
ls

```


Transcriptomes tend to be big files, and for each read, the algorithm needs to search thus, efficient searches are needed. So the fasta file needs to be converted to an easily searchable file, a transcriptome index. 
To do this for kallisto, we need to invoke the program within kallisto:

```bash
kallisto index
```
Specifying the program alone without any arguments and it will tell you what arguments are required. Here you'll see that it is required to specify -i: the name and location of the index file to be created. 
The other requirement is that you specify the fasta file (transcriptome to be used as the reference to build the index).

```bash
kallisto index -i MAVI_kallisto_index MAVI_Trinity.fasta 
```
This step takes approximately 10 minutes
now let's make a directory to put our results into

```bash

mkdir ../kallisto_results
```
<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>

## Quantify reads

Now that you've created the index, you can quantify transcript abundance.
Below is an example for a single file. You can run for each read file individually (using tab to autocomplete for correct file names), or you can build it into a loop. 


```bash
cd MAVI_reads/
kallisto quant -i ../MAVI_kallisto_index -o ../kallisto_results/7_MAVI_SH_JB1 -b 30 --rf-stranded  ./7_MAVI_SH_JB1_F_val_1.fq ./7_MAVI_SH_JB1_R_val_2.fq
```
**Question 1:**  -o specifies what the output will be called and where it will go. What do the other options specify? (e.g., -b). Run ```kallisto quant``` or see the [manual online](https://pachterlab.github.io/kallisto/manual). 


**Try it: ** Using the skills from the previous two modules, run all the Manakin reads through Kallisto as a loop. **Hint** Use ```sed``` to remove/replace bits of the filename for input into ```kallisto```. Don't spend too long on this - we'll discuss in group later. 


Running these commands should generate 6 results directories with 3 output files for samples 7 through 12 (abundances.h5,  abundances.tsv, run_info.json).

Now, move back into your main directory and make a tab-delimited text file file that describes each of the samples as described in the ```MAVI_sample_key.xlsx```

```bash
cd ..

pico
#pico is a text editor program within the commmand line. You could use another like vi.

sample	tissue
7_MAVI_SH_JB1	SH
10_MAVI_PEC_JB1	PEC
11
12...
#Make sure that each word is separated by a tab, unless you are going to a new line
```

Press ctrl X to exit, and save the files as ```MAVI_samples.txt```

Tomorrow we will analyse the results from the full dataset.

<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>

# Splice-aware mapping with STAR

If you are working on species with a genome, you can also do splice aware mapping.
Some commonly used splice-aware aligners are [STAR](https://github.com/alexdobin/STAR) and [HISAT2](http://daehwankimlab.github.io/hisat2/about/).
Here we are going to use STAR (Spliced Transcripts Alignment to a Reference) to map our *Manacus* reads to the reference genome (Yes, I lied above, we have a genome assembly). 


Just like ```kallisto```, the ```STAR``` aligner needs to make a genome index so that it can efficiently access the genome. 
However, STAR is a little slow at doing the index. So genome index has already been made for you here ```~/Bioinformatics_Workshop/mavi_genome```

Here's the code used to arrive at this point. 

```bash
cd ~/Bioinformatics_Workshop/
mkdir mavi_genome/
cd mavi_genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/715/985/GCF_001715985.3_ASM171598v3/GCF_001715985.3_ASM171598v3_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/715/985/GCF_001715985.3_ASM171598v3/GCF_001715985.3_ASM171598v3_genomic.gtf.gz

gunzip *.gz #uncompress files because STAR doesn't like compressed files. 


STAR --runThreadN 50 \
--runMode genomeGenerate \
--genomeDir mavi_index \
--genomeFastaFiles GCF_001715985.3_ASM171598v3_genomic.fna \
--sjdbGTFfile GCF_001715985.3_ASM171598v3_genomic.gtf \
--sjdbOverhang 99  \
--genomeSAindexNbases 13 

#this took about 20 min
```

When running 

**Question 2:** If you wanted to download the transcripts for a transcriptome based quantification (like with Kallisto) which file would you download from the [GenBank RefSeq Assembly](https://www.ncbi.nlm.nih.gov/assembly/GCF_001715985.3)?

**Question 3:** What does the ```--sjdbOverhang``` argument do?

**Question 4:** What do the \ do in the code?

Now, you can run one of the trimmed fastqs from before.

```bash
mkdir ~/<YOURDIR>/alignments/

STAR --genomeDir ~/Bioinformatics_Workshop/mavi_genome/mavi_index/ \
--runThreadN 6 \
--readFilesIn 7_MAVI_SH_JB1_F_val_1.fq 7_MAVI_SH_JB1_R_val_2.fq \
--outFileNamePrefix ../alignments/7_MAVI_SH_JB1 \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes Standard \
--quantMode GeneCounts &
```

**Question 5:** What do the different arguments do in this? refer to the STAR manual


Now, let's move into our ```alignments/``` folder and look at what we've produced. There should be a number of files there including those ending in ```Aligned.sortedByCoord.out.bam```, ```Log.final.out``` and ```PerGene.out.tab```

First we can look at the bam file. This is our reads aligned. Because they are compressed in BGZF format you need to use a special program to make it human readable. 


```bash
samtools view 7_MAVI_SH_JB1Aligned.sortedByCoord.out.bam | head
```
This command shows the first few lines of the alignment section of the file. if you want to view the header you need to use the ```-h``` flag.

**Question 6:** What do the first 5 columns of the alignment tell you?

Now, let's open the ```Log.final.out``` file. 

**Question 7:** What is the percent of unmapped, multi-mapped and uniquely mapped reads? What could affect these numbers?
**Question 8:** What is the Number of splices? What does this mean?

Now let's have a look at the contents of  ```ReadsPerGene.out.tab```

```bash
less -S 7_MAVI_SH_JB1ReadsPerGene.out.tab
```

hit ```q``` to exit.

**Question 9:** What does this file tell us? What do the three columns mean?

Don't worry that there are a lot of zeros - don't forget you're using a truncated fastq file. The differential expression example tomorrow we will use a full dataset. 

**Question 10:** Other aligners, such as ```HISAT2``` do not provide any counts for reads against the gene models in the reference. How would you go about getting these counts?

<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>

## Alignment QC

Now, let's use `qualimap` to have a look at our alignment. 

Let's make a new directory for our quality reports and then run the program. 

```bash
mkdir qc

qualimap rnaseq \
-outdir qc/7_MAVI_SH_JB1_F_quali \
-a proportional \
-bam 7_MAVI_SH_JB1Aligned.sortedByCoord.out.bam \
-p strand-specific-reverse \
-gtf ~/Bioinformatics_Workshop/mavi_genome/GCF_001715985.3_ASM171598v3_genomic.gtf \
--java-mem-size=8G
```
Download the all the results to your desktop to view. This time you need to download the whole folder.

```bash
scp -r -P 1200 ngsclass@<IP.ADRESS>:~/<YOURDIR>/alignments/qc/7_MAVI_SH_JB1_F_quali/ .
```

**Question 11:** what does the `-r` flag do?

**Question 12:** Think about the quality of the mapping. Do we have good mapping rates? Why/why not? Do we have good representation in exonic regions?


# Your Assignment

Try your new knowledge on the dataset provided to you. 

<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>