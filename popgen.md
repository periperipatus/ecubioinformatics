# Population Genomics

# Table of contents

* [Mapping and SNP calling](#mapping-and-SNP-calling)
* [Association Testing](#association-testing)

# Mapping and SNP calling

This tutorial is (heavily) based on [Samtools primer](https://github.com/ecerami/samtools_primer) tutorial by Ethan Cerami. 


Within your directory on the server, make a new directory for this tutorial and then enter it.

First, we are going to get the *E. coli* genome, then subset it to make it smaller for the purposes of a speedy tutorial.

```bash
wget https://raw.githubusercontent.com/ecerami/samtools_primer/master/tutorial/genomes/NC_008253.fna

~/programs/wgsim/wgsim -N10000 -S1 NC_008253.1K.fna NC_008253_1K.read1.fq NC_008253_1K.read2.fq
```

Just like in the main RNAseq tutorials, before you can do an alignment you need to index the genome. Here we are using the software [bwa](https://github.com/lh3/bwa). Which stands for Burrows-Wheeler Aligner. Recall the Burrows-Wheeler index methodology from the lecture on mapping.

## Alignment

```bash
bwa index NC_008253.1K.fna
```

Now we are going to map these simulated reads back to the genome. 

```bash
bwa mem NC_008253.1K.fna NC_008253_1K.read1.fq NC_008253_1K.read2.fq > NC_008253_1K.sam

less -S NC_008253_1K.sam
```

Recall our lecture that went through the SAM file specifications. 
Now we want to convert that SAM file into a BAM file which is a binary format that saves a lot of space. 


```bash
samtools view -b -S  -o NC_008253_1K.bam NC_008253_1K.sam
```

**Question 1:** What do the -b and -S arguments mean?

Now, let's look at the BAM... 

```bash
less -S NC_008253_1K.bam
```

Errr? 
bamfiles are binary and thus not human readable.
what if you want to look at a bam file? you can convert it back to sam format, or you can do this:


```bash
samtools view NC_008253_1K.bam | less -S
```

Samtools includes a bunch of tools that can be used to interrogate and filter SAM and BAM files
For example the `-f` argument enables you to interrogate what information is in the ['flag'](https://broadinstitute.github.io/picard/explain-flags.html).
For example the flag 0x4 (or 4) indicates unmapped reads. 

```bash
samtools view -c -f 4 NC_008253_1K.bam
```

**Question 2:** How many reads are unmapped? What does the `-c` argument do?

```bash
samtools view -c -q 42 NC_008253_1K.bam
```

**Question 3:** How many reads are aligned with a quality of 42 (max quality score from bwa)?

**Question 4:** Consider how you might use these functions to filter out low quality alignments.

<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>

## Indexing and Sorting the BAM

For some applications, such as SNP calling, sam files have to be sorted by genomic position

```bash
samtools sort -o NC_008253_1K.bam.sorted NC_008253_1K.bam
```

Like genomes, BAM and SAM files can be indexed also, and is required by many downstream applications. 
This command will create another file with the same name with a `.bai` appended to the end. It is important that these file names have the same prefix (everything before .bai) otherwise programs that require an indexed BAM wont know where to look.

```bash
samtools index NC_008253_1K.bam.sorted
```

Sorted and Indexed files are required for viewing the alignment in IGV. The sorted and indexed alignment can also be viewed in the CLI using the following command:


```bash
samtools tview NC_008253_1K.bam.sorted NC_008253.1000.fna
```


<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>

## Variant calling

There are a number of tools for variant calling `samtools mpileup` is one of them. 
Another commonly used software package is the [Genome analysis Took Kit (GATK)](https://gatk.broadinstitute.org/hc/en-us)

```bash
samtools mpileup -g -f NC_008253.fna NC_008253_1K.bam.sorted > sim_variants.bcf
bcftools call -v -c sim_variants.bcf > sim_variants.vcf
```

Look at the VCF file. I'm not going to tell you how to do that. 

For a full specification on the VCF file format please go [here](https://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/).
This has some core information stored in the header (everything preceded with `##`) about what the reference genome was, and what software was used to call. 
The final line of the header will be the column names for the data included below (indicated with a single `#`).
Now the file is essentially 1 row per variant, with information on the scaffold/chromosome, precise position and what the ref and alternate bases are.
When you have multiple individuals in the file, they will be in columns. 
To merge multiple VCFs from the same experiment you can use `bcftools merge` or `PicardTools MergeVcfs`.

**Question 5:** Now that we know header information is specified with a `#`, find out how many sites have been called as variant sites.


<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>

# Association Testing

This part of the tutorial is derived from the documentation on [association testing with GPAT++ and vcflib](https://github.com/zeeev/vcflib/wiki/Association-testing-with-GPAT)

Here we will estimate pFst from variant call data, which incorporates information from the Genotype Likelihood information in the VCF file. 

Within the `vcflib` directories there are some sample data in the "samples" directory.
We will do a quick analysis of `scaffold612.vcf`  to demonstrate what you can do with variants based on the frequency in different populations.

The VCF is here: `~/programs/vcflib/samples/scaffold612.vcf`

Calculate the pFST across this scaffold. The numbers refer to the index position of the individuals in the vcf file in the two populations being compared.

```bash
pFst --target 1,20,25,29,30,38,43,46 --background 2,3,4,5,6,7,21,22,23,24,26,26,28,31,32,33,34,35,36,37,39,40,41,42,44,45  --deltaaf 0.0 --file ~/programs/vcflib/samples/scaffold612.vcf --counts --type PL > 612.counts
```

**Question 6:** What does the `--deltaaaf` argument do? In what circumstances would you want to set this to >0?

Then you can use a premade plotting script. You can invoke it in the unix CLI like so:

```bash
R --vanilla < ~/programs/vcflib/bin/plotPfst.R --args 612.counts 
```

There are other cool plotting and analysis scripts in GPAT/vcflib feel free to mess around with them. 
