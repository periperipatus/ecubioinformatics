# Annotations and GO

Here, we will use our Trinity assembled *Manacus vitellinus* transcriptome and annotate it based on the Zebra finch (*Taeniopygia guttata*) genome. This species has the best genome of the closely related birds. 

This is a similar process to what would happen if you did an RNA-seq experiment *de novo*.

Please make sure you understand what each of the commands do as you follow along with the tutorial. 

Make sure you are in your directory.

We are going to download Ensembl Zebra finch cDNA file (Release 89), that contains known and predicted transcripts from the genome assembly. Then we need to unzip it.

```bash
wget ftp://ftp.ensembl.org/pub/release-89/fasta/taeniopygia_guttata/cdna/Taeniopygia_guttata.taeGut3.2.4.cdna.all.fa.gz
gunzip Taeniopygia_guttata.taeGut3.2.4.cdna.all.fa.gz

```

Make a BLAST database of the cDNA using ```makeblastdb```

```bash
makeblastdb -in Taeniopygia_guttata.taeGut3.2.4.cdna.all.fa -out ZFdb -dbtype 'nucl' -parse_seqids -hash_index &
```

**Question 1:** What do the flags mean?

Now, we BLAST our transcriptome against the Zebra finch BLAST database.


```bash
blastn -db ZFdb -query MAVI_Trinity.fasta -num_threads 1 -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > MAVIvsZEFI.blast.txt
```

**Question 2:** What is the e-value?
**Question 3:** What do the flags ```-max_target_seqs``` and ```-outfmt 6``` do?

Now we have the output from our differential transcript expression analysis, and this file. We need to associate the BLAST results with the differential expression results.

Let's have a look at the BLAST results first. 

```bash

head MAVIvsZEFI.blast.txt

```

**Question 4:** What do each of the columns mean? Use the NCBI documentation. 



