# Differential Expression of Transcripts from Kallisto

Yesterday you quantified transcripts from our experiment on different muscle types to our *de novo* assembly of *Manacus vitellinus*.
Today we are going to explore differential expression of these *transcripts* using ```DESeq2``` a standard method of evaluation differential expression.


## Move files to your desktop


Now, download ```~/<YOUR DIR>/MAVI_samples.txt``` from your directrory,  and the whole folder ```~/Bolton/kallisto_results/``` to your desktop

To do this, open a new terminal window and navigate to to the folder you made for the tutorial.

copy the file and folder across. The example below is for the folder. 
```bash
scp -r -P 1200 ngsclass@IP.ADRESS:~/Bolton/kallisto_results/ .
```

**Question 1:** What is different about this command vs. what we have explored earlier for ```scp```


## Differential Expression

Open R Studio on your desktop, and load the following packages into your environment:

```r
library(rhdf5)
library(tximport)
library(DESeq2)
```

Then, make sure you're in your correct working directory. I made my tutorial folder on my desktop so it's stored here: 

```r
setwd("C:/Users/perif/Desktop/BioinformaticsWorkshop")
```

### Tximport

Read in the sample metadata
```r
samples <- read.table(file.path("MAVI_samples.txt"), header = TRUE, stringsAsFactors=FALSE)
samples<- samples[order(samples$sample),] ## make sure the sample names are in the order a computer would read them...
```


Now, let's make a named vector that contains the path to all the files we want to read in.

```r
files<- file.path("kallisto_results", samples$sample, "abundance.h5")
names(files) <- samples$sample
```

Then we can use ```tximport``` to load in the data. 

```r
txi <- tximport(files, type = "kallisto", txOut = TRUE)
```

**Question 2:** What does the ```type=``` and ```txOut=``` arguments specify here?

Let's get a snapshot of our data

```r
names(txi)
head(txi$counts)
```
If you have something that looks like this ####################################### you are good to go.

### DESeq2

The sample metadata table needs to have row names that match the colum names for the ```txi``` object.

```r
rownames(samples)<- samples$sample
all.equal(rownames(samples), colnames(txi$counts)) #check that they're the same
```

If the last line returned ```TRUE``` you are good to go. 

```r
dds <- DESeqDataSetFromTximport(txi, samples, ~tissue)
dds<- DESeq(dds) #runs differential expression
```

Ok, now we have to extract the differential expression results from the complicated ```DESeq2``` object called ```dds``` using the wrapper function ```results()```


```r
res <- results(dds, alpha=0.05)
res<- res[order(res$padj),] #orders results so the most significant are at the top.
res
```

**Question 3:** What does the ```alpha=``` argument specify? 
**Question 4:** Why do the results need a pvalue adjustmnet?
**Question 5:** What does the ```log2FoldChange``` column mean?

For a summary of the output we can use the ```summary()``` function on the results object.

```r
summary(res)
````

**Question 6:** How many transcripts are differentially expressed? 

Now, let's plot a volcano plot. 

```r
plotMA(res)
```

**Question 7:** What do the dots represent? What are the blue ones? 



# Differential Expression of Genes from *Pipra filicauda* tissues

