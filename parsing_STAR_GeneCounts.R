library(plyr)
samples <- read.table(file.path("MAVI_samples.txt"), header = TRUE, stringsAsFactors=FALSE)
files<- list.files(path="STAR_results", pattern="*\\.tab")
files<- file.path("STAR_results",paste0(samples$sample, "ReadsPerGene.out.tab"))
data_list = lapply(files, read.table, col.names=c("gene_name", "unstr_counts", "sense", "antisense"), sep="\t")


samplenames<- sub("ReadsPerGene.out.tab", "", files)
samplenames<- sub("STAR_results/", "", samplenames)

#exclude the top part of the data
data_list<- lapply(data_list, function(x)x[5:length(x[,1]),])

#include only the gene names and the antisense counts. we want antisense because we used stranded RNA library prep and this column is equivalent to
data_list<- lapply(data_list, function(x)x[,c(1,4)])

#prepare data for inclusion in a matrix gene by sample.
data_list2<- list()
for(i in 1:length(data_list)){
  sub<- data_list[[i]]
  fname<- samplenames[i]
  sub<- rename(sub, c("antisense"=fname))
  data_list2[[i]]<- sub
}
data<- do.call("cbind", data_list2)

#remove duplicate gene_name columns
data<- data[, -seq(3, length(colnames(data)), 2)] 


#now we have a data frame that we can use in DEseq! Let's write it to file.
write.csv(data, "MAVI_STAR_results_compiled.csv", row.names=TRUE) 