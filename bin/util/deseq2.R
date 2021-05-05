#!/usr/bin/env Rscript
library(DESeq2)
library(ggplot2)
### Load DESeq2Constructor
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
other.name <- file.path(script.basename, "/../DESeq2Constructor.R")
source(other.name)
# source("~/git/IRFinder/bin/DESeq2Constructor.R")
### Read args
# setwd("~/test/IRFinder2/Diff/sing/")
# args=c("./groups.tsv", "0.05", "0" , "0" ,"0")

args <- commandArgs(trailingOnly = T)
groups=read.table(args[1], stringsAsFactors = F, header = T)
out_folder=dirname(args[1])
IRratio_thr=as.numeric(args[2])
warning_filter=args[3]
cooks_cutoff=args[4]=="1"
if (cooks_cutoff ){
  print("cooks_cutoff enabled")
} else {
  print("cooks_cutoff disabled")
}

independentFiltering= args[5]=="1"
if (independentFiltering ){
  print("independentFiltering enabled")
} else {
  print("independentFiltering disabled")
}

paths = as.vector(groups$Files)
experiment = groups[,c("SampleName", "Condition")]

experiment$Condition=factor(experiment$Condition) 
rownames(experiment)=NULL

metaList=DESeqDataSetFromIRFinder(filePaths=paths, designMatrix=experiment, designFormula=~1, irratio_thr=IRratio_thr, warning_filter=warning_filter )

dds = metaList$DESeq2Object
design(dds) = ~Condition + Condition:IRFinder     
conditions=levels(experiment$Condition)
dds = DESeq(dds)
resultsNames(dds)

for ( i in 1:(length(conditions)-1) ){
  for (j in (i+1):length(conditions)){
    contrast_name=paste0(conditions[i], "_", conditions[j])
    res = results(dds, contrast=list(paste0("Condition", conditions[i] ,".IRFinderIR"),paste0("Condition", conditions[j] ,".IRFinderIR")), cooksCutoff=cooks_cutoff, independentFiltering=independentFiltering)  
    res$padj[is.na(res$padj)]=1
    if ( sum(res$padj < 0.05 ) > 0 ){
      pdf(paste0(out_folder, "/", contrast_name, "_plot.pdf"))
      nn_counts = counts(dds, normalized=F)
      for ( name in rownames(res)[res$padj < 0.05]){
        dat = data.frame( name = experiment$SampleName, grp= experiment$Condition,
                          intron_depth = nn_counts[name, paste0("intronDepth.", experiment$SampleName)] ,   
                          max_splice= nn_counts[name, paste0("maxSplice.", experiment$SampleName)])
        dat$IRratio = dat$intron_depth / ( dat$intron_depth+dat$max_splice)
        print(ggplot(dat)+geom_boxplot(aes(x=grp, fill=grp, y=IRratio )) + ggtitle(paste0(name, "\n", res[name, "padj"])))
      }
      dev.off()
    }
    write.table(res, file = paste0(out_folder, "/", contrast_name, "_DESeq2.tsv") ,sep="\t", quote = F)
  }
}

quit(save = "no")



