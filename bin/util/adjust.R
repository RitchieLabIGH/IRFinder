#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = T)
dat=read.table(args[1], stringsAsFactors = F, header = T)
cols=colnames(dat)
for (cn in cols[grepl(pattern ="p.val", x = cols)]  ){
  dat[,paste0(cn,"_BH_adjusted")]=p.adjust(dat[,cn], method = "BH")
}
write.table(x=dat, file = paste0(args[1], "_adjusted.tsv"), row.names = F, col.names = T, quote = F, sep="\t")

