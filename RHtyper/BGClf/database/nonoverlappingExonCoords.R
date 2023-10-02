#!/usr/bin/env Rscript
suppressMessages(library("dplyr"))

args=commandArgs(trailingOnly=TRUE)


file=args[1]


df <- read.table(file, header=F, sep="\t")
colnames(df) <- c("chr","start","end","id","frame","strand")
df <- df[grep("^NM_", df$id),]


df.collapse <- unique(df[,c('chr', 'start', 'end')])
attach(df.collapse)
df.collapse <- df.collapse[order(chr, start, end),]

detach(df.collapse)
write.table(df.collapse, file=gzfile(paste(file, 'nodup', 'gz', sep='.')), sep="\t", row.names=F, quote=F)


