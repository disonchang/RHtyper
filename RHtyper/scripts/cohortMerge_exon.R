#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("dplyr"))

source("/home/tchang1/software/bin/R/paper_theme.R")

parser <- ArgumentParser(description='Merge bloodtype result for a cohort and generate summary plots')
parser$add_argument("-d", '--input_directory', type="character", help='folder with the CNV results')
parser$add_argument("-e", "--extension", type="character", help='file extension')


args <- parser$parse_args()

###
mergeFile <- function(path, extension, r=FALSE){
	out=data.frame()
	files=list.files(path=path, pattern=extension, recursive=r, full.names=TRUE)
        #print(files)

	for(f in files){
                fn=basename(f)
                sample=gsub("\\..*","",fn)
		temp=read.table(f, sep="\t", header=T, stringsAsFactors=FALSE)
                temp$sample <- sample
		out <- rbind(out, temp)
		rm(temp)
	}
	return(out)

}
###
indata <- mergeFile(args$input_directory, args$extension)

agg_indata <- indata %>% group_by(exon) %>%
	summarize(RHDmed=median(RHD_cov),
                  RHCEmed=median(RHCE_cov)
                 )

print(agg_indata)
indata$exon <- as.factor(indata$exon)
p1 <- ggplot(indata, aes(x=exon, y=RHD_cov))+geom_boxplot(outlier.shape = NA)+paper_theme+geom_jitter(shape=16, position=position_jitter(0.2), size=0.2)+
	xlab("Exon")+ylab("RHD coverage")
p2 <- ggplot(indata, aes(x=exon, y=RHCE_cov))+geom_boxplot(outlier.shape = NA)+paper_theme+geom_jitter(shape=16, position=position_jitter(0.2), size=0.2)+
	xlab("Exon")+ylab("RHCE coverage")

pdf('RH_exon_cov.pdf', width=7, height=7)
multiplot(p1,p2,cols=1)
dev.off()


