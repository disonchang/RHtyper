#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("dplyr"))

source("/home/tchang1/software/bin/R/paper_theme.R")

parser <- ArgumentParser(description='Merge bloodtype result for a cohort anf generate summary plots')
parser$add_argument("-d", '--input_directory', type="character", help='folder with the prediction results')
parser$add_argument("-e", "--extension", type="character", help='file extension')
parser$add_argument("-g", "--gene", type="character", help='gene')


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
		temp <- subset(temp, ! ALT %in% c('-'))
		out <- rbind(out, temp)
		rm(temp)
	}
	return(out)

}
###
indata <- mergeFile(args$input_directory, args$extension)
indata$REF <- gsub("TRUE", "T", indata$REF)
indata$ALT <- gsub("TRUE", "T", indata$ALT)
indata$REFaa <- gsub("TRUE", "T", indata$REFaa)
indata$ALTaa <- gsub("TRUE", "T", indata$ALTaa)
indata$REFaa <- gsub("FALSE", "F", indata$REFaa)
indata$ALTaa <- gsub("FALSE", "F", indata$ALTaa)

indata$class <- ifelse(indata$ALTaa==indata$REFaa, "silent", ifelse(grepl("\\*", indata$ALTaa), 'nonsense', 'missense'))
indata$acc <- gsub("m$","",indata$acc)
indata$aapos <- as.integer(indata$aapos)
indata$aachange <- paste("p.", indata$REFaa, indata$aapos, indata$ALTaa, sep="")
#patient	Chromosome	start	class	refseq	gene	aachange	REF	ALT	VAF	Disease
#SJBALL042232_R2	16	28944770	frameshift	NM_001178098	CD19	p.Y259fs	T	TGT	0.52	Relapse2


write.table(indata, file=paste(args$gene, 'merged.var.data.txt', sep='.'), quote=F, row.names=F, sep='\t', na="")

subdata <- indata[,c('sample','chr','gpos','class','acc', 'gene', 'aachange', 'REF','ALT')]
colnames(subdata) <- c('patient','Chromosome','start','class','refseq','gene','aachange','REF','ALT')

write.table(subdata, file=paste(args$gene, 'merged.var.proteinpaint.txt', sep='.'), quote=F, row.names=F, sep='\t', na="")


####
unique_var <- subdata %>% group_by(gene, Chromosome,start,REF,ALT,class,aachange) %>%
	summarize(sample_n=n())
write.table(unique_var, file=paste(args$gene, 'merged.var.unique.txt', sep='.'), quote=F, row.names=F, sep='\t', na="")



