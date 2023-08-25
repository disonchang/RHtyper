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
		temp=read.table(f, sep="\t", header=T)
		out <- rbind(out, temp)
		rm(temp)
	}
	return(out)

}
###
indata <- mergeFile(args$input_directory, args$extension)
write.table(indata, file=paste(args$gene, 'merged.data.txt', sep='.'), quote=F, row.names=F, sep='\t', na="")

#print(head(indata))

indata.sum <- indata %>% group_by(ID) %>%
	summarize(allele_n=n())

problem_sample=indata.sum[which(indata.sum$allele_n != 2),]

if (nrow(problem_sample) > 0){
	message("Problematic samples")
        print(problem_sample)
	#print(indata[which(indata$ID %in% problem_sample$ID),])
	
	remove_alias <- c('RHCE*ce48C, 733G','RHCE*Ce941C')
	indata <- subset(indata, ! (indata$ID %in% problem_sample$ID & indata$Alias %in% remove_alias))

        index2rm <- indata[which(indata$ID=='SJSCD045677_G1' & indata$Alias=='RHCE*CeNR'),'Index']
        indata <- subset(indata, ! (indata$ID=='SJSCD045677_G1' & indata$Index==index2rm))

        index2rm=max(indata[which(indata$ID=='SJSCD045383_G1'),'Index'])
	indata <- subset(indata, ! (indata$ID=='SJSCD045383_G1' & indata$Index==index2rm))	

	
	indata.sum <- indata %>% group_by(ID) %>% summarize(allele_n=n())
	problem_sample=indata.sum[which(indata.sum$allele_n != 2),]
	if (nrow(problem_sample) > 0){
		print("problem remains")
                print(problem_sample)
                print(indata[which(indata$ID %in% problem_sample$ID),])
	} else {
		print(nrow(indata.sum )) 
	}
	
}


write.table(indata, file=paste(args$gene, 'merged.cleaned.data.txt', sep='.'), quote=F, row.names=F, sep='\t', na="")

indata$Alias <- as.character(indata$Alias)
indata$Bloodtype <- as.character(indata$Bloodtype)
indata$Alias <- ifelse(is.na(indata$Alias), indata$Bloodtype, indata$Alias)
indata$Alias <- ifelse(indata$Alias=='', indata$Bloodtype, indata$Alias)


#order <- allele_freq[order(allele_freq$Allele_n),'Alias']
#allele_freq$Alias <- factor(as.character(allele_freq$Alias), levels=order)

allele_freq <- data.frame(indata %>% group_by(Bloodtype) %>% 
	summarize(Allele_n=n(), Alias=paste(unique(Alias),collapse='/'))
	)
allele_freq$Alias <- gsub(" ","",allele_freq$Alias)
allele_freq$Allele_p <- round(allele_freq$Allele_n/sum(allele_freq$Allele_n),4) * 100

order <- allele_freq[order(allele_freq$Allele_n),'Alias']
allele_freq$Alias <- factor(as.character(allele_freq$Alias), levels=order)
order2 <- allele_freq[order(allele_freq$Allele_n),'Bloodtype']
allele_freq$Bloodtype <- factor(as.character(allele_freq$Bloodtype), levels=order2)

print(allele_freq)
print(sum(allele_freq$Allele_n))

write.table(allele_freq, file=paste(args$gene, 'allele_freq','txt', sep='.'), row.names=F, quote=F, sep="\t")


width=nrow(allele_freq)/3

#colors=getPalette(5, 'Set2')
colors=c('#F61161','#0D78BD')
names(colors) <- c('RHD','RHCE')
## #80B1D3
## #B80D48


bloodtype_n=length(allele_freq$Alias)
bloodtype_label=allele_freq$Bloodtype

p <- ggplot(allele_freq, aes(x=Alias, y=Allele_n))+
	geom_bar(colour="black", fill=colors[args$gene], stat = "identity")+
	geom_text(aes(x=Alias, y=Allele_n, label=paste(Allele_n, "(",Allele_p,"%)", sep="")), hjust=-0.05)+
	#annotate("text", x = 1:length(allele_freq$Bloodtype), y = -50, label = allele_freq$Bloodtype, size = 1)+
        paper_theme_flip+
	#coord_cartesian(ylim = c(0, max(allele_freq$Allele_n)+250), expand = FALSE, clip = "off") +
        coord_flip()+
     	ylim(0, max(allele_freq$Allele_n)+300) +
	xlab("")+
	ylab("Allele frequency")
p2 <- ggplot(allele_freq, aes(x=Bloodtype, y=Allele_n))+
        geom_bar(colour="black", fill=colors[args$gene], stat = "identity")+
        geom_text(aes(x=Bloodtype, y=Allele_n, label=paste(Allele_n, "(",Allele_p,"%)", sep="")), hjust=-0.05)+
        #annotate("text", x = 1:length(allele_freq$Bloodtype), y = -50, label = allele_freq$Bloodtype, size = 1)+
        paper_theme_flip+
        #coord_cartesian(ylim = c(0, max(allele_freq$Allele_n)+250), expand = FALSE, clip = "off") +
        coord_flip()+
        ylim(0, max(allele_freq$Allele_n)+300) +
        xlab("")+
        ylab("Allele frequency")

ggsave(p, file=paste(args$gene, 'allele_freq','pdf', sep='.'), width=7)
ggsave(p2, file=paste(args$gene, 'allele_freq2','pdf', sep='.'), width=7)



