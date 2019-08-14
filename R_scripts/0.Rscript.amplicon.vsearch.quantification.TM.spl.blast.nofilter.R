#Rscript
#Quantification of number of reads and abundance per species per sample for TIERMIX SAMPLES only

setwd("/home/mailys/Sequencing.results/Amplicon.analysis/3.Vsearch/3.7.Assignement")

dataa=read.csv("../3.6.Clustering/bfr.all.otutab.txt", stringsAsFactors=T, sep="\t", h=T)
files=list.files(".", pattern="TM.*.bfr.otus.nofilter.blastn.id97.qc200.out") #.out files with assignment (align folder)

name=vector()
target=vector()
count=vector()

for (i in files) {
assg=read.csv(i, h=F, sep="\t")
fname=strsplit(i, "[.]")[[1]]
nname=fname[2]
fname2=paste("TierMix",nname,"bfr", sep="")
message(paste("Processing", fname2, sep=" "))
nspl=which(colnames(dataa)==fname2)
subdataa=dataa[,c(1,nspl)]
subdataa=subset(subdataa, subdataa[,2]!=0)
message(paste(nrow(subdataa),"OTUs"))
subassg=subset(assg, assg$V1 %in% subdataa$X.OTU.ID)
	for (j in unique(subassg$V2)) {
	subtarget=subset(subassg, subassg$V2==j)
	message(paste(fname2, j, nrow(subtarget),"OTUs"))
	subtab=subset(subdataa, subdataa$X.OTU.ID %in% subtarget$V1)
	quanti=sum(subtab[,2])
	name=append(name, fname2)
	target=append(target,j)
	count=append(count, quanti)
	}
	message("")
}

ttable=as.data.frame(cbind(name, target, count))
ttable$count=as.numeric(as.character(ttable$count))

ttable$target=gsub("_.*_Assembly", "", ttable$target)

ttable2=cast(ttable, target~name, value="count")
ttable2[is.na(ttable2)]<-0
write.csv(ttable2, file="MC2.all.quantification.matrix.blastn.nofilter.csv", row.names=F)
