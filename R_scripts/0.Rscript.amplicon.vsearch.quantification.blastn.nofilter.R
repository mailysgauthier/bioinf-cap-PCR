#Rscript
#Quantification of number of reads and abundance per species per sample

args=commandArgs(trailingOnly = TRUE)

if (length(args)==0) {
  stop("Error : please give at least two arguments 
  At least one argument must be provided : 
  args[1] = .txt file with abundance per OTU ($clust folder)
  args[2] = .paf files with assignment (align folder)
  Optional arguments
  args[3] = directory for args[2]/paf file and output table", call.=FALSE)
  } else if (length(args)==1) {
  stop("Error : please give at least two arguments 
  At least one argument must be provided : 
  args[1] = .txt file with abundance per OTU ($clust folder)
  args[2] = .paf files with assignment (align folder)
  Optional arguments
  args[3] = directory for args[2]/out file and output table", call.=FALSE)
  } else if (length(args)==2){
  args[3]=getwd()
  message("Current directory will be used")} 

setwd(args[3])

dataa=read.csv(args[1], stringsAsFactors=T, sep=" ", h=T) #.txt file with abundance per OTU ($clust folder)
assg=read.csv(args[2], h=F, sep="\t") #.out files with assignment (align folder)

uniq.spl=names(dataa)[2:ncol(dataa)]

name=vector()
target=vector()
count=vector()

for (i in uniq.spl) {
message(paste("Processing", i, sep=" "))
n=which(colnames(dataa)==i)
subdataa=dataa[,c(1,n)]
subdataa=subset(subdataa, subdataa[,2]!=0)
message(paste(nrow(subdataa),"OTUs"))
subassg=subset(assg, assg$V1 %in% subdataa$X.OTU.ID)
	for (j in unique(subassg$V2)) {
	subtarget=subset(subassg, subassg$V2==j)
	message(paste(i, j, nrow(subtarget),"OTUs"))
	subtab=subset(subdataa, subdataa$X.OTU.ID %in% subtarget$V1)
	quanti=sum(subtab[,2])
	name=append(name, i)
	target=append(target,j)
	count=append(count, quanti)
	}
	message("")
}

ttable=as.data.frame(cbind(name, target, count))
ttable$count=as.numeric(as.character(ttable$count))
fname=strsplit(args[2], "[.]")[[1]]
fname2=paste(fname[1:2], collapse=".")
write.table(ttable, file=paste(fname2, "quantification.blastn.nofilter.txt", sep="."))

library(reshape) #cast
ttable2=cast(ttable, target~name, value="count")
ttable2[is.na(ttable2)]<-0
write.table(ttable2, file=paste(fname2,"quantification.matrix.blastn.nofilter.txt", sep="."), row.names=F)



