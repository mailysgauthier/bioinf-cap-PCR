#Rscript
#For the capture vsearch pipeline
#Count the abundance per species and put the results in a matrix 
args=commandArgs(trailingOnly = TRUE)

print("Output file = All.samples.affil.count.blastn.COI.csv")

if (length(args)!=2) {
  stop("Two argument must be provided : directory of input files and output files", call.=FALSE)
}

dirin=args[1] 
dirout=args[2]

library(reshape) #cast
setwd(dirin)

sample_name=vector()
sp_name=vector()
sp_numb=vector()

list_name=list.files(pattern=".aln.COI.blastn.uniq.id97.qc250.clean.out")

for (i in list_name) {
	dataa=read.csv(i, sep="\t", h=F)
	temp=paste(".aln.COI.blastn.uniq.id97.qc250.clean.out", sep="")
	name <- gsub(temp,"",i[grep(temp,i)])
	uniq_sp=unique(dataa[,2])
		for (j in uniq_sp) {
		count=length(which(dataa[,2]==j))
		sample_name=append(sample_name,name)
		sp_name=append(sp_name,j)
		sp_numb=append(sp_numb, count)}
	}
dataaf=as.data.frame(cbind(sample_name,sp_name,sp_numb))
dataaf2=cast(dataaf, sample_name~sp_name, value="sp_numb")

setwd(dirout)
write.csv(dataaf2, file="All.samples.affil.count.blastn.COI.csv", row.names=F)














