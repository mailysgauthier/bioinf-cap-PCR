#Rscript
#For the amplicon vsearch pipeline
#Make one table for TM samples and MC samples (BFR amplicon)to select OTUs present in these samples
args=commandArgs(trailingOnly = TRUE)

if (length(args)==0) {
  stop("Error : please give one or two arguments 
  At least one argument must be provided : 
  args[1] = input table / OTUs table.txt
  Optional arguments
  args[2] = directory for input and output table", call.=FALSE)
  } else if (length(args)==1){
  args[2]=getwd()
  message("Current directory will be used")} 
 
dir=args[2]

setwd(dir)
dataa=read.table(args[1], h=T)

subMC=dataa[,1:27]
indexMC=which(rowSums(subMC[,2:ncol(subMC)])==0)
subMC=subMC[-indexMC,]
write.table(subMC, file="bfr.MC.otutab.iddef1.txt", row.names=F)
MCname=as.data.frame(subMC$X.OTU.ID) 
write.table(MCname, file="bfr.MC.seqnames", row.names=F)

subTM=dataa[,c(1,28:37)]
indexTM=which(rowSums(subTM[,2:ncol(subTM)])==0)
subTM=subTM[-indexTM,]
write.table(subTM, file="bfr.TM.otutab.iddef1.txt", row.names=F)
TMname=as.data.frame(subTM$X.OTU.ID) 
write.table(TMname, file="bfr.TM.seqnames", row.names=F)


