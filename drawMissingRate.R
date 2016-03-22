# read in the imiss and lmiss files and plot the missingness for both
args=(commandArgs(TRUE))
fileName = as.character(args[1])
iMiss = read.table(paste0(fileName,'.imiss'),header=T,stringsAsFactors=F)
lMiss = read.table(paste0(fileName,'.lmiss'),header=T,stringsAsFactors=F)

png(paste0(fileName,'_IndMissing.png'),width=741,height=812)
hist(iMiss$F_MISS,labels=T,xlab='Missing Rate',main='Indidual Missing Rate')
dev.off()

png(paste0(fileName,'_SNPMissing.png'),width=741,height=812)
hist(lMiss$F_MISS,labels=T,xlab='Missing Rate',main='SNP Missing Rate')
dev.off()
