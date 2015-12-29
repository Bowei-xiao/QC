# Using plink input to do Heterozygousity Check
args=(commandArgs(TRUE))
fileName = as.character(args[1])
p1_auto = read.table(paste0(fileName,'_autoHet.het'),header=T,stringsAsFactors = F)
p1_sex = read.table(paste0(fileName,'_sexHet.het'),header=T,stringsAsFactors = F)
p1_auto$het = 1-p1_auto$O.HOM./p1_auto$N.NM.
p1_sex$het = 1-p1_sex$O.HOM./p1_sex$N.NM.

plotHet = function(p_auto,p_sex,file){
pOut = p_auto[p_auto$het<=mean(p_auto$het)-3*sd(p_auto$het) | p_auto$het>=mean(p_auto$het)+3*sd(p_auto$het),]
png(paste0(file,'_HetPlot.png'),width=741, height=812)
plot(p_auto$het,p_sex$het,xlab='Autosomal Chromosome heterozygousity',ylab='Sex Chromosome heterozygousity',main='Heterozygousity check')
abline(v=mean(p_auto$het))
abline(v=mean(p_auto$het)-3*sd(p_auto$het),lty=2)
abline(v=mean(p_auto$het)+3*sd(p_auto$het),lty=2)
points(pOut$het,p_sex$het[as.numeric(rownames(pOut))],col='red')
dev.off()
return(pOut)
}
p1Out = plotHet(p1_auto,p1_sex,fileName)
write.table(p1Out[,c('FID','IID')],paste0(fileName,'_indNeedremove.txt'),quote=F,col.names=F,row.names=F)
