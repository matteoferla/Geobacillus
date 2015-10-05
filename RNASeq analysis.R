#In Excel save as tsv to avoid Win new line issue. No spaces in header
datafile=file.choose()
cds = read.table( datafile, header=TRUE, row.names=1)
#Issue?
#try with ,fill=TRUE as an extra argument then
#> head(cds)

#show scatterplot of replicate vs. replicate. I need to imporve this as I did it matlab. Potentially moving after normalisation.
library(Hmisc)
par(mfrow=c(2,3))
for (i in 1:6) {
  r=rcorr(t(rbind(cds[,i*2-1]+1, cds[,i*2]+1)),type="pearson")
  plot(log2(cds[,i*2-1]+1),log2(cds[,i*2]+1),main=des$condition[i*2],sub=paste("ρ=",toString(r$r[1,2])), xlab="log2(cds+1)",ylab="log2(cds+1)",col="blue",pch=20, ylim=c(0,15), xlim=c(0,15))
  lines(c(0,20),c(0,20),col="black")
}


des=data.frame(row.names=colnames( cds),condition=c("UM5","UM5","UM1","UM1","UR5","UR5","UR1","UR1","AM5","AM5","AM1","AM1"),aerobicity=c("unaero","unaero","unaero","unaero","unaero","unaero","unaero","unaero","aero","aero","aero","aero"),richness=c("min","min","min","min","rich","rich","rich","rich","min","min","min","min"),speed=c("0.05/h","0.05/h","0.1/h","0.1/h","0.05/h","0.05/h","0.1/h","0.1/h","0.05/h","0.05/h","0.1/h","0.1/h"))
library("DESeq2")
#if not installed, google it or bioconductor.
#load all just for fun.
dds <- DESeqDataSetFromMatrix(countData = cds, colData = des,design = ~ aerobicity)
rld <- rlog(dds)
plotPCA(rld, intgroup=c("aerobicity", "richness", "speed"))
write.csv(as.data.frame(assay(rld)),file="F:\\RNASeq\\R analysis\\1_rld.csv")

library(corrgram)
distsRL <- as.matrix(dist(t(assay(rld))))

corrgram(distsRL, order=NULL, lower.panel=panel.shade,
         upper.panel=NULL, panel.txt=colnames(dds),
         main="Euclidean distance between samples")


#load parts for actual analysis
#Richness
RvM <- DESeqDataSetFromMatrix(countData = cds[,1:8], colData = des[1:8,],design = ~ richness)
RvM$aerobicity <- droplevels(RvM$aerobicity)
design(RvM) <-formula(~ richness + speed)
RvM <-DESeq(RvM)

#Alt code. Pointless test.
#RvM <-DESeq(RvM, test="LRT",full=design(RvM),reduced=formula(~ richness),betaPrior=TRUE)
#the manual says test nbinomWaldTest or nbinomLRT

resR <-results(RvM, contrast =c("richness","rich","min"))
write.csv(as.data.frame(resR),file="F:\\RNASeq\\R analysis\\1_rich.csv")
resR2 <-results(RvM, contrast =c("speed","0.05/h","0.1/h"))
write.csv(as.data.frame(resR2),file="F:\\RNASeq\\R analysis\\1_rspeed.csv")
#Aerobicity
AvU <- DESeqDataSetFromMatrix(countData = cds[,c(1:4,9:12)], colData = des[c(1:4,9:12),],design = ~ aerobicity)
AvU$richness <- droplevels(RvM$richness)
design(AvU) <-formula(~ aerobicity + speed)
AvU <-DESeq(AvU)
resA <-results(AvU, contrast =c("aerobicity","aero","unaero"))
write.csv(as.data.frame(resA),file="F:\\RNASeq\\R analysis\\1_aero.csv")
resA2 <-results(AvU, contrast =c("speed","0.05/h","0.1/h"))
write.csv(as.data.frame(resA2),file="F:\\RNASeq\\R analysis\\1_aspeed.csv")

#MA plots
par(mfrow=c(2,2))
plotMA(resR, main="Rich vs. minimal (unaerobic only)")
plotMA(resR2, main="0.1/h vs. 0.05/h in richness subset")
plotMA(resA, main="aerobic vs. unaerobic (minimal only)")
plotMA(resR2, main="0.1/h vs. 0.05/h in aerobicity subset")

#Histograms of p values
par(mfrow=c(2,3))
hist(resR$pvalue,main="Richness", xlab="p-value")
hist(resR$padj,main="Richness", xlab="BH-adj p-value")
hist(resR2$pvalue,main="Speed (Richness sub)", xlab="p-value")
hist(resA$pvalue,main="Aerobicity", xlab="p-value")
hist(resA$padj,main="Aerobicity", xlab="BH-adj p-value")
hist(resA2$pvalue,main="Speed (aerobicity sub)", xlab="p-value")

#let's compare the speeds in a single condition
Sp1 <- DESeqDataSetFromMatrix(countData = cds[,1:4], colData = des[1:4,],design = ~ speed)
Sp1$aerobicity <- droplevels(Sp1$aerobicity)
Sp1$richness <- droplevels(Sp1$richness)
Sp1 <-DESeq(Sp1) 
resS1 <-results(Sp1, contrast =c("speed","0.05/h","0.1/h"))
Sp2 <- DESeqDataSetFromMatrix(countData = cds[,5:8], colData = des[5:8,],design = ~ speed)
Sp2$aerobicity <- droplevels(Sp2$aerobicity)
Sp2$richness <- droplevels(Sp2$richness)
Sp2 <-DESeq(Sp2) 
resS2 <-results(Sp2, contrast =c("speed","0.05/h","0.1/h"))
Sp3 <- DESeqDataSetFromMatrix(countData = cds[,9:12], colData = des[9:12,],design = ~ speed)
Sp3$aerobicity <- droplevels(Sp3$aerobicity)
Sp3$richness <- droplevels(Sp3$richness)
Sp3 <-DESeq(Sp3) 
resS3 <-results(Sp3, contrast =c("speed","0.05/h","0.1/h"))
par(mfrow=c(1,3))
plotMA(resS1, main="0.1/h vs. 0.05/h in UM subset")
plotMA(resS2, main="0.1/h vs. 0.05/h in UR subset")
plotMA(resS3, main="0.1/h vs. 0.05/h in AM subset")

write.csv(as.data.frame(resS1),file="F:\\RNASeq\\6 norm & sig\\Sp1.csv")
write.csv(as.data.frame(resS2),file="F:\\RNASeq\\6 norm & sig\\Sp2.csv")
write.csv(as.data.frame(resS3),file="F:\\RNASeq\\6 norm & sig\\Sp3.csv")

