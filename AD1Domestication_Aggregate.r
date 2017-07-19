options(stringsAsFactors = FALSE)
library(DESeq2)

######
#for consolidated reads
#read in raw counts
file.names <- list.files(pattern = "*.counts")
counts.raw <- data.frame()
for (i in 1:length(file.names)){
  #open file
  temp <- read.table(file.names[i], header = FALSE, sep = "\t")
  temp <- temp[grep("Gorai",temp[,1]),]
  #merge into data.frame so all counts are together
  if (dim(counts.raw)[1] == 0){
    counts.raw <- temp
  }
  else{
    counts.raw <- merge(counts.raw, temp, by = "V1")
  }
}
row.names(counts.raw) <- counts.raw[,1]; counts.raw <- counts.raw[,-1]
colnames(counts.raw) <- sapply(strsplit(file.names,split = "[.]"),"[[",1)

#make colMat - metadata file
accession <- c(rep("AD1.CRB252",4),rep("AD1.Maxxa",4),rep("AD1.TM1",4),rep("AD1.TX2095",4),
               rep("AD1.TX665",4),rep("AD1.Yuc",4))
timepoint <- rep(c(10,15,20,5),6)
cult <- c(rep("D",12),rep("W",12))

colMat <- data.frame(accession = accession, timepoint = factor(timepoint), cult = factor(cult))
colMat$sample0 <- paste(colMat$cult,colMat$timepoint,sep = ".")

#build DDS object
AD1.ddsPair <- DESeqDataSetFromMatrix(countData = counts.raw, colData = colMat, design = ~ cult + timepoint)
AD1.ddsPair <- AD1.ddsPair[rowSums(counts(AD1.ddsPair))/24 >= 1,]
AD1.ddsPair <- AD1.ddsPair[rowSums(counts(AD1.ddsPair))/24 >= 12,]#average if half of samples had at least one read for each gene
#removing more genes actually adds variation, so stick with average one read per gene

rld <- rlog(AD1.ddsPair, blind=FALSE)
AD1.ddsPair <- estimateSizeFactors(AD1.ddsPair)
plot(assay(rld)[,1:2],pch=16, cex=0.3)

sampleDists <- dist(t(assay(rld)))
library("pheatmap")
library("RColorBrewer")

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste( rld$accession, rld$cult, rld$timepoint, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

plotPCA(rld, intgroup = c("cult", "timepoint"), ntop = 20000)
##remove Yuc 10 dpa - does not match with other 10 dpa

AD1.ddsPair <- DESeq(AD1.ddsPair)

AD1.dds0 <- DESeqDataSetFromMatrix(countData = counts.raw[,-21], colData = colMat[-21,], design = ~ sample0)
AD1.dds0 <- AD1.dds0[rowSums(counts(AD1.dds0))/24 >= 1,]
AD1.dds0 <- DESeq(AD1.dds0)

#wilds
summary(results(AD1.dds.wTX2094, contrast = c("sample0","W.5","W.10")), alpha =0.05)
summary(results(AD1.dds.wTX2094, contrast = c("sample0","W.10","W.15")), alpha =0.05)
summary(results(AD1.dds.wTX2094, contrast = c("sample0","W.15","W.20")), alpha =0.05)
#doms
summary(results(AD1.dds.wTX2094, contrast = c("sample0","D.5","D.10")), alpha =0.05)
summary(results(AD1.dds.wTX2094, contrast = c("sample0","D.10","D.15")), alpha =0.05)
summary(results(AD1.dds.wTX2094, contrast = c("sample0","D.15","D.20")), alpha =0.05)
#between wild and domesticated
summary(results(AD1.dds.wTX2094, contrast = c("sample0","W.5","D.5")), alpha =0.05)
summary(results(AD1.dds.wTX2094, contrast = c("sample0","W.10","D.10")), alpha =0.05)
summary(results(AD1.dds.wTX2094, contrast = c("sample0","W.15","D.15")), alpha =0.05)
summary(results(AD1.dds.wTX2094, contrast = c("sample0","W.20","D.20")), alpha =0.05)
#expressed genes
table(rowSums(counts(AD1.dds.wTX2094[,AD1.dds.wTX2094$sample0 %in% "W.5"])) > 1)
table(rowSums(counts(AD1.dds.wTX2094[,AD1.dds.wTX2094$sample0 %in% "W.10"])) > 1)
table(rowSums(counts(AD1.dds.wTX2094[,AD1.dds.wTX2094$sample0 %in% "W.15"])) > 1)
table(rowSums(counts(AD1.dds.wTX2094[,AD1.dds.wTX2094$sample0 %in% "W.20"])) > 1)
table(rowSums(counts(AD1.dds.wTX2094[,AD1.dds.wTX2094$sample0 %in% "D.5"])) > 1)
table(rowSums(counts(AD1.dds.wTX2094[,AD1.dds.wTX2094$sample0 %in% "D.10"])) > 1)
table(rowSums(counts(AD1.dds.wTX2094[,AD1.dds.wTX2094$sample0 %in% "D.15"])) > 1)
table(rowSums(counts(AD1.dds.wTX2094[,AD1.dds.wTX2094$sample0 %in% "D.20"])) > 1)

#Now to get out JJ's genes
gois <- read.table(file = "yang_D5_IDs.txt", header = FALSE)[,1]
#rldMat <- rlog(as.matrix(counts.raw))
#rownames(rldMat) <- rownames(counts.raw)
#goiMat <- rldMat[rownames(rldMat) %in% gois,]
#write.table(goiMat,file = "GenesOfInterest.rlog.txt", sep = "\t",row.names = TRUE, col.names = TRUE)
#gois <- unique(gois)

#wild
r.w5v10 <- results(AD1.dds.wTX2094, contrast = c("sample0","W.5","W.10"))
r.w5v10.goi <- r.w5v10[rownames(r.w5v10) %in% allgois,]
write.table(r.w5v10.goi,file = "DESeq2.wild5v10.goi.txt",sep="\t")
summary(r.w5v10.goi)

r.w10v15 <- results(AD1.dds.wTX2094, contrast = c("sample0","W.10","W.15"))
r.w10v15.goi <- r.w10v15[rownames(r.w10v15) %in% allgois,]
write.table(r.w10v15.goi,file = "DESeq2.wild10v15.goi.txt",sep="\t")
summary(r.w10v15.goi)

r.w15v20 <- results(AD1.dds.wTX2094, contrast = c("sample0","W.15","W.20"))
r.w15v20.goi <- r.w15v20[rownames(r.w15v20) %in% allgois,]
write.table(r.w15v20.goi,file = "DESeq2.wild15v20.goi.txt",sep="\t")
summary(r.w15v20.goi)

#domesticated
r.d5v10 <- results(AD1.dds.wTX2094, contrast = c("sample0","D.5","D.10"))
r.d5v10.goi <- r.d5v10[rownames(r.d5v10) %in% allgois,]
write.table(r.d5v10.goi,file = "DESeq2.dom5v10.goi.txt",sep="\t")
summary(r.d5v10.goi)

r.d10v15 <- results(AD1.dds.wTX2094, contrast = c("sample0","D.10","D.15"))
r.d10v15.goi <- r.d10v15[rownames(r.d10v15) %in% allgois,]
write.table(r.d10v15.goi,file = "DESeq2.dom10v15.goi.txt",sep="\t")
summary(r.d10v15.goi)

r.d15v20 <- results(AD1.dds.wTX2094, contrast = c("sample0","D.15","D.20"))
r.d15v20.goi <- r.d15v20[rownames(r.d15v20) %in% allgois,]
write.table(r.w15v20.goi,file = "DESeq2.dom15v20.goi.txt",sep="\t")
summary(r.d15v20.goi)

#Between wild and domesticated
r.w5vd5 <- results(AD1.dds.wTX2094, contrast = c("sample0","W.5","D.5"))
r.w5vd5.goi <- r.w5vd5[rownames(r.w5vd5) %in% allgois,]
write.table(r.w5vd5.goi,file = "DESeq2.w5vd5.goi.txt",sep="\t")
summary(r.w5vd5.goi)

r.w10vd10 <- results(AD1.dds.wTX2094, contrast = c("sample0","W.10","D.10"))
r.w10vd10.goi <- r.w10vd10[rownames(r.w10vd10) %in% allgois,]
write.table(r.w10vd10.goi,file = "DESeq2.w10vd10.goi.txt",sep="\t")
summary(r.w10vd10.goi)

r.w15vd15 <- results(AD1.dds.wTX2094, contrast = c("sample0","W.15","D.15"))
r.w15vd15.goi <- r.w15vd15[rownames(r.w15vd15) %in% allgois,]
write.table(r.w15vd15.goi,file = "DESeq2.w15vd15.goi.txt",sep="\t")
summary(r.w15vd15.goi)

r.w20vd20 <- results(AD1.dds.wTX2094, contrast = c("sample0","W.20","D.20"))
r.w20vd20.goi <- r.w20vd20[rownames(r.w20vd20) %in% allgois,]
write.table(r.w20vd20.goi,file = "DESeq2.w20vd20.goi.txt",sep="\t")
summary(r.w20vd20.goi)

###
#use edgeR to get RPKM for the the partitioned reads
#already counted with HTseq based on BamBam-made bam files
source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
library(edgeR)

#read in count files, as for unpartitioned reads
file.names <- list.files(pattern = "*.A.counts")
counts.A.raw <- data.frame()
for (i in 1:length(file.names)){
  #open file
  temp <- read.table(file.names[i], header = FALSE, sep = "\t")
  temp <- temp[grep("Gorai",temp[,1]),]
  #merge into data.frame so all counts are together
  if (dim(counts.A.raw)[1] == 0){
    counts.A.raw <- temp
  }
  else{
    counts.A.raw <- merge(counts.A.raw, temp, by = "V1")
  }
}
row.names(counts.A.raw) <- counts.A.raw[,1]; counts.A.raw <- counts.A.raw[,-1]
colnames(counts.A.raw) <- sapply(strsplit(file.names,split = "[.]"),"[[",1)
#for A and D!

DGE.A.r <- DGEList(counts.A.raw)
groups <- factor(colMat$sample0)
DGE.A.r <- DGEList(counts.A.raw, group = groups)
DGE.A.r <- calcNormFactors(DGE.A.r, na.rm = TRUE)
geneLengths <- read.table("D5geneLengths.txt",header = FALSE)

#this is giving me trouble, but I forgot how easy it is to just calculate RPKM
#We can just write the code 
#Number of reads mapped divided by length of gene in kbases divided by library size over 1 million

#calculate "per million reads mapped"
libSizeD <- colSums(counts.D.raw)
milLibSizeD <- libSizeD/10^6
libSizeA <- colSums(counts.A.raw)
milLibSizeA <- libSizeA/10^6

#calculate "reads per killobase"
geneLengths.kb <- geneLengths/10^3
gL.kb <- geneLengths.kb[,1]
rpkm.A <- sweep((counts.A.raw/gL.kb),2,milLibSizeA,"/")
rpkm.D <- sweep((counts.D.raw/gL.kb),2,milLibSizeD,"/")

#Peng Liu suggested that Fisher's exact test comparing the RPKM of each homoeolog
#was statistically appropriate.
#So FET is given by fisher.test()
#need to make a 2 by 2 matrix for each comparison
#
#	[A homoeolog count]	[D homoeolog count]
#
#	[A library size]		[D library size]
#
testMat <- matrix(nrow=2,ncol=2)
#testing the matrix function
matrix(data = c(1,2,3,4),nrow=2,ncol=2)
#makes
# 1	3
# 2	4
#
#make a for loop to go through all the comparisons
#or should I do this with something like sweep
#on reddit they suggest using apply like so

#try this again with counts
fisherMat <- matrix()
for(i in 1:dim(rpkm.A)[2]){
	matTemp <- cbind(rpkm.A[,i],rep(milLibSizeA[i],times=dim(rpkm.A)[1]),rpkm.D[,i],rep(milLibSizeD[i],times=dim(rpkm.D)[1]))
	row.names(matTemp) <- NULL
	fisherMat <- t(apply(matTemp,1, function(x) unlist(fisher.test(matrix(x,nr=2))[c(3,1)])))
	print(c(colnames(rpkm.A)[i],sum(fisherMat[,2] < 0.05, na.rm = TRUE)))
	write.table(fisherMat, file = paste(colnames(rpkm.A)[i],".homoeologBiasFisher.txt",sep = ''), row.names = rownames(counts.A.raw), col.names = FALSE)
}

##Get total number of DE genes
rnames <- vector("list", length = 6)
for(i in 1:6){
	temp = comparisons[[i]][!is.na(comparisons[[i]]$padj),];
	rnames[[i]] <- rownames(temp[temp$padj < 0.05,])
}
allSigRows <- unlist(rnames)
unique(allSigRows)

##Get list of DE genes sorted by adj pvalue
results.list <- vector("list",length = 10)
results.list[[1]] <- r.w5v10[!is.na(r.w5v10$padj),]
results.list[[2]] <- r.w10v15[!is.na(r.w10v15$padj),]
results.list[[3]] <- r.w15v20[!is.na(r.w15v20$padj),]
results.list[[4]] <- r.d5v10[!is.na(r.d5v10$padj),]
results.list[[5]] <- r.d10v15[!is.na(r.d10v15$padj),]
results.list[[6]] <- r.d15v20[!is.na(r.d15v20$padj),]
results.list[[7]] <- r.w5vd5[!is.na(r.w5vd5$padj),]
results.list[[8]] <- r.w10vd10[!is.na(r.w10vd10$padj),]
results.list[[9]] <- r.w15vd15[!is.na(r.w15vd15$padj),]
results.list[[10]] <- r.w20vd20[!is.na(r.w20vd20$padj),]

for(k in 1:length(results.list)){
	tempSig <- results.list[[k]][results.list[[k]]$padj < 0.05,]
	upSig <- tempSig[tempSig$log2FoldChange > 0,]
	downSig <- tempSig[tempSig$log2FoldChange < 0,]
	assign(paste("results_",k,"_upreg",sep=""),upSig)
	assign(paste("results_",k,"_downreg",sep=""),downSig)
	write.table(x = upSig,file = paste("results_",k,"_upSig.txt",sep=""))
	write.table(x = downSig,file = paste("results_",k,"_downSig.txt",sep=""))
}

##setting up for topGO
library(topGO)
#read in mapping of GO terms
geneID2GO <- readMappings(file = "G.raimondii_JGI_221_v2.1_genes2GO.map", sep = " ", IDsep =",")
#read in gene universe
de.results <- list.files(pattern = "results")
geneNames <- names(geneID2GO)
GOresults<-data.frame()

for(m in 1:length(de.results)){
	de.table <- read.table(file = de.results[m])#starts with results 10, then 1,2,3,...; results 10 is w20v d20, results 1 is w5vw10
	sig.list <- rownames(de.table)
	geneList <- factor(as.integer(geneNames %in% sig.list))
	names(geneList) <- geneNames
	
	pdf(file=paste("topGO/DEresult",m,".pdf", sep=""))
    # topGO analysis
    remove(enrich)
	for(on in c("MF","BP","CC"))
    {
		GOdata <- new("topGOdata", ontology = on, allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
		result <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
		results.table <- GenTable(GOdata, result, topNodes = length(result@score))
		# reduce results to GO terms passing Benjamini-Hochberg multiple hypothesis corrected pval <= 0.05, FDR <= 5%
		results.table$qval.bh<-p.adjust(results.table[,"result1"],method="BH")
		# label ontology type
		results.table$ontology<-on
		# reduce results to GO terms passing Benjamini-Hochberg multiple hypothesis corrected pval <= 0.05, consider FDR <= 5% in future
		keep <- results.table[as.numeric(results.table[,"qval.bh"])<0.05,]
		if(exists("enrich")) {enrich<- rbind(enrich, keep)}
        if(!exists("enrich")) {enrich<- keep}
        
        # draw figure for GO terms pval<=0.05 before FDR correction
        if(is.na(sigNo<-length(keep$ontology))){next}
        showSigOfNodes(GOdata, score(result), firstSigNodes = sigNo, useInfo = "all")
        mtext(on, line=-1)
    }
	dev.off()
    if(dim(enrich)[1]>0)
    {
        enrichDE <- enrich
		enrichDE$DE <- m
        GOresults<-rbind(GOresults,enrichDE)
	}
}

write.table(GOresults, file="topGO/GOresults.txt", sep="\t", row.names=FALSE)

##Use this for orders of results
#> de.results
# [1] "results_10_downSig.txt" "results_10_upSig.txt"   "results_1_downSig.txt"
# [4] "results_1_upSig.txt"    "results_2_downSig.txt"  "results_2_upSig.txt"
# [7] "results_3_downSig.txt"  "results_3_upSig.txt"    "results_4_downSig.txt"
#[10] "results_4_upSig.txt"    "results_5_downSig.txt"  "results_5_upSig.txt"
#[13] "results_6_downSig.txt"  "results_6_upSig.txt"    "results_7_downSig.txt"
#[16] "results_7_upSig.txt"    "results_8_downSig.txt"  "results_8_upSig.txt"
#[19] "results_9_downSig.txt"  "results_9_upSig.txt"

DE.aggr.devAll <- list(length = 4)
DE.aggr.devAll[[1]] <- results(AD1.aggr.wTX.ddsCultAndDpa, contrast = c("timepoint","5","10"))
DE.aggr.devAll[[2]] <- results(AD1.aggr.wTX.ddsCultAndDpa, contrast = c("timepoint","10","15"))
DE.aggr.devAll[[3]] <- results(AD1.aggr.wTX.ddsCultAndDpa, contrast = c("timepoint","15","20"))
DE.aggr.devAll[[4]] <- results(AD1.aggr.wTX.ddsCultAndDpa, contrast = c("cult","W","D"))

#Add for modules of interest
GOresults<-data.frame()
shortLabels = c("15","26")

for(m in 1:length(modsInterest)){
	sig.list <- modsInterest[[m]]
	geneList <- factor(as.integer(geneNames %in% sig.list))
	names(geneList) <- geneNames
	pdf(file=paste("topGO/ModResult/moduleOfInterest",shortLabels[m],".pdf", sep=""))
    # topGO analysis
    if(exists("enrich")) {remove(enrich)}
	for(on in c("MF","BP","CC"))
    {
		GOdata <- new("topGOdata", ontology = on, allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
		result <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
		results.table <- GenTable(GOdata, result, topNodes = length(result@score))
		# reduce results to GO terms passing Benjamini-Hochberg multiple hypothesis corrected pval <= 0.05, FDR <= 5%
		results.table$qval.bh<-p.adjust(results.table[,"result1"],method="BH")
		# label ontology type
		results.table$ontology<-on
		# reduce results to GO terms passing Benjamini-Hochberg multiple hypothesis corrected pval <= 0.05, consider FDR <= 5% in future
		keep <- results.table[as.numeric(results.table[,"qval.bh"])<0.05,]
		if(exists("enrich")) {enrich<- rbind(enrich, keep)}
        if(!exists("enrich")) {enrich<- keep}
        
        # draw figure for GO terms pval<=0.05 before FDR correction
        if(is.na(sigNo<-length(keep$ontology))){next}
        showSigOfNodes(GOdata, score(result), firstSigNodes = sigNo, useInfo = "all")
        mtext(on, line=-1)
    }
	dev.off()
    if(dim(enrich)[1]>0)
    {
        enrichMod <- enrich
		enrichMod$mod <- shortLabels[[m]]
		GOresults<-rbind(GOresults,enrichMod)
	}
}

write.table(GOresults, file="topGO/ModResult/GOresults.txt", sep="\t", row.names=FALSE)

#making a better PCA graph
library(genefilter)
sumPCA<-
function (x, intgroup = "condition", ntop = 500)
{
    rv = rowVars(assay(x))
    select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca = prcomp(t(assay(x)[select, ]))
    fac = factor(apply(as.data.frame(colData(x)[, intgroup, drop = FALSE]),
    1, paste, collapse = " : "))
    return(pca)
}
data<-sumPCA(rld, intgroup = c("cult", "timepoint"))
summary(data)

#Importance of components:
#                           PC1     PC2      PC3     PC4     PC5     PC6     PC7    PC8     PC9    PC10    PC11    PC12    PC13    PC14    PC15    PC16    PC17
#Standard deviation     53.1896 15.5827 12.20257 10.1120 9.07308 7.99356 6.28646 5.8493 5.67010 4.79848 4.53138 3.91018 3.73555 3.29760 3.13308 2.95674 2.78700
#Proportion of Variance  0.7608  0.0653  0.04004  0.0275 0.02214 0.01718 0.01063 0.0092 0.00865 0.00619 0.00552 0.00411 0.00375 0.00292 0.00264 0.00235 0.00209
#Cumulative Proportion   0.7608  0.8261  0.86613  0.8936 0.91576 0.93294 0.94357 0.9528 0.96142 0.96761 0.97313 0.97724 0.98099 0.98392 0.98656 0.98891 0.99100
#                          PC18    PC19    PC20    PC21    PC22    PC23    PC24      PC25
#Standard deviation     2.64286 2.52994 2.29554 2.15574 2.04137 1.79675 1.66628 1.009e-14
#Proportion of Variance 0.00188 0.00172 0.00142 0.00125 0.00112 0.00087 0.00075 0.000e+00
#Cumulative Proportion  0.99288 0.99460 0.99601 0.99726 0.99839 0.99925 1.00000 1.000e+00
library(ggplot2)
qplot(PC1, PC2, main="PCA of top 20000 gene counts", color=colMat.wTX2094$timepoint, shape=colMat.wTX2094$cult, data=as.data.frame(data$x)) +
geom_point(size = 3) +
xlab(paste0("PC1: 76.1% variance")) +
ylab(paste0("PC2: 6.53% variance")) +
scale_colour_discrete(name = "DPA", breaks = c(5,10,15,20), labels = c("5","10","15","20")) +
scale_shape_discrete(name = "Status", breaks = c("W","D"), labels = c("wild","domesticated"))
ggsave("AD1.rld.allSamples.png", width = 10, height = 5)

#looking at homoeolog bias in modules of aggregate network
#GOAL: homoeolog bias by module!
Allnet12Modules <- data.frame(Gene = colnames(multiExpr[[3]]$data),Modules = Allnet12$colors)
temp <- log(counts.raw.ratio)
temp[is.nan(as.matrix(temp))] <- NA
temp[is.infinite(as.matrix(temp))] <- NA
rowMeans(temp)
rowMeans(temp, na.rm = TRUE)
all.ADRatio.means <- rowMeans(temp, na.rm =TRUE)
all.ADRatio.means.matched <- all.ADRatio.means[names(all.ADRatio.means) %in% Allnet12Modules$Gene]
AggrAllModsAndHomoeoBias <- cbind(Allnet12Modules,all.ADRatio.means.matched)

#Use code from the all
for(i in 1:(length(unique(AggrAllModsAndHomoeoBias$Modules))-1)){
	temp <- subset(AggrAllModsAndHomoeoBias,Modules %in% i)
	ggplot(temp,aes(x = Gene, y = LogADRatio)) +
	geom_point() +
	ggtitle(paste("Module",i,"Average Homoeolog Bias",sep = " ")) +
	geom_hline(yintercept = 0, colour = "red")
	ggsave(paste("Module",i,"AverageLogHomoeoBias.pdf",sep = ""),device = "pdf")
}
#problem is that this averages, loses any bias in wild/dom

#Is it DE in both A and D copies
table(table(gsub(pattern = "_[AD]", rep="", unique(c(DEgenes[[1]],DEgenes[[2]],DEgenes[[3]])))) == 2)