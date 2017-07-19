options(stringsAsFactors = FALSE)
library(DESeq2)
##for partitioned reads
#read in raw counts
file.names.A <- list.files(pattern = "*.A.counts")
file.names.D <- list.files(pattern = "*.D.counts")
counts.A.raw <- data.frame()
counts.D.raw <- data.frame()
##for A genome
for (i in 1:length(file.names.A)){
  #open file
  temp <- read.table(file.names.A[i], header = FALSE, sep = "\t")
  temp <- temp[grepl("Gorai.0",temp[,1]),]
  #merge into data.frame so all counts are together
  if (dim(counts.A.raw)[1] == 0){
    counts.A.raw <- temp
  }
  else{
    counts.A.raw <- merge(counts.A.raw, temp, by = "V1")
  }
}
row.names(counts.A.raw) <- counts.A.raw[,1]; counts.A.raw <- counts.A.raw[,-1]
colnames(counts.A.raw) <- sapply(strsplit(file.names.A,split = "[.]"),"[[",1)

##for D genome
for (i in 1:length(file.names.D)){
  #open file
  temp <- read.table(file.names.D[i], header = FALSE, sep = "\t")
  temp <- temp[grepl("Gorai.0",temp[,1]),]
  #merge into data.frame so all counts are together
  if (dim(counts.D.raw)[1] == 0){
    counts.D.raw <- temp
  }
  else{
    counts.D.raw <- merge(counts.D.raw, temp, by = "V1")
  }
}
row.names(counts.D.raw) <- counts.D.raw[,1]; counts.D.raw <- counts.D.raw[,-1]
colnames(counts.D.raw) <- sapply(strsplit(file.names.D,split = "[.]"),"[[",1)

##put them together
#reorder first
counts.A.raw <- counts.A.raw[,c(4,1:3,8,5:7,12,9:11,17,14:16,21,18:20,24,13,22,23)]
counts.D.raw <- counts.D.raw[,c(4,1:3,8,5:7,12,9:11,17,14:16,21,18:20,24,13,22,23)]
#rename with appropriate homoeolog identifiers
A.names <- paste(rownames(counts.A.raw),"A",sep="_")
D.names <- paste(rownames(counts.D.raw),"D",sep="_")

rownames(counts.A.raw) <- A.names
rownames(counts.D.raw) <- D.names
#put em together
counts.AD.raw <- rbind(counts.A.raw,counts.D.raw)

#make colMat - metadata file
#we will use Yuc for TX2094, it is true, but just note that
accession <- c(rep("AD1.CRB252",4),rep("AD1.Maxxa",4),rep("AD1.TM1",4),rep("AD1.TX2095",4),
               rep("AD1.TX665",4),rep("AD1.Yuc",4))
timepoint <- rep(c(5,10,15,20),6)
cult <- c(rep("D",12),rep("W",12))
colMat <- data.frame(accession = accession, timepoint = factor(timepoint), cult = factor(cult))
colMat$sample0 <- paste(colMat$cult,colMat$timepoint,sep = ".")

#Start DESeq2 process
AD1.AD.dds0 <- DESeqDataSetFromMatrix(countData = counts.AD.raw, colData = colMat, design = ~ sample0)
AD1.AD.ddsCultDev <- DESeqDataSetFromMatrix(countData = counts.AD.raw, colData = colMat, design = ~ cult + timepoint)
#filter out genes with no counts
AD1.AD.dds0 <- AD1.AD.dds0[rowSums(counts(AD1.AD.dds0))/24 >= 1,]
AD1.AD.ddsCultDev <- AD1.AD.ddsCultDev[rowSums(counts(AD1.AD.ddsCultDev))/24 >= 1,]
AD1.AD.dds0 <- estimateSizeFactors(AD1.AD.dds0)
AD1.AD.ddsCultDev <- estimateSizeFactors(AD1.AD.ddsCultDev)

#build PCA to confirm proper samples
AD.rld <- rlog(AD1.AD.dds0, blind=FALSE)
pdf("AD1.partitionedReads.rlog.plots.pdf")
plot(assay(AD.rld)[,1:2],pch=16, cex=0.3)

sampleDists <- dist(t(assay(AD.rld)))
library("pheatmap")
library("RColorBrewer")

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste( AD.rld$accession, AD.rld$cult, AD.rld$timepoint, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

plotPCA(AD.rld, intgroup = c("cult", "timepoint"), ntop = 20000)
dev.off()

##ok, everything looks fine, the same as with aggregate
AD1.AD.dds0 <- DESeq(AD1.AD.dds0)
AD1.AD.ddsCultDev <- DESeq(AD1.AD.ddsCultDev)

#wilds
summary(results(AD1.AD.dds0, contrast = c("sample0","W.5","W.10")), alpha =0.05)
summary(results(AD1.AD.dds0, contrast = c("sample0","W.10","W.15")), alpha =0.05)
summary(results(AD1.AD.dds0, contrast = c("sample0","W.15","W.20")), alpha =0.05)
#doms
summary(results(AD1.AD.dds0, contrast = c("sample0","D.5","D.10")), alpha =0.05)
summary(results(AD1.AD.dds0, contrast = c("sample0","D.10","D.15")), alpha =0.05)
summary(results(AD1.AD.dds0, contrast = c("sample0","D.15","D.20")), alpha =0.05)
#between wild and domesticated
summary(results(AD1.AD.dds0, contrast = c("sample0","W.5","D.5")), alpha =0.05)
summary(results(AD1.AD.dds0, contrast = c("sample0","W.10","D.10")), alpha =0.05)
summary(results(AD1.AD.dds0, contrast = c("sample0","W.15","D.15")), alpha =0.05)
summary(results(AD1.AD.dds0, contrast = c("sample0","W.20","D.20")), alpha =0.05)
#expressed genes
table(rowSums(counts(AD1.AD.dds0[,AD1.AD.dds0$sample0 %in% "W.5"])) > 1)
table(rowSums(counts(AD1.AD.dds0[,AD1.AD.dds0$sample0 %in% "W.10"])) > 1)
table(rowSums(counts(AD1.AD.dds0[,AD1.AD.dds0$sample0 %in% "W.15"])) > 1)
table(rowSums(counts(AD1.AD.dds0[,AD1.AD.dds0$sample0 %in% "W.20"])) > 1)
table(rowSums(counts(AD1.AD.dds0[,AD1.AD.dds0$sample0 %in% "D.5"])) > 1)
table(rowSums(counts(AD1.AD.dds0[,AD1.AD.dds0$sample0 %in% "D.10"])) > 1)
table(rowSums(counts(AD1.AD.dds0[,AD1.AD.dds0$sample0 %in% "D.15"])) > 1)
table(rowSums(counts(AD1.AD.dds0[,AD1.AD.dds0$sample0 %in% "D.20"])) > 1)

#get parsed counts
#wild
table(rowSums(counts(AD1.AD.dds0[grep("_A",rownames(AD1.AD.dds0)),AD1.AD.dds0$sample0 %in% "W.5"])) > 1)
table(rowSums(counts(AD1.AD.dds0[grep("_D",rownames(AD1.AD.dds0)),AD1.AD.dds0$sample0 %in% "W.5"])) > 1)
table(rowSums(counts(AD1.AD.dds0[grep("_A",rownames(AD1.AD.dds0)),AD1.AD.dds0$sample0 %in% "W.10"])) > 1)
table(rowSums(counts(AD1.AD.dds0[grep("_D",rownames(AD1.AD.dds0)),AD1.AD.dds0$sample0 %in% "W.10"])) > 1)
table(rowSums(counts(AD1.AD.dds0[grep("_A",rownames(AD1.AD.dds0)),AD1.AD.dds0$sample0 %in% "W.15"])) > 1)
table(rowSums(counts(AD1.AD.dds0[grep("_D",rownames(AD1.AD.dds0)),AD1.AD.dds0$sample0 %in% "W.15"])) > 1)
table(rowSums(counts(AD1.AD.dds0[grep("_A",rownames(AD1.AD.dds0)),AD1.AD.dds0$sample0 %in% "W.20"])) > 1)
table(rowSums(counts(AD1.AD.dds0[grep("_D",rownames(AD1.AD.dds0)),AD1.AD.dds0$sample0 %in% "W.20"])) > 1)
#domesticated
table(rowSums(counts(AD1.AD.dds0[grep("_A",rownames(AD1.AD.dds0)),AD1.AD.dds0$sample0 %in% "D.5"])) > 1)
table(rowSums(counts(AD1.AD.dds0[grep("_D",rownames(AD1.AD.dds0)),AD1.AD.dds0$sample0 %in% "D.5"])) > 1)
table(rowSums(counts(AD1.AD.dds0[grep("_A",rownames(AD1.AD.dds0)),AD1.AD.dds0$sample0 %in% "D.10"])) > 1)
table(rowSums(counts(AD1.AD.dds0[grep("_D",rownames(AD1.AD.dds0)),AD1.AD.dds0$sample0 %in% "D.10"])) > 1)
table(rowSums(counts(AD1.AD.dds0[grep("_A",rownames(AD1.AD.dds0)),AD1.AD.dds0$sample0 %in% "D.15"])) > 1)
table(rowSums(counts(AD1.AD.dds0[grep("_D",rownames(AD1.AD.dds0)),AD1.AD.dds0$sample0 %in% "D.15"])) > 1)
table(rowSums(counts(AD1.AD.dds0[grep("_A",rownames(AD1.AD.dds0)),AD1.AD.dds0$sample0 %in% "D.20"])) > 1)
table(rowSums(counts(AD1.AD.dds0[grep("_D",rownames(AD1.AD.dds0)),AD1.AD.dds0$sample0 %in% "D.20"])) > 1)

#wild
r.w5v10 <- results(AD1.AD.dds0, contrast = c("sample0","W.5","W.10"))
r.w5v10.goi <- r.w5v10[rownames(r.w5v10) %in% allgois,]
write.table(r.w5v10.goi,file = "DESeq2.partitioned.wild5v10.goi.txt",sep="\t")
summary(r.w5v10.goi, alpha = 0.05)

r.w10v15 <- results(AD1.AD.dds0, contrast = c("sample0","W.10","W.15"))
r.w10v15.goi <- r.w10v15[rownames(r.w10v15) %in% allgois,]
write.table(r.w10v15.goi,file = "DESeq2.partitioned.wild10v15.goi.txt",sep="\t")
summary(r.w10v15.goi, alpha = 0.05)

r.w15v20 <- results(AD1.AD.dds0, contrast = c("sample0","W.15","W.20"))
r.w15v20.goi <- r.w15v20[rownames(r.w15v20) %in% allgois,]
write.table(r.w15v20.goi,file = "DESeq2.partitioned.wild15v20.goi.txt",sep="\t")
summary(r.w15v20.goi, alpha = 0.05)

#domesticated
r.d5v10 <- results(AD1.AD.dds0, contrast = c("sample0","D.5","D.10"))
r.d5v10.goi <- r.d5v10[rownames(r.d5v10) %in% allgois,]
write.table(r.d5v10.goi,file = "DESeq2.partitioned.dom5v10.goi.txt",sep="\t")
summary(r.d5v10.goi, alpha = 0.05)

r.d10v15 <- results(AD1.AD.dds0, contrast = c("sample0","D.10","D.15"))
r.d10v15.goi <- r.d10v15[rownames(r.d10v15) %in% allgois,]
write.table(r.d10v15.goi,file = "DESeq2.partitioned.dom10v15.goi.txt",sep="\t")
summary(r.d10v15.goi, alpha = 0.05)

r.d15v20 <- results(AD1.AD.dds0, contrast = c("sample0","D.15","D.20"))
r.d15v20.goi <- r.d15v20[rownames(r.d15v20) %in% allgois,]
write.table(r.w15v20.goi,file = "DESeq2.partitioned.dom15v20.goi.txt",sep="\t")
summary(r.d15v20.goi, alpha = 0.05)

#Between wild and domesticated
r.w5vd5 <- results(AD1.AD.dds0, contrast = c("sample0","W.5","D.5"))
r.w5vd5.goi <- r.w5vd5[rownames(r.w5vd5) %in% allgois,]
write.table(r.w5vd5.goi,file = "DESeq2.partitioned.w5vd5.goi.txt",sep="\t")
summary(r.w5vd5.goi, alpha = 0.05)

r.w10vd10 <- results(AD1.AD.dds0, contrast = c("sample0","W.10","D.10"))
r.w10vd10.goi <- r.w10vd10[rownames(r.w10vd10) %in% allgois,]
write.table(r.w10vd10.goi,file = "DESeq2.partitioned.w10vd10.goi.txt",sep="\t")
summary(r.w10vd10.goi, alpha = 0.05)

r.w15vd15 <- results(AD1.AD.dds0, contrast = c("sample0","W.15","D.15"))
r.w15vd15.goi <- r.w15vd15[rownames(r.w15vd15) %in% allgois,]
write.table(r.w15vd15.goi,file = "DESeq2.partitioned.w15vd15.goi.txt",sep="\t")
summary(r.w15vd15.goi, alpha = 0.05)

r.w20vd20 <- results(AD1.AD.dds0, contrast = c("sample0","W.20","D.20"))
r.w20vd20.goi <- r.w20vd20[rownames(r.w20vd20) %in% allgois,]
write.table(r.w20vd20.goi,file = "DESeq2.partitioned.w20vd20.goi.txt",sep="\t")
summary(r.w20vd20.goi, alpha = 0.05)

##set up for topGO
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
geneID2GO <- readMappings(file = "D5_genes2GO.AD.map", sep = " ", IDsep =",")##make sure to reset this file so that gene names have _A and _D after them
#read in gene universe
de.results <- list.files(pattern = "results")
geneNames <- names(geneID2GO)
GOresults<-data.frame()

for(m in 1:length(de.results)){
	de.table <- read.table(file = de.results[m])
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

#This will get DE partitioned into A and D with direction
#length(grep("_A",rownames(subset(r.w5v10, log2FoldChange > 0 & padj < 0.05 & baseMean > 3))))
#[1] 118
#> length(grep("_A",rownames(subset(r.w5v10, log2FoldChange < 0 & padj < 0.05 & baseMean > 3))))
#[1] 589
#> length(grep("_D",rownames(subset(r.w5v10, log2FoldChange > 0 & padj < 0.05 & baseMean > 3))))
#[1] 139
#> length(grep("_D",rownames(subset(r.w5v10, log2FoldChange < 0 & padj < 0.05 & baseMean > 3))))
#[1] 576
#> length(grep("_A",rownames(subset(r.w10v15, log2FoldChange > 0 & padj < 0.05 & baseMean > 3))))
#[1] 215
#> length(grep("_A",rownames(subset(r.w10v15, log2FoldChange < 0 & padj < 0.05 & baseMean > 3))))
#[1] 367
#> length(grep("_D",rownames(subset(r.w10v15, log2FoldChange > 0 & padj < 0.05 & baseMean > 3))))
#[1] 206
#> length(grep("_D",rownames(subset(r.w10v15, log2FoldChange < 0 & padj < 0.05 & baseMean > 3))))
#[1] 392
#> length(grep("_A",rownames(subset(r.w15v20, log2FoldChange > 0 & padj < 0.05 & baseMean > 3))))
#[1] 1580
#> length(grep("_A",rownames(subset(r.w15v20, log2FoldChange < 0 & padj < 0.05 & baseMean > 3))))
#[1] 1371
#> length(grep("_D",rownames(subset(r.w15v20, log2FoldChange > 0 & padj < 0.05 & baseMean > 3))))
#[1] 1625
#> length(grep("_D",rownames(subset(r.w15v20, log2FoldChange < 0 & padj < 0.05 & baseMean > 3))))
#[1] 1380
#> length(grep("_A",rownames(subset(r.d5v10, log2FoldChange > 0 & padj < 0.05 & baseMean > 3))))
#[1] 256
#> length(grep("_A",rownames(subset(r.d5v10, log2FoldChange < 0 & padj < 0.05 & baseMean > 3))))
#[1] 667
#> length(grep("_D",rownames(subset(r.d5v10, log2FoldChange > 0 & padj < 0.05 & baseMean > 3))))
#[1] 258
#> length(grep("_D",rownames(subset(r.d5v10, log2FoldChange < 0 & padj < 0.05 & baseMean > 3))))
#[1] 653
#> length(grep("_A",rownames(subset(r.d10v15, log2FoldChange > 0 & padj < 0.05 & baseMean > 3))))
#[1] 328
#> length(grep("_A",rownames(subset(r.d10v15, log2FoldChange < 0 & padj < 0.05 & baseMean > 3))))
#[1] 368
#> length(grep("_D",rownames(subset(r.d10v15, log2FoldChange > 0 & padj < 0.05 & baseMean > 3))))
#[1] 323
#> length(grep("_D",rownames(subset(r.d10v15, log2FoldChange < 0 & padj < 0.05 & baseMean > 3))))
#[1] 354
#> length(grep("_A",rownames(subset(r.d15v20, log2FoldChange > 0 & padj < 0.05 & baseMean > 3))))
#[1] 1796
#> length(grep("_A",rownames(subset(r.d15v20, log2FoldChange < 0 & padj < 0.05 & baseMean > 3))))
#[1] 2048
#> length(grep("_D",rownames(subset(r.d15v20, log2FoldChange > 0 & padj < 0.05 & baseMean > 3))))
#[1] 1895
#> length(grep("_D",rownames(subset(r.d15v20, log2FoldChange < 0 & padj < 0.05 & baseMean > 3))))
#[1] 2101
#> length(grep("_A",rownames(subset(r.w5vd5, log2FoldChange > 0 & padj < 0.05 & baseMean > 3))))
#[1] 290
#> length(grep("_D",rownames(subset(r.w5vd5, log2FoldChange > 0 & padj < 0.05 & baseMean > 3))))
#[1] 264
#> length(grep("_A",rownames(subset(r.w5vd5, log2FoldChange < 0 & padj < 0.05 & baseMean > 3))))
#[1] 164
#> length(grep("_D",rownames(subset(r.w5vd5, log2FoldChange < 0 & padj < 0.05 & baseMean > 3))))
#[1] 161
#> length(grep("_A",rownames(subset(r.w10vd10, log2FoldChange > 0 & padj < 0.05 & baseMean > 3))))
#[1] 362
#> length(grep("_D",rownames(subset(r.w10vd10, log2FoldChange > 0 & padj < 0.05 & baseMean > 3))))
#[1] 297
#> length(grep("_A",rownames(subset(r.w10vd10, log2FoldChange < 0 & padj < 0.05 & baseMean > 3))))
#[1] 188
#> length(grep("_D",rownames(subset(r.w10vd10, log2FoldChange < 0 & padj < 0.05 & baseMean > 3))))
#[1] 183
#> length(grep("_A",rownames(subset(r.w15vd15, log2FoldChange > 0 & padj < 0.05 & baseMean > 3))))
#[1] 274
#> length(grep("_D",rownames(subset(r.w15vd15, log2FoldChange > 0 & padj < 0.05 & baseMean > 3))))
#[1] 265
#> length(grep("_A",rownames(subset(r.w15vd15, log2FoldChange < 0 & padj < 0.05 & baseMean > 3))))
#[1] 189
#> length(grep("_D",rownames(subset(r.w15vd15, log2FoldChange < 0 & padj < 0.05 & baseMean > 3))))
#[1] 170
#> length(grep("_A",rownames(subset(r.w20vd20, log2FoldChange > 0 & padj < 0.05 & baseMean > 3))))
#[1] 453
#> length(grep("_D",rownames(subset(r.w20vd20, log2FoldChange > 0 & padj < 0.05 & baseMean > 3))))
#[1] 449
#> length(grep("_A",rownames(subset(r.w20vd20, log2FoldChange < 0 & padj < 0.05 & baseMean > 3))))
#[1] 439
#> length(grep("_D",rownames(subset(r.w20vd20, log2FoldChange < 0 & padj < 0.05 & baseMean > 3))))
#[1] 396
##need to set the same basemean threshold across all, it is fluctuating if you just run it with results() and summary()

#Working on adding lines to Fig 2; we need homoeolog biases (all, dom and wild?)
library(lattice)
counts.raw.ratio <- counts.A.raw / counts.D.raw
counts.gg.ratio <- data.frame()
keep <- data.frame()
for(i in 1:dim(counts.raw.ratio)[2]){
	temp_cols <- data.frame(rownames(counts.raw.ratio),as.double(counts.raw.ratio[,i]),rep(colMat$cult[i],times=length(rownames(counts.raw.ratio))),
	rep(colMat$timepoint[i],times=length(rownames(counts.raw.ratio))))
	temp_cols <- temp_cols[!is.nan(as.numeric(temp_cols[,2])),]
	temp_cols <- temp_cols[!is.infinite(as.numeric(temp_cols[,2])),]
	keep <- rbind(keep, temp_cols)
}
counts.gg.ratio <- keep
names(counts.gg.ratio) <- c("Gene","AD.Ratio","Cult","DPA")
pdf("HomoeologExprBias.pdf")
ggplot(counts.gg.ratio,aes(x = Gene, y = AD.Ratio, colour = Cult)) + geom_point()
dev.off()

pdf("HomoeologExprBias_LogScale.pdf")
ggplot(counts.gg.ratio,aes(x = Gene, y = AD.LogRatio, colour = Cult)) + geom_point()
dev.off()
#ttest for all gene ratios
t.test(counts.gg.ratioAndMods$AD.LogRatio)
#Error in if (stderr < 10 * .Machine$double.eps * abs(mx)) stop("data are essentially constant") :
#  missing value where TRUE/FALSE needed
#Numbers are too small, might be why I need to use Fisher's exact test
#Error is in -Inf values that weren't removed

counts.gg.ratioAndMods <- merge(counts.gg.ratio,aggrColorDat,by="Gene",all.x =TRUE)
for(i in 1:(length(unique(counts.gg.ratioAndMods$AggrAllMod))-1)){
	temp <- subset(counts.gg.ratioAndMods,AggrAllMod %in% i)
	ggplot(temp,aes(x = Gene, y = AD.LogRatio, colour = Cult, shape = DPA)) +
	geom_point() +
	ggtitle(paste("Module",i,"Homoeolog Bias",sep = " ")) +
	geom_hline(yintercept = 0, colour = "yellow")
	ggsave(paste("Module",i,"LogHomoeoBias.pdf",sep = ""),device = "pdf")
}

counts.gg.ratioAndMods$Condition <-paste(counts.gg.ratioAndMods$Cult,counts.gg.ratioAndMods$DPA,sep=".")

#Try to redo without introducing divide by 0 errors
counts.AD.logratio <- log2(counts.A.raw)-log2(counts.D.raw)
rownames(counts.AD.logratio) <- gsub(pat= "_A",rep = "",rownames(counts.AD.logratio))
pdf("homoeoBoxplot.pdf")
boxplot(counts.AD.logratio)
dev.off()
#everything centered on zero = good!
counts.AD.finite.logratio <- counts.AD.logratio[apply(counts.AD.logratio,2,is.finite)]
apply(counts.AD.finite.logratio,2,quantile)

apply(counts.AD.logratio,2,function(x)quantile(x,c(0.05, 0.95)) )



n=24
stat<-data.frame(mean=rowMeans(counts.AD.logratio), sd=apply(counts.AD.logratio,1,sd) )
stat$t<-stat$mean/(stat$sd/sqrt(n))
stat$p <- 2*pt(-abs(stat$t),df=n-1)
stat$p.bh<-p.adjust(stat$p, method="BH")
stat$p.bonferroni<-p.adjust(stat$p, method="bonferroni")
# volcano plot
plot(stat$mean,log2(stat$p ))
abline(h=log2(0.05))  # pvalue
abline(v=2.8)
abline(v=-2.8)

#Forget this right now - difficulties with undefined values
sf <- sizeFactors(AD1.AD.dds0)#these are from the full library
#get count tables ready to combine for this test
rownames(counts.D.raw) <- gsub(pat = "_D",rep = "",rownames(counts.D.raw))
rownames(counts.A.raw) <- gsub(pat = "_A",rep = "",rownames(counts.A.raw))
colnames(counts.D.raw) <- paste(colnames(counts.D.raw),"_D",sep = "")
colnames(counts.A.raw) <- paste(colnames(counts.A.raw),"_A",sep = "")
counts.AD.bySample.raw <- cbind(counts.A.raw,counts.D.raw)
coldata.HomoeoBias <- data.frame(sample = colnames(counts.AD.bySample.raw),Subgenome = c(rep("A",24),rep("D",24)),Cult = c(rep("D",12),rep("W",12),rep("D",12),rep("W",12)),DPA = rep(c(5,10,15,20),times=12))
coldata.HomoeoBias$Sample.Raw <- gsub(pat = "_A|_D",rep = "",coldata.HomoeoBias$sample)
coldata.HomoeoBias <- cbind(coldata.HomoeoBias,rep(sf,2))
names(coldata.HomoeoBias)[7] <- "sizeFactors"
for( i in c("D","W") )
{
    select<-which(coldata.HomoeoBias$Cult==i)
    AD.homoeoBias.dds <- DESeqDataSetFromMatrix( countData = counts.AD.bySample.raw[,select], colData = coldata.HomoeoBias[select,], design = ~ Subgenome + DPA)
	AD.homoeoBias.dds <- AD.homoeoBias.dds[rowSums(counts(AD.homoeoBias.dds)) > 24,]
    # instead of estimating size factors from partitioned homoeologs reads, giving full library factors to both At and Dt bins
    sizeFactors(AD.homoeoBias.dds)<-coldata.HomoeoBias$sizeFactors[select]
    AD.homoeoBias.dds <- estimateDispersions(AD.homoeoBias.dds)
    AD.homoeoBias.dds <- nbinomWaldTest(AD.homoeoBias.dds)
    res <- results(AD.homoeoBias.dds, c("Subgenome", "A","D"))
    print(i)
    print(res)
    summary(res, alpha=0.05)
    assign(paste(i,"homoeobias",sep="."),res)
}
#[1] "D"
#log2 fold change (MAP): Subgenome A vs D
#Wald test p-value: Subgenome A vs D
#DataFrame with 26315 rows and 6 columns
#                   baseMean log2FoldChange     lfcSE        stat       pvalue
#                  <numeric>      <numeric> <numeric>   <numeric>    <numeric>
#Gorai.001G000300   1.844922    0.641872554 0.5054923  1.26979686 2.041570e-01
#Gorai.001G000400   7.623851    3.795013337 0.4198230  9.03955464 1.573114e-19
#Gorai.001G000500  30.241549    0.005638241 0.1876536  0.03004601 9.760304e-01
#Gorai.001G000800   6.421469    0.415908951 0.4028546  1.03240466 3.018826e-01
#Gorai.001G000900 239.457631   -0.301300521 0.1323294 -2.27689841 2.279229e-02
#...                     ...            ...       ...         ...          ...
#Gorai.013G272300 307.951516      -0.375612 0.1465125   -2.563687 1.035670e-02
#Gorai.013G272400   6.648970      -1.836069 0.3454621   -5.314819 1.067635e-07
#Gorai.013G272500  10.642212       1.982496 0.3246169    6.107187 1.014023e-09
#Gorai.013G272600   3.383006      -2.127926 0.4188912   -5.079901 3.776313e-07
#Gorai.013G272700  13.328595       4.908505 0.4001607   12.266334 1.373317e-34
#                         padj
#                    <numeric>
#Gorai.001G000300 2.739930e-01
#Gorai.001G000400 2.104727e-18
#Gorai.001G000500 9.821729e-01
#Gorai.001G000800 3.813847e-01
#Gorai.001G000900 3.963408e-02
#...                       ...
#Gorai.013G272300 1.949210e-02
#Gorai.013G272400 4.616855e-07
#Gorai.013G272500 5.643496e-09
#Gorai.013G272600 1.525214e-06
#Gorai.013G272700 4.740762e-33
#
#out of 26315 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)     : 7375, 28%
#LFC < 0 (down)   : 8059, 31%
#outliers [1]     : 252, 0.96%
#low counts [2]   : 0, 0%
#(mean count < 1)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results
#LFC > 0.5
#out of 16478 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)     : 6507, 39%
#LFC < 0 (down)   : 7106, 43%
#outliers [1]     : 213, 1.3%
#low counts [2]   : 0, 0%
#(mean count < 1)


#[1] "W"
#log2 fold change (MAP): Subgenome A vs D
#Wald test p-value: Subgenome A vs D
#DataFrame with 26242 rows and 6 columns
#                   baseMean log2FoldChange     lfcSE      stat       pvalue
#                  <numeric>      <numeric> <numeric> <numeric>    <numeric>
#Gorai.001G000300   1.918624      0.6650898 0.5078812  1.309538 1.903521e-01
#Gorai.001G000400   9.948167      3.1960036 0.4792307  6.669029 2.575016e-11
#Gorai.001G000500  23.304728      0.8702340 0.5595722  1.555177 1.199038e-01
#Gorai.001G000800   8.659669      1.1743274 0.3584406  3.276212 1.052097e-03
#Gorai.001G000900 235.282854     -0.1660611 0.1377794 -1.205268 2.280997e-01
#...                     ...            ...       ...       ...          ...
#Gorai.013G272300 298.978772     -0.4829393 0.1836809 -2.629229 8.557869e-03
#Gorai.013G272400  11.351113     -2.0967986 0.3170583 -6.613290 3.758712e-11
#Gorai.013G272500   7.272771      3.7705590 0.6060465  6.221567 4.922133e-10
#Gorai.013G272600   2.545222     -0.7021493 0.4173158 -1.682537 9.246476e-02
#Gorai.013G272700   8.075435      3.9967365 0.4573644  8.738627 2.359605e-18
#                         padj
#                    <numeric>
#Gorai.001G000300 2.746395e-01
#Gorai.001G000400 2.558154e-10
#Gorai.001G000500 1.864415e-01
#Gorai.001G000800 2.936853e-03
#Gorai.001G000900 3.177600e-01
#...                       ...
#Gorai.013G272300 1.904292e-02
#Gorai.013G272400 3.674647e-10
#Gorai.013G272500 4.153855e-09
#Gorai.013G272600 1.496214e-01
#Gorai.013G272700 4.883966e-17
# LFC > 0
#out of 26242 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)     : 6307, 24%
#LFC < 0 (down)   : 6960, 27%
#outliers [1]     : 452, 1.7%
#low counts [2]   : 0, 0%
#(mean count < 1)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results
# LFC > 0.5
#out of 15629 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)     : 5598, 36%
#LFC < 0 (down)   : 6230, 40%
#outliers [1]     : 344, 2.2%
#low counts [2]   : 0, 0%
#(mean count < 1)



for( i in unique(coldata.HomoeoBias$Sample.Inter) )
{
    select<-which(coldata.HomoeoBias$Sample.Inter==i)
    AD.homoeoBias.dds <- DESeqDataSetFromMatrix(countData = counts.AD.bySample.raw[,select], colData = coldata.HomoeoBias[select,], design = ~ Subgenome)
	AD.homoeoBias.dds <- AD.homoeoBias.dds[rowSums(counts(AD.homoeoBias.dds)) > 24,]
    # instead of estimating size factors from partitioned homoeologs reads, giving full library factors to both At and Dt bins
    sizeFactors(AD.homoeoBias.dds)<-coldata.HomoeoBias$sizeFactors[select]
    AD.homoeoBias.dds <- estimateDispersions(AD.homoeoBias.dds)
    AD.homoeoBias.dds <- nbinomWaldTest(AD.homoeoBias.dds)
    res <- results(AD.homoeoBias.dds)
    print(i)
    print(res)
    summary(res, alpha=0.05)
    assign(paste(i,"homoeobias",sep = "."),res)
}
#D.5.homoeobias
#log2 fold change (MAP): Subgenome D vs A
#Wald test p-value: Subgenome D vs A
#DataFrame with 22077 rows and 6 columns
#                   baseMean log2FoldChange     lfcSE       stat       pvalue
#                  <numeric>      <numeric> <numeric>  <numeric>    <numeric>
#Gorai.001G000400   7.141610     -3.6959883 0.9216119 -4.0103521 6.062827e-05
#Gorai.001G000500  38.783127      0.5193597 0.4064832  1.2776905 2.013586e-01
#Gorai.001G000800   4.218152      0.2243703 0.9436278  0.2377742 8.120563e-01
#Gorai.001G000900 318.233819      0.1776133 0.2487016  0.7141626 4.751267e-01
#Gorai.001G001000  29.286889      0.3499824 0.4496722  0.7783058 4.363888e-01
#...                     ...            ...       ...        ...          ...
#Gorai.013G272300 254.911426      0.5536029 0.2356029   2.349729 1.878707e-02
#Gorai.013G272400  10.123316      1.9857946 0.8093737   2.453495 1.414753e-02
#Gorai.013G272500  11.418509     -1.7791605 0.7510504  -2.368896 1.784126e-02
#Gorai.013G272600   4.744419      2.3915938 1.0408521   2.297727 2.157734e-02
#Gorai.013G272700  13.554742     -4.4685184 0.7544984  -5.922502 3.170792e-09
#                         padj
#                    <numeric>
#Gorai.001G000400  0.000330346
#Gorai.001G000500  0.342438810
#Gorai.001G000800  0.886557569
#Gorai.001G000900  0.629302134
#Gorai.001G001000  0.594512792
#...                       ...
#Gorai.013G272300 5.099316e-02
#Gorai.013G272400 4.011700e-02
#Gorai.013G272500 4.875090e-02
#Gorai.013G272600 5.723160e-02
#Gorai.013G272700 3.666082e-08
#LFC > 0
#out of 22077 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)     : 4264, 19%
#LFC < 0 (down)   : 3812, 17%
#outliers [1]     : 86, 0.39%
#low counts [2]   : 0, 0%
#(mean count < 2)
#LFC > 0.5
#out of 13660 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)     : 4261, 31%
#LFC < 0 (down)   : 3810, 28%
#outliers [1]     : 69, 0.51%
#low counts [2]   : 0, 0%
#(mean count < 2)


#> D.10.homoeobias
#log2 fold change (MAP): Subgenome D vs A
#Wald test p-value: Subgenome D vs A
#DataFrame with 22101 rows and 6 columns
#                   baseMean log2FoldChange     lfcSE       stat       pvalue
#                  <numeric>      <numeric> <numeric>  <numeric>    <numeric>
#Gorai.001G000300   3.420375     -0.2737666 0.9597397 -0.2852509 7.754519e-01
#Gorai.001G000400   9.787293     -3.3791201 0.8047340 -4.1990521 2.680347e-05
#Gorai.001G000500  34.088044     -0.1491153 0.4989583 -0.2988533 7.650520e-01
#Gorai.001G000800   5.811792     -0.2587905 0.8343108 -0.3101847 7.564205e-01
#Gorai.001G000900 269.065729      0.3787379 0.2301116  1.6458876 9.978690e-02
#...                     ...            ...       ...        ...          ...
#Gorai.013G272300 257.447948      0.3101447 0.2228413   1.391774 1.639909e-01
#Gorai.013G272400   7.208552      1.2210893 0.8130879   1.501792 1.331507e-01
#Gorai.013G272500  14.269342     -1.8972945 0.6327628  -2.998429 2.713751e-03
#Gorai.013G272600   4.671027      2.1321377 1.0010571   2.129886 3.318101e-02
#Gorai.013G272700  13.260367     -4.4125332 0.8626130  -5.115310 3.132260e-07
#                         padj
#                    <numeric>
#Gorai.001G000300 0.8747527134
#Gorai.001G000400 0.0001958352
#Gorai.001G000500 0.8681304001
#Gorai.001G000800 0.8624768605
#Gorai.001G000900 0.2109829194
#...                       ...
#Gorai.013G272300 3.056185e-01
#Gorai.013G272400 2.621157e-01
#Gorai.013G272500 1.098979e-02
#Gorai.013G272600 8.857121e-02
#Gorai.013G272700 3.537895e-06

#LFC > 0
#out of 22101 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)     : 3860, 17%
#LFC < 0 (down)   : 3425, 15%
#outliers [1]     : 87, 0.39%
#low counts [2]   : 0, 0%
#(mean count < 3)
#LFC > 0.5
#out of 13314 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)     : 3858, 29%
#LFC < 0 (down)   : 3422, 26%
#outliers [1]     : 77, 0.58%
#low counts [2]   : 0, 0%
#(mean count < 3)


#> D.15.homoeobias
#log2 fold change (MAP): Subgenome D vs A
#Wald test p-value: Subgenome D vs A
#DataFrame with 22889 rows and 6 columns
#                   baseMean log2FoldChange     lfcSE       stat       pvalue
#                  <numeric>      <numeric> <numeric>  <numeric>    <numeric>
#Gorai.001G000400   7.296319     -2.8600467 0.8877562 -3.2216579  0.001274512
#Gorai.001G000500  24.998707     -0.3814671 0.4659312 -0.8187200  0.412946191
#Gorai.001G000800   6.677462     -1.2076799 0.8039011 -1.5022741  0.133026335
#Gorai.001G000900 260.766638      0.3338812 0.2237820  1.4919932  0.135700927
#Gorai.001G001000  34.668948      0.2090346 0.4171439  0.5011091  0.616294366
#...                     ...            ...       ...        ...          ...
#Gorai.013G271700 908.519559     -0.2299640 0.9117293 -0.2522283 8.008646e-01
#Gorai.013G272300 237.783366      0.2242715 0.2627491  0.8535576 3.933502e-01
#Gorai.013G272400   6.407346      1.5439499 0.7857463  1.9649471 4.942036e-02
#Gorai.013G272500  11.187159     -1.7446598 0.6429466 -2.7135377 6.656900e-03
#Gorai.013G272700  16.875239     -4.5176203 0.7583806 -5.9569300 2.570203e-09
#                         padj
#                    <numeric>
#Gorai.001G000400  0.005653166
#Gorai.001G000500  0.589076049
#Gorai.001G000800  0.261315248
#Gorai.001G000900  0.264765571
#Gorai.001G001000  0.759187769
#...                       ...
#Gorai.013G271700 8.866882e-01
#Gorai.013G272300 5.699962e-01
#Gorai.013G272400 1.217656e-01
#Gorai.013G272500 2.301925e-02
#Gorai.013G272700 4.252127e-08
#LFC > 0
#out of 22889 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)     : 3999, 17%
#LFC < 0 (down)   : 3590, 16%
#outliers [1]     : 108, 0.47%
#low counts [2]   : 0, 0%
#(mean count < 3)
#LFC > 0.5
#out of 13609 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)     : 3996, 29%
#LFC < 0 (down)   : 3588, 26%
#outliers [1]     : 91, 0.67%
#low counts [2]   : 0, 0%
#(mean count < 3)


#> D.20.homoeobias
#log2 fold change (MAP): Subgenome D vs A
#Wald test p-value: Subgenome D vs A
#DataFrame with 20401 rows and 6 columns
#                   baseMean log2FoldChange     lfcSE        stat       pvalue
#                  <numeric>      <numeric> <numeric>   <numeric>    <numeric>
#Gorai.001G000400   6.270183    -3.28777388 1.2572961 -2.61495596  0.008923896
#Gorai.001G000500  23.096318    -0.03818103 0.6002158 -0.06361217  0.949279043
#Gorai.001G000900 109.764336     0.33012711 0.5189809  0.63610649  0.524707016
#Gorai.001G001000  17.609325     0.07387180 0.6673755  0.11069000  0.911862178
#Gorai.001G001100 103.449203     0.51076435 0.4859948  1.05096680  0.293273839
#...                     ...            ...       ...         ...          ...
#Gorai.013G271600 165.062298     -0.8998704 0.3190944  -2.8200756 0.0048012333
#Gorai.013G271700 218.607224      0.1863904 0.9689418   0.1923649 0.8474563588
#Gorai.013G272300 481.663324      0.3380781 0.2732648   1.2371812 0.2160198199
#Gorai.013G272500   5.693839     -1.9549243 1.2050557  -1.6222689 0.1047457933
#Gorai.013G272700   9.624032     -4.5739961 1.1347995  -4.0306645 0.0000556194
#                         padj
#                    <numeric>
#Gorai.001G000400   0.03657706
#Gorai.001G000500   0.97832138
#Gorai.001G000900   0.72153115
#Gorai.001G001000   0.96065526
#Gorai.001G001100   0.50831625
#...                       ...
#Gorai.013G271600 0.0219190025
#Gorai.013G271700 0.9259964169
#Gorai.013G272300 0.4156994173
#Gorai.013G272500 0.2474066994
#Gorai.013G272700 0.0004855384
#LFC > 0
#out of 20401 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)     : 2783, 14%
#LFC < 0 (down)   : 2580, 13%
#outliers [1]     : 122, 0.6%
#low counts [2]   : 0, 0%
#(mean count < 3)
#LFC > 0.5
#out of 12687 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)     : 2783, 22%
#LFC < 0 (down)   : 2580, 20%
#outliers [1]     : 112, 0.88%
#low counts [2]   : 0, 0%
#(mean count < 3)


#> W.5.homoeobias
#log2 fold change (MAP): Subgenome D vs A
#Wald test p-value: Subgenome D vs A
#DataFrame with 20558 rows and 6 columns
#                   baseMean log2FoldChange     lfcSE       stat       pvalue
#                  <numeric>      <numeric> <numeric>  <numeric>    <numeric>
#Gorai.001G000500  29.931996     -0.5936749 0.9186304 -0.6462609    0.5181104
#Gorai.001G000800   5.296194     -0.7531966 0.8675312 -0.8682069    0.3852811
#Gorai.001G000900 277.612139      0.3114601 0.2352856  1.3237533    0.1855850
#Gorai.001G001000  32.483912     -0.1167765 0.5388356 -0.2167201    0.8284265
#Gorai.001G001100 229.412205     -0.3174040 0.4161518 -0.7627120    0.4456352
#...                     ...            ...       ...        ...          ...
#Gorai.013G272300 243.758622      0.4299627 0.3507694  1.2257701 2.202852e-01
#Gorai.013G272400  13.499646      2.3549525 0.7303618  3.2243641 1.262528e-03
#Gorai.013G272500   8.449244     -3.6797153 1.0566111 -3.4825635 4.966374e-04
#Gorai.013G272600   4.511320      0.1803387 0.8516347  0.2117559 8.322975e-01
#Gorai.013G272700   9.491514     -4.0679648 0.9264864 -4.3907442 1.129634e-05
#                         padj
#                    <numeric>
#Gorai.001G000500    0.6895649
#Gorai.001G000800    0.5695568
#Gorai.001G000900    0.3480033
#Gorai.001G001000    0.9053698
#Gorai.001G001100    0.6249907
#...                       ...
#Gorai.013G272300 0.3921022540
#Gorai.013G272400 0.0063536956
#Gorai.013G272500 0.0028630458
#Gorai.013G272600 0.9077063665
#Gorai.013G272700 0.0001076763
#LFC > 0
#out of 20558 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)     : 3145, 15%
#LFC < 0 (down)   : 2785, 14%
#outliers [1]     : 131, 0.64%
#low counts [2]   : 0, 0%
#(mean count < 4)
#LFC > 0.5
#out of 12115 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)     : 3143, 26%
#LFC < 0 (down)   : 2784, 23%
#outliers [1]     : 118, 0.97%
#low counts [2]   : 0, 0%
#(mean count < 4)


#> W.10.homoeobias
#log2 fold change (MAP): Subgenome D vs A
#Wald test p-value: Subgenome D vs A
#DataFrame with 21631 rows and 6 columns
#                    baseMean log2FoldChange     lfcSE         stat       pvalue
#                   <numeric>      <numeric> <numeric>    <numeric>    <numeric>
#Gorai.001G000400    8.130867   -3.406878437 0.9589863 -3.552582960 0.0003814687
#Gorai.001G000500   24.813797   -0.761597746 0.9048301 -0.841702474 0.3999545132
#Gorai.001G000800    9.473208   -0.197938624 0.7843461 -0.252361330 0.8007617876
#Gorai.001G000900  252.734604    0.186032587 0.2231112  0.833810962 0.4043875169
#Gorai.001G001000   46.574460   -0.001946323 0.3571737 -0.005449233 0.9956521630
#...                      ...            ...       ...          ...          ...
#Gorai.013G271700 2284.678222     -0.8819219 0.2898904    -3.042260 0.0023480932
#Gorai.013G272300  231.987900      0.4758979 0.4164208     1.142829 0.2531094770
#Gorai.013G272400   12.288775      2.3909263 0.7679364     3.113443 0.0018491816
#Gorai.013G272500   12.064450     -2.6736145 1.0426601    -2.564224 0.0103406679
#Gorai.013G272700    7.967918     -3.3431059 1.0017422    -3.337292 0.0008459913
#                        padj
#                   <numeric>
#Gorai.001G000400 0.003232319
#Gorai.001G000500 0.633101535
#Gorai.001G000800 0.907191884
#Gorai.001G000900 0.637626306
#Gorai.001G001000 0.997883837
#...                      ...
#Gorai.013G271700 0.014374536
#Gorai.013G272300 0.482716252
#Gorai.013G272400 0.011850996
#Gorai.013G272500 0.048427178
#Gorai.013G272700 0.006263371

#LFC > 0
#out of 21631 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)     : 2485, 11%
#LFC < 0 (down)   : 2141, 9.9%
#outliers [1]     : 168, 0.78%
#low counts [2]   : 0, 0%
#(mean count < 4)
#LFC > 0.5
#out of 11986 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)     : 2485, 21%
#LFC < 0 (down)   : 2140, 18%
#outliers [1]     : 150, 1.3%
#low counts [2]   : 0, 0%
#(mean count < 4)

#> W.15.homoeobias
#log2 fold change (MAP): Subgenome D vs A
#Wald test p-value: Subgenome D vs A
#DataFrame with 23199 rows and 6 columns
#                   baseMean log2FoldChange     lfcSE       stat       pvalue
#                  <numeric>      <numeric> <numeric>  <numeric>    <numeric>
#Gorai.001G000400   7.894739     -3.0416480 0.9945668 -3.0582642  0.002226232
#Gorai.001G000500  22.127526     -1.2738155 0.9280735 -1.3725373  0.169896239
#Gorai.001G000800   9.531078     -1.5509390 0.7024015 -2.2080520  0.027240646
#Gorai.001G000900 244.538436      0.1294900 0.4028635  0.3214241  0.747889021
#Gorai.001G001000  36.996390     -0.7543667 0.3951012 -1.9092997  0.056223439
#...                     ...            ...       ...        ...          ...
#Gorai.013G272100   3.204374     -2.4089486 1.2224465  -1.970596 4.877006e-02
#Gorai.013G272300 271.750345      0.5967948 0.2590482   2.303799 2.123395e-02
#Gorai.013G272400  12.350435      2.3290510 0.5937627   3.922528 8.762468e-05
#Gorai.013G272500   4.297953     -3.1461478 1.0328791  -3.045998 2.319092e-03
#Gorai.013G272700   9.168184     -3.0167287 0.7226497  -4.174538 2.985909e-05
#                         padj
#                    <numeric>
#Gorai.001G000400   0.01560149
#Gorai.001G000500   0.39299651
#Gorai.001G000800   0.10982867
#Gorai.001G000900   0.88441794
#Gorai.001G001000   0.18495882
#...                       ...
#Gorai.013G272100 0.1668472744
#Gorai.013G272300 0.0913204519
#Gorai.013G272400 0.0010827395
#Gorai.013G272500 0.0160736087
#Gorai.013G272700 0.0004347475

#LFC > 0
#out of 23199 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)     : 2379, 10%
#LFC < 0 (down)   : 1990, 8.6%
#outliers [1]     : 500, 2.2%
#low counts [2]   : 0, 0%
#(mean count < 2)
#LFC > 0.5
#out of 13008 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)     : 2379, 18%
#LFC < 0 (down)   : 1990, 15%
#outliers [1]     : 368, 2.8%
#low counts [2]   : 0, 0%
#(mean count < 2)


#> W.20.homoeobias
#log2 fold change (MAP): Subgenome D vs A
#Wald test p-value: Subgenome D vs A
#DataFrame with 21015 rows and 6 columns
#                   baseMean log2FoldChange     lfcSE        stat       pvalue
#                  <numeric>      <numeric> <numeric>   <numeric>    <numeric>
#Gorai.001G000400   19.97107    -2.15630752 0.7773084 -2.77406934  0.005535988
#Gorai.001G000500   16.34559    -0.40135530 1.0474148 -0.38318658  0.701581428
#Gorai.001G000800   10.33820    -2.07545622 0.8887892 -2.33515015  0.019535578
#Gorai.001G000900  166.24624     0.01697051 0.3917583  0.04331883  0.965447378
#Gorai.001G001000   38.18488    -0.76708105 0.4515895 -1.69862450  0.089389957
#...                     ...            ...       ...         ...          ...
#Gorai.013G271600 106.299793     -0.8720694 0.4247577  -2.0530985 0.0400630279
#Gorai.013G271700 353.699734     -0.1464523 0.7908602  -0.1851810 0.8530870657
#Gorai.013G272300 448.418221      0.4219206 0.4354513   0.9689270 0.3325815884
#Gorai.013G272400   7.265597      0.2700389 0.7914951   0.3411757 0.7329713041
#Gorai.013G272700   5.674126     -3.9646616 1.1480008  -3.4535354 0.0005532897
#                       padj
#                  <numeric>
#Gorai.001G000400 0.02798159
#Gorai.001G000500 0.84861921
#Gorai.001G000800 0.07605648
#Gorai.001G000900 0.98657928
#Gorai.001G001000 0.23560940
#...                     ...
#Gorai.013G271600 0.13098039
#Gorai.013G271700 0.93110182
#Gorai.013G272300 0.56127881
#Gorai.013G272400 0.86793827
#Gorai.013G272700 0.00419081
#LFC > 0
#adjusted p-value < 0.05
#LFC > 0 (up)     : 2498, 12%
#LFC < 0 (down)   : 2209, 11%
#outliers [1]     : 284, 1.4%
#low counts [2]   : 0, 0%
#(mean count < 3)
#LFC > 0.5
#out of 12388 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)     : 2498, 20%
#LFC < 0 (down)   : 2207, 18%
#outliers [1]     : 244, 2%
#low counts [2]   : 0, 0%
#(mean count < 3)

AD.homoeoBias.dds <- DESeqDataSetFromMatrix(countData = counts.AD.bySample.raw, colData = coldata.HomoeoBias, design = ~ Subgenome)
AD.homoeoBias.dds <- AD.homoeoBias.dds[rowSums(counts(AD.homoeoBias.dds)) > 24,]
sizeFactors(AD.homoeoBias.dds)<-coldata.HomoeoBias$sizeFactors
AD.homoeoBias.dds <- estimateDispersions(AD.homoeoBias.dds)
AD.homoeoBias.dds <- nbinomWaldTest(AD.homoeoBias.dds)
res <- results(AD.homoeoBias.dds)
for( i in sort(as.numeric(unique(aggrColorDat$AggrAllMod))) )
{
    select<-which(aggrColorDat$AggrAllMod==i)
	selectGenesFromAggr <- aggrColorDat$Gene[select]
	selectGenes <- rownames(res) %in% selectGenesFromAggr
	temp.res <- res[selectGenes & !is.na(res$padj),]
	up.weak <- dim(temp.res[temp.res$padj < 0.05 & temp.res$log2FoldChange > 0,])[1]
	down.weak <- dim(temp.res[temp.res$padj < 0.05 & temp.res$log2FoldChange < 0,])[1]
	up.strong <- dim(temp.res[temp.res$padj < 0.05 & temp.res$log2FoldChange > 0.5,])[1]
	down.strong <- dim(temp.res[temp.res$padj < 0.05 & temp.res$log2FoldChange < -0.5,])[1]
    print(paste("Module",i,down.weak,down.strong,up.weak,up.strong,sep = " "))
}



#Will need to redo to get results if I want those!

