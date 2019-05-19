##R script for differential expression analysis of G. hirsutum (AD1) fiber expression data
library(DESeq2)

#read in pre-generated files with read counts for each gene
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

#reorganize so that samples are in correct order (5,10,15,20 dpa)
counts.A.raw <- counts.A.raw[,c(4,1:3,8,5:7,12,9:11,17,14:16,21,18:20,25,13,23,24)]
counts.D.raw <- counts.D.raw[,c(4,1:3,8,5:7,12,9:11,17,14:16,21,18:20,25,13,23,24)]

#rename genes so that they are A or D copy
A.names <- paste(rownames(counts.A.raw),"A",sep="_")
D.names <- paste(rownames(counts.D.raw),"D",sep="_")

rownames(counts.A.raw) <- A.names
rownames(counts.D.raw) <- D.names

#this creates the final raw read count table with A and D read counts together
counts.AD.raw <- rbind(counts.A.raw,counts.D.raw)

#make colMat - metadata file
#Note: Yuc(atanense) == TX2094
accession <- c(rep("AD1.CRB252",4),rep("AD1.Maxxa",4),rep("AD1.TM1",4),rep("AD1.TX2095",4),
               rep("AD1.TX665",4),rep("AD1.Yuc",4))
timepoint <- rep(c(5,10,15,20),6)
cult <- c(rep("D",12),rep("W",12))
colMat <- data.frame(accession = accession, timepoint = factor(timepoint), cult = factor(cult))
colMat$sample0 <- paste(colMat$cult,colMat$timepoint,sep = ".")

#Build DESeq data table
#used average of one read per gene per sample for consistency with coexpression network
AD1.part.dds <- DESeqDataSetFromMatrix(countData = counts.AD.raw, colData = colMat, design = ~ cult + timepoint + cult:timepoint)
AD1.part.dds <- AD1.part.dds[rowSums(counts(AD1.part.dds))/24 >= 1,]
AD1.part.dds <- estimateSizeFactors(AD1.part.dds)

#Build heatmap and PCA of expression data
AD.part.rlog <- rlog(AD1.part.dds, blind=TRUE)
pdf("AD1.part.blind.rlog.pdf")
plot(assay(AD.part.rlog)[,1:2],pch=16, cex=0.3)
dev.off()
pdf("AD1.part.blind.PCA.pdf")
plotPCA(AD.part.rlog, intgroup = c("cult", "timepoint"), ntop = 0)
dev.off()

AD1.part.dds <- DESeq(AD1.part.dds)

#Because of the interaction effect, there is a bit of art to writing contrasts
#see http://rpubs.com/ge600/deseq2 for help devising contrasts

#get resultsNames for contrasts. i.e pairwise comparisons
resultsNames(AD1.part.dds)
#"Intercept"         "cult_W_vs_D"       "timepoint_10_vs_5" "timepoint_15_vs_5"
#"timepoint_20_vs_5" "cultW.timepoint10" "cultW.timepoint15" "cultW.timepoint20"

#Wild vs dom at each timepoint
#5 dpa
results.wvd.05dpa <- results(AD1.part.dds,contrast=c("cult","W","D"))
#10 dpa
results.wvd.10dpa <- results(AD1.part.dds,contrast=list(c("cult_W_vs_D","cultW.timepoint10")))
#15 dpa
results.wvd.15dpa <- results(AD1.part.dds,contrast=list(c("cult_W_vs_D","cultW.timepoint15")))
#20 dpa
results.wvd.20dpa <- results(AD1.part.dds,contrast=list(c("cult_W_vs_D","cultW.timepoint20")))
#Then 5 dpa vs 10 dpa vs 15 dpa vs 20 dpa
#5dpa vs 10dpa in dom
results.dom.05v10 <- results(AD1.part.dds,contrast=c("timepoint","10","5"))
#10dpa vs 15dpa in dom
results.dom.10v15 <- results(AD1.part.dds,contrast=c("timepoint","15","10"))
#could be written "results(AD1.part.dds,contrast=c(0,0,-1,1,0,0,0,0))"
#15dpa vs 20dpa in dom
results.dom.15v20 <- results(AD1.part.dds,contrast=c("timepoint","20","15"))
#5dpa vs 10dpa in wild
results.wild.05v10 <- results(AD1.part.dds,contrast=list(c("timepoint_10_vs_5","cultW.timepoint10")))
#10dpa vs 15dpa in wild
results.wild.10v15 <- results(AD1.part.dds,contrast=c(0,0,-1,1,0,-1,1,0))
#15dpa vs 20dpa in wild
results.wild.15v20 <- results(AD1.part.dds,contrast=c(0,0,0,-1,1,0,-1,1))
#Interaction effects (difference in effect of wild v dom at each time point difference)
#Wild v Dom at 10 dpa vs 5 dpa
results.wvd.05v10dpa <- results(AD1.part.dds,name="cultW.timepoint10")
#Wild v Dom at 15 vs 10 dpa
results.wvd.10v15dpa <- results(AD1.part.dds,contrast=c(0,0,0,0,0,-1,1,0))
#Wild v Dom at 20 dpa vs 15 dpa
results.wvd.15v20dpa <- results(AD1.part.dds,contrast=c(0,0,0,0,0,0,-1,1))

#write output tables
#wild vs dom DEG lists
write.table(results.wvd.05dpa,file="DEG.WvD.05dpa.txt")
write.table(results.wvd.10dpa,file="DEG.WvD.10dpa.txt")
write.table(results.wvd.15dpa,file="DEG.WvD.15dpa.txt")
write.table(results.wvd.20dpa,file="DEG.WvD.20dpa.txt")
#consecutive DPA in wild and dom
write.table(results.dom.05v10,file="DEG.D.05v10dpa.txt")
write.table(results.dom.10v15,file="DEG.D.10v15dpa.txt")
write.table(results.dom.15v20,file="DEG.D.15v20dpa.txt")
write.table(results.wild.05v10,file="DEG.W.05v10dpa.txt")
write.table(results.wild.10v15,file="DEG.W.10v15dpa.txt")
write.table(results.wild.15v20,file="DEG.W.15v20dpa.txt")
#consecutive interaction terms
write.table(results.wvd.05v10dpa,file="DEG.WvD.05v10dpa.txt")
write.table(results.wvd.10v15dpa,file="DEG.WvD.10v15dpa.txt")
write.table(results.wvd.15v20dpa,file="DEG.WvD.15v20dpa.txt")

#Note: all cult comparisons are W - D
#All timepoint comparisons are the higher number - the lower (e.g. 10 - 5)
#So when you you see the LFC up term, that means those are genes upregulated in
#10 dpa as compared to 5 dpa

##GO analysis of DE genes with topGO
library(topGO)

#make objects with DE results
deg.res = ls(pattern = 'results')
for(j in 1:length(deg.res)){
	deg.set <- get(deg.res[j])
	assign(paste(deg.res[j],".up",sep=""),deg.set[deg.set$log2FoldChange > 0,])
	assign(paste(deg.res[j],".down",sep=""),deg.set[deg.set$log2FoldChange < 0,])	
}

#set up ID and GO term mapping
geneID2GO <- readMappings(file = "D5_genes2GO.AD.map", sep = " ", IDsep =",")
geneNames <- names(geneID2GO)
GOresults<-data.frame()

#set up table with different types of DGE terms
deg.full.res <- c(paste(deg.res,".up",sep=""),paste(deg.res,".down",sep=""))

#perform GO analysis over all types of DGE
for(m in 1:length(deg.full.res)){
	deg.set <- get(deg.full.res[m])
	sig.list <- rownames(deg.set[!is.na(deg.set$padj) & deg.set$padj < 0.05,])
	geneList <- factor(as.integer(geneNames %in% sig.list))
	if(length(levels(geneList)) != 2) next
	names(geneList) <- geneNames
	pdf(file=paste("GoTerms/",deg.full.res[m],".pdf", sep=""))
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
	#add enriched terms to whole GO enrichment results table
	if(dim(enrich)[1]>0)
	{
		enrichDE <- enrich
		enrichDE$DE <- deg.full.res[m]
		GOresults<-rbind(GOresults,enrichDE)
	}
}

#write out table of GO results
write.table(GOresults, file="GoTerms/GOresults.txt", sep="\t", row.names=FALSE)

#LRT to test effect of whole terms (development, domestication, interaction)

##LRT test for domestication status and timepoint terms
design(AD1.part.dds) <- ~cult + timepoint
AD1.part.domLRT.dds <- DESeq(AD1.part.dds, test="LRT", reduced = ~ timepoint)
AD1.part.timeLRT.dds <- DESeq(AD1.part.dds, test="LRT", reduced = ~ cult)

timeDEG.tbl <- results(AD1.part.timeLRT.dds)
domDEG.tbl <- results(AD1.part.domLRT.dds)
timeDEG.sig <- timeDEG.tbl$padj < 0.05
is.na(timeDEG.sig) <- "FALSE"
domDEG.sig <- domDEG.tbl$padj < 0.05
is.na(domDEG.sig) <- "FALSE"
part.timeDEG.list <- rownames(timeDEG.tbl)[timeDEG.sig]
part.domDEG.list <- rownames(domDEG.tbl)[domDEG.sig]
part.timeDEG.noNA.list<-part.timeDEG.list[!is.na(part.timeDEG.list)]
part.domDEG.noNA.list<-part.domDEG.list[!is.na(part.domDEG.list)]

write.table(file = "partitioned_time_DEG.txt",as.data.frame(part.timeDEG.noNA.list),quote=F)
write.table(file = "partitioned_dom_DEG.txt",as.data.frame(part.domDEG.noNA.list),quote=F)

#repeat GO enrichment analysis for 
remove(enrich)
GOresults <- data.frame()

go.list <- c("part.domDEG.noNA.list","part.timeDEG.noNA.list")
for(m in 1:length(go.list)){
	deg.set <- get(go.list[m])
	geneList <- factor(as.integer(geneNames %in% deg.set))
	if(length(levels(geneList)) != 2) next
	names(geneList) <- geneNames
	
	pdf(file=paste("GoTerms/",go.list[m],".pdf", sep=""))
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
		enrichDE$DE <- go.list[m]
		GOresults<-rbind(GOresults,enrichDE)
	}
}

write.table(GOresults, file="GoTerms/GOresults.part_LRT_dom_dev.txt", sep="\t", row.names=FALSE)

#Interaction effect
AD1.part.LRT.dds <- DESeq(AD1.part.dds, test="LRT", reduced = ~ cult + timepoint)
LRT.res <- results(AD1.part.LRT.dds)
write.table(LRT.res,file="DEG.LRTinteraction.txt")

#do enrichment for LRT interaction genes
LRT.res.up <- LRT.res[LRT.res$log2FoldChange > 0,]
LRT.res.down <- LRT.res[LRT.res$log2FoldChange < 0,]
deg.LRT.res <- c("LRT.res.up","LRT.res.down")

for(m in 1:length(deg.LRT.res)){
	deg.set <- get(deg.LRT.res[m])
	sig.list <- rownames(deg.set[!is.na(deg.set$padj) & deg.set$padj < 0.05,])
	geneList <- factor(as.integer(geneNames %in% sig.list))
	if(length(levels(geneList)) != 2) next
	names(geneList) <- geneNames
	
	pdf(file=paste("GoTerms/",deg.LRT.res[m],".pdf", sep=""))
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
		enrichDE$DE <- deg.full.res[m]
        	GOresults<-rbind(GOresults,enrichDE)
	}
}

write.table(GOresults, file="GoTerms/GOresults.wLRT.txt", sep="\t", row.names=FALSE)
