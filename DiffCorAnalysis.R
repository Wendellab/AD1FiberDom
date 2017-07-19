######
##DiffCor analysis for AD1 fiber
######


library(flashClust);
library(RColorBrewer);
library(ggplot2);
library(WGCNA);
options(stringsAsFactors = FALSE);

load("R-01-dataInput.RData")

for(i in 1:2)
{
    gg <- shortLabels[i]
    print(gg)
    subDat    <-  multiExpr[[i]]$data
    subDat   <-  apply(subDat,2,as.numeric)
    adj = adjacency(subDat, power = 12, type = "signed")
    assign(paste(gg,"adj",sep=""),adj)
}

# percentage of rewired edges R/C
# C = n(n-1)/2
# R = sum(abs(Xadj-Yadj))/2

n = checkSets(multiExpr)$nGenes #29706 in Aggr, 50996 in Part
C = n*(n-1)/2 #441208365 in Aggr, 1,300,270,510 in part
rate = sum(abs(Domadj-Wildadj))/2/C
print(rate)
# [1] 0.02949936 in Aggr, [1] 0.02785087 in Part

save(Wildadj,Domadj,n,C,rate,file = "rewiring.Dom_Wild.RData")#added *part.RData for part
#actually did not save partitioned, too big!
library(DiffCorr)


X<-t(multiExpr[[1]]$data)
Y<-t(multiExpr[[2]]$data)
#for large files, split into separate comps
X1 <-  X[1:25498,]
X2 <- X[25499:50996,]
Y1 <- Y[1:25498,]
Y2 <- Y[25499:50996,]
    
outfile <- paste0(shortLabels[1],"vs",shortLabels[2],".res.txt" )
    
comp.2.cc.fdr(output.file=outfile, X, Y, threshold=0.05)
##So it turns out this is too big of a matrix to execute the program
#Instead, we will take the list of genes called as DC and look at them specifically with DGCA

#We have this table now, now we need to interpret the data present therein
#Let's look at DiffCorr in consensus modules

lnames <- load("s8.moduleConsensus.Dom_Wild.RData")
lnames
# [1] "consMEs"               "unmergedLabels"        "moduleColors"
# [4] "moduleLabels"          "moduleLabels.adjusted" "consTree"
probes = colnames(multiExpr[[3]]$data) 
names(moduleLabels.adjusted) <- probes
mSizes<-as.data.frame(table(moduleLabels.adjusted))



x<-read.table("DomvsWild.res.txt",header=TRUE, sep="\t")
p= 2*nrow(x)/nGenes/(nGenes-1)

x<-x[,1:2]
x$X<-moduleLabels.adjusted[x$molecule.X]
x$Y<-moduleLabels.adjusted[x$molecule.Y]
# get rewired connnections considering consensus modules
a <- as.matrix(xtabs(~X+Y, data=x))

# sum up both directions
new <- a
new[lower.tri(new)]<-NA
new[upper.tri(new)] <- a[upper.tri(a)] + t(a)[upper.tri(a)]
# make result table
n<- as.data.frame(new)
n<-n[!is.na(n$Freq),]
levels(n$X) <- levels(n$Y)
n$type <- ifelse(n$X==n$Y,"within","between")

n$X.size <- mSizes$Freq[match(n$X, mSizes$moduleLabels.adjusted)]
n$Y.size <- mSizes$Freq[match(n$Y, mSizes$moduleLabels.adjusted)]
n$edge.total <- n$X.size*n$Y.size
n$edge.total[n$type=="within"] <- (n$X.size[n$type=="within"] -1) * n$X.size[n$type=="within"] /2
n$rewired.percentage<-n$Freq/n$edge.total

#
n$p = p
# Probability of observed number >= n$Freq based on the overall p value, if <0.05, significantly high rewiring
n$binom.pvalue <- pbinom(n$Freq-1, size=n$edge.total, prob=p, lower.tail=FALSE)
n$qvalue <- p.adjust(n$binom.pvalue, "BH")

outfile <- "DomvsWild.rewiringByCMod.txt"
write.table(n, outfile, row.names=FALSE, sep="\t" )

#Next, specific DCGs
print(paste0("Differential coexpression pairs: ", nrow(x)," (", p,")" ) )
#[1] "Differential coexpression pairs: 15364 (3.48225492052944e-05)"

ks <- as.data.frame(table(c(as.character(x$molecule.X), as.character(x$molecule.Y))) )
ks$P <- apply(ks, 1, function(x) pbinom(as.numeric(x[2]), size=nGenes, prob=p, lower.tail=FALSE) )
ks$q.bh<-p.adjust(ks$P, "BH")
ks.sig<-ks[ks$P<0.05 & ks$q.bh<0.05,]
n.dcg <- nrow(ks.sig)
p.dcg <- n.dcg/nGenes
print(paste0("Differential coexpression geness: ", n.dcg," (", p.dcg,")" ) )
#[1] "Differential coexpression geness: 1273 (0.0428532956305124)"

write.table(ks.sig, "DomvsWild.gene.dcg.txt", row.names=FALSE, sep="\t" )

# move all files into the "DC" folder
library(DESeq2)
lnames <- load("aggregateDEresults.RData")
#[1] "DE.aggr.devAll"   "DE.aggr.devByCon"

degByCon <- list()
for(i in 1:length(DE.aggr.devByCon)){
	check <- DE.aggr.devByCon[[i]]$padj < 0.05
	check[is.na(check)] = FALSE
	degs <- rownames(DE.aggr.devByCon[[i]][check,])
	degByCon[[i]] <- degs
}
comp.df<-as.data.frame( sapply(degByCon,length) )

dcgWildvDom <- ks.sig[,1]
sapply(dcgWildvDom,length)

# compare DE with DC
names(comp.df)<-"DE_No"
comp.df$DE_perc <- comp.df$DE_No/29706

for(i in 1:length(degByCon)){
	comp.df[i,"overlap"]<-length(intersect(degByCon[[i]], dcgWildvDom))
}

comp.df$overlap_perc<-comp.df$overlap/length(dcgWildvDom)
#comp.df
#   DE_No    DE_perc overlap overlap_perc
#1   1150 0.03871272       9  0.007069914
#2    971 0.03268700      30  0.023566379
#3   4316 0.14529051     139  0.109190888
#4   1279 0.04305528      16  0.012568735
#5   1049 0.03531273      24  0.018853103
#6   5667 0.19076954     230  0.180675570
#7    486 0.01636033      69  0.054202671
#8    630 0.02120784      46  0.036135114
#9    492 0.01656231      44  0.034564022
#10  1251 0.04211270     126  0.098978790

save(comp.df,dcgWildvDom,degByCon,ks,file = "DiffCorr.WildvsDom.DEandDC.RData")

degAllCon <- unique(sort(unlist(degByCon)))
degXdcg <-intersect(dcgWildvDom,degAllCon)
length(degXdcg) #407



#save and load DE/DC overlap and DC only list to PARTITIONED session 
#test DC genes first
testdcg.A <- paste0(dcg.aggr,"_A")
testdcg.A <- paste0(dcg.aggr,"_D")
testdcg.A <- paste0(dcg.aggr,"_A")
testdcg.D <- paste0(dcg.aggr,"_D")
testdcg <- c(testdcg.A,testdcg.D)
length(testdcg)
#[1] 2546

design_mat_AD <- cbind(Wild = c(rep(0,times=12),rep(1,times=12)),Dom = c(rep(1,times=12),rep(0,times=12)))
#test if DC genes from aggr are also DC in partitioned
ddCor_res <- ddcorAll(inputMat = inputDat,splitSet = testdcg, design = as.matrix(design_mat_AD),compare = c("Wild","Dom"),adjust = "fndr",nPerm = 0)
ddCor_resSig <- ddcorFindSignificant(ddCor_res,adjusted=TRUE,classes=FALSE)
length(intersect(testdcg,unique(sort(unlist(ddCor_resSig)))))
#[1] 1957 overlap with aggr dcgs
overlapInPart <- data.frame(Gene = degXdcg.aggr,DiffCorr.A = as.numeric(rep("0",times = length(degXdcg.aggr))),DiffCorr.D = as.numeric(rep("0",times = length(degXdcg.aggr))))
overlapInPart$DiffCorr.A[paste0(overlapInPart$Gene,"_A") %in% unique(sort(unlist(ddCor_resSig)))] = 1
overlapInPart$DiffCorr.D[paste0(overlapInPart$Gene,"_D") %in% unique(sort(unlist(ddCor_resSig)))] = 1
write.table(overlapInPart,"overlapDEGxDCGAggrAndPart.txt",sep = "\t",row.names = FALSE, quote = FALSE)

#
library(topGO)
geneID2GO <- readMappings(file = "G.raimondii_JGI_221_v2.1_genes2GO.map", sep = " ", IDsep =",")
geneNames <- colnames(multiExpr[[3]]$data)

sig.list <- dcgWildvDom
geneList <- factor(as.integer(geneNames %in% sig.list))
names(geneList) <- geneNames

pdf(file="DE-DCresult.WildvDom.pdf")
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
    enrichDC <- enrich
    GOresults<-rbind(GOresults,enrichDE)
}


##run again with DExDC
load("D5.annotation.RData")

DEDCannot <- annotation221[annotation221$gene %in% degXdcg,]
#There are doubles from multiple transcripts
DEDCannot[!duplicated(DEDCannot$gene),]
DExDC_modules <- as.data.frame(cbind(degXdcg,moduleLabels.adjusted[match(degXdcg,colnames(multiExpr[[3]]$data))]))
table(as.numeric(DExDC_modules$V2))
# 1  2  3  4  5  7  8  9 10 11 13 14 15 16 17 18 20 21 22 23 24 25 26 27 29 30
#75 23  1 19 15  1  9 11 14  4 18 13  1  2  6 10  3  5  6  8  8  7 11  3  7  5
#31 32 33 36 37 38 40 43 44 45 47
# 2 12  3  3 87  2  2  1  2  4  4

DExDCbyCon<-list()
for(i in 1:length(degByCon)){bool = degByCon[[i]] %in% dcgWildvDom; DExDCbyCon[[i]] = degByCon[[i]][bool]}
DExDCmelt <- rep("",times = length(degXdcg))
for(i in 1:length(degXdcg))
{
	for(j in 1:length(DExDCbyCon)){if(any(grepl(degXdcg[i],DExDCbyCon[[j]]))){DExDCmelt[i] <- paste(DExDCmelt[i],names(DExDCbyCon)[j],sep=",")}}
}

DExDCmelt<- gsub(pattern = "^,", rep ="", x= DExDCmelt, perl = TRUE)
DExDCmelt <- cbind(degXdcg,DExDCmelt)
DExDCmelt <- as.data.frame(DExDCmelt)
colnames(DExDCmelt) <- c("Gene","sig.DE")
DExDC_mod_de <- merge(DExDC_modules,DExDCmelt)
m <- merge(DEDCannot,DExDC_mod_de,by = "Gene")
DEDCannotwcModDE <- m
write.table(DEDCannotwcModDE,file = "DExDCannotation.cModDE.txt",sep = "\t",quote= FALSE,row.names = FALSE)

#Need to get directionality .. doh!
#DE.aggr.devByCon has DE
degByConDir <- list()
for(i in 1:length(DE.aggr.devByCon)){
	bool.up <- DE.aggr.devByCon[[i]]$padj < 0.05 & DE.aggr.devByCon[[i]]$log2FoldChange > 0
	bool.up[is.na(bool.up)] <- FALSE
	comp.up <- rownames(DE.aggr.devByCon[[i]])[bool.up]
	bool.down <- DE.aggr.devByCon[[i]]$padj < 0.05 & DE.aggr.devByCon[[i]]$log2FoldChange < 0
	bool.down[is.na(bool.down)] <- TRUE
	comp.down <- rownames(DE.aggr.devByCon[[i]])[bool.down]
	comp.name <- strsplit(attributes(DE.aggr.devByCon[[i]])[5]$elementMetadata[2,2],": ")[[1]][2]
	comp.name <- gsub(pattern  = "sample0 ",rep = "", comp.name)
	degByConDir[[i]] <- list(comp.up, comp.down)
	names(degByConDir)[i] <- comp.name
}

DExDCbyConDir<-list()
for(i in 1:length(degByConDir)){
	bool.up = degByConDir[[i]][[1]] %in% dcgWildvDom
	bool.down = degByConDir[[i]][[2]] %in% dcgWildvDom
	DExDCbyConDir[[i]] = list(degByConDir[[i]][[1]][bool.up],degByConDir[[i]][[2]][bool.down])
	names(DExDCbyConDir)[i] <- names(degByConDir)[i]
}

DExDCmeltUp <- rep("",times = length(degXdcg))
DExDCmeltDown <- rep("",times = length(degXdcg))
for(i in 1:length(degXdcg))
{
	for(j in 1:length(DExDCbyConDir)){
		if(any(grepl(degXdcg[i],DExDCbyConDir[[j]][[1]]))){DExDCmeltUp[i] <- paste(DExDCmeltUp[i],names(DExDCbyConDir)[j],sep=",")}
		if(any(grepl(degXdcg[i],DExDCbyConDir[[j]][[2]]))){DExDCmeltDown[i] <- paste(DExDCmeltDown[i],names(DExDCbyConDir)[j],sep=",")}
	}
}

DExDCmeltUp<- gsub(pattern = "^,", rep ="", x= DExDCmeltUp, perl = TRUE);
DExDCmeltDown<- gsub(pattern = "^,", rep ="", x= DExDCmeltDown, perl = TRUE);
DExDCmeltUp <- cbind(degXdcg,DExDCmeltUp)
DExDCmeltDown <- cbind(degXdcg,DExDCmeltDown)
DExDCmeltDir <- merge(as.data.frame(DExDCmeltUp),as.data.frame(DExDCmeltDown),by = "degXdcg")
colnames(DExDCmeltDir) <- c("Gene","sigUpDE","sigDownDE")

dcgWildvDomDir <- data.frame(Gene = dcgWildvDom, Mean = 0, UpCount = 0, DownCount = 0)
for(i in 1:length(dcgWildvDom)){
	mol1 = grep(dcgWildvDom[i],x$molecule.X)
	mol2 = grep(dcgWildvDom[i],x$molecule.Y)
	avg = mean(c(x[mol1,]$X.r1.r2.,x[mol2,]$X.r1.r2))
	count.tab <- table(sign(c(x[mol1,]$X.r1.r2.,x[mol2,]$X.r1.r2)))
	dcgWildvDomDir$Mean[i] = avg
	dcgWildvDomDir$UpCount[i] = as.numeric(count.tab["1"])
	dcgWildvDomDir$DownCount[i] = as.numeric(count.tab["-1"])
}
dcgWildvDomDir$UpCount[is.na(dcgWildvDomDir$UpCount)] = 0
dcgWildvDomDir$DownCount[is.na(dcgWildvDomDir$DownCount)] = 0
degxdcgDirDat <- dcgWildvDomDir[match(degXdcg,dcgWildvDomDir$Gene),]


DExDC_mod_de <- merge(DExDC_modules,DExDCmeltDir)
DExDC_mod_de_diffcorr <- merge(DExDC_mod_de,degxdcgDirDat)
DEDCannotwcModDECor <- merge(DEDCannot,DExDC_mod_de,by = "Gene")
write.table(DEDCannotwcModDECor,file = "DExDCannotation.cModDECor.txt",sep = "\t",quote= FALSE,row.names = FALSE)

