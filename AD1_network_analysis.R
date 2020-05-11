##AD1Network.R
##AD1 network analysis script
library(DESeq2)
library(flashClust)
library(WGCNA)
library(RColorBrewer)
library(ggplot2)
library(Rmisc)
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads=8)

counts.AD.expr <- counts.AD.raw[!(rowSums(counts.AD.raw) %in% 0),]
dim(counts.AD.expr)

#normalize read counts with rlog
AD1.dds <- DESeqDataSetFromMatrix(countData = counts.AD.expr, colData = colMat, design = ~ sample0)
AD.rld <- rlog(AD1.dds,blind=FALSE)
dim(assay(AD.rld))
library(WGCNA)
enableWGCNAThreads(nThreads=8)
#Euclidean distance matrix calculation for samples
A <- adjacency((assay(AD1.rld)),type="distance")
#calculate connectivity
k=as.numeric(apply(A,2,sum))-1
#standardized connectivity
Z.k=scale(k)

#set threshold for distance matrix
thresholdZ.k=-5
outlierColor1=ifelse(Z.k<thresholdZ.k,"red","black")

thresholdZ.k= -2.5
outlierColor2=ifelse(Z.k<thresholdZ.k ,"red","black")

#create the clustering tree
sampleTree = flashClust(as.dist(1-A), method = "average")
datColors=data.frame(outlierC_5=outlierColor1,outlierC_2.5=outlierColor2)

# Plot the sample dendrogram and the colors underneath
pdf("AD.part.dendrogram_and_trait_heatmap.rlog.pdf")
plotDendroAndColors(sampleTree,groupLabels=names(datColors),
colors=datColors,main="Sample dendrogram and outlier heatmap")
dev.off()

#Normalized expression is in AD.rld
#Set accession specific data sets
AD.rldDat <- as.data.frame(assay(AD.rld))
ADrldT <- t(AD.rldDat)#need to samples in rows, genes in columns for WGCNA
Tx2095data <- ADrldT[grep("TX2095",rownames(ADrldT)),]
Tx665data <- ADrldT[grep("TX665",rownames(ADrldT)),]
Yucdata <- ADrldT[c(grep("Yuc",rownames(ADrldT)),grep("TX2094",rownames(ADrldT))),]
Maxxadata <- ADrldT[grep("Maxxa",rownames(ADrldT)),]
CRBdata <- ADrldT[grep("CRB252",rownames(ADrldT)),]
TM1data <- ADrldT[grep("TM1",rownames(ADrldT)),]

#Now make the three data sets
dataDom <- rbind(CRBdata,Maxxadata,TM1data)
dataWild <- rbind(Tx2095data,Tx665data,Yucdata)
dataAll <- ADrldT

nSets = 3
setLabels = c("Domesticated G. hirsutum", "Wild G. hirsutum", "All G. hirsutum")
shortLabels = c("Dom", "Wild", "All")
multiExpr = vector(mode = "list",length = nSets)
multiExpr[[1]] = list(data = dataDom)
multiExpr[[2]] = list(data = dataWild)
multiExpr[[3]] = list(data = dataAll)

#make metadata object
accession <- c(rep("AD1.CRB252",4),rep("AD1.Maxxa",4),rep("AD1.TM1",4),rep("AD1.TX2095",4),
               rep("AD1.TX665",4),rep("AD1.Yuc",4))
timepoint <- rep(c(5,10,15,20),6)
cult <- c(rep("D",12),rep("W",12))

colMat <- data.frame(accession = accession, timepoint = factor(timepoint), cult = factor(cult))
colMat$sample0 <- paste(colMat$cult,colMat$timepoint,sep = ".")

#Confirm sample and gene quality
gsg <- goodSamplesGenesMS(multiExpr, verbose =3)

#check on data sets and clustering
pdf(file = "SampleClustering.pdf", width = 12, height = 12);
par(mfrow=c(nSets,1))
par(mar = c(0, 4, 2, 0))
sampleTrees = list()
for (set in 1:nSets)
{
    sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}
for (set in 1:nSets)
plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
xlab="", sub="", cex = 0.7);
dev.off()

save(multiExpr,nSets, setLabels, shortLabels, colMat, file = "R-01-dataInput.RData")
#check point save

##Part 2: Thresholding and Power

#load prior data set if not already loaded
nGenes<-checkSets(multiExpr)$nGenes

# Choose a set of soft-thresholding powers
types<-c("signed")
for (type in types)
{
    powers = c(c(1:10), seq(from = 12, to=40, by=2))
    # Initialize a list to hold the results of scale-free analysis
    powerTables = vector(mode = "list", length = nSets);
    # Call the network topology analysis function for each set in turn
    for(set in 1:nSets){
        powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers, verbose = 2, networkType = type)[[2]])
	}
    collectGarbage()
    
    # Plot the results:
    colors=brewer.pal(5,"Set1")
    # Will plot these columns of the returned scale free analysis tables
    plotCols = c(2,5,6,7)
    colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity", "Max connectivity");
    # Get the minima and maxima of the plotted points
    ylim = matrix(NA, nrow = 2, ncol = 4);
    for (set in 1:nSets)
    {
        for (col in 1:length(plotCols))
        {
            ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
            ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
        }
    }
    
    # Plot the quantities in the chosen columns vs. the soft thresholding power
    sizeGrWindow(8, 6)
    pdf(paste("s2.ChooseSoftThresholdPower_withMJTx2094_",gsub(".* ","", type), ".pdf", sep="") )
    par(mfcol = c(2,2));
    par(mar = c(4.2, 4.2 , 2.2, 0.5))
    cex1 = 0.7;
    for (col in 1:length(plotCols)) for (set in 1:nSets)
    {
        if (set==1)
        {
            plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2], xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col], main = colNames[col]);
            addGrid()
        }
    if (col==1)
    {
        text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
        labels=powers,cex=cex1,col=colors[set]);
    } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]], labels=powers,cex=cex1,col=colors[set]);
    if (col==1)
    {
        legend("bottomright", legend = shortLabels, col = colors, pch = 20) ;
    } else
    legend("topright", legend = shortLabels, col = colors, pch = 20) ;
    }
    
    dev.off()
    assign(paste("powerTables.",gsub(".* ","", type),sep=""),powerTables)
}

save(powerTables.signed , file = "R-02-choosePower.wTX2094.RData")
#A power of 12, as suggested by Horvath, seems good

##Network Construction

#for aggregate
powers = c(12,16,20)
#for partitioned
powers = 12

#building network
for (set in 1:nSets )
{
    # Extract total read counts for each genome
    subDat    <-  multiExpr[[set]]$data
    subDat   <-  apply(subDat,2,as.numeric)  # important, otherwise report error
    for (j in powers )
    {
        softPower = j
        print(paste("Start building network for ",shortLabels[set]," using soft threshold ",j,"......",sep=""))
        # Network construction
        
        net = blockwiseModules(
             # Input data
             subDat,
             # Data checking options
             checkMissingData = TRUE,
             
             # Options for splitting data into blocks
             blocks = NULL,
             randomSeed = 12345,
             maxBlockSize = 20000,  # 5000 for 4G memory, 20000 for 16G, 30000 for 32 G
             
             # Network construction arguments: correlation options, use bicor instead of default pearson
             corType = "pearson",
             # Adjacency and topology overlap function options
             power = j, networkType = "signed", TOMType = "signed",
             
             # Saving or returning TOM
             saveTOMs = TRUE,
             saveTOMFileBase = paste(shortLabels[set],"_power",j,"_TOM",sep=""),
             
             # Basic tree cut options
             deepSplit = 2,  #default, known to reasonable
             minModuleSize = min(30, ncol(subDat)/2 ), #default 20, use 30 for transcriptome
             pamStage = TRUE, pamRespectsDendro = TRUE, #default, known to reasonable
             
             # Thredhold to merge modules: a height cut of 0.25 corresponding to correlation of 0.75
             mergeCutHeight = 0.25,
             
             # others
             reassignThreshold = 0,
             numericLabels = TRUE,
             verbose = 3)
             
        assign(paste(shortLabels[set],"net",j,sep=""), net)
        }
}
save(list=grep(".+net.+",ls(), value=TRUE), file = "R-03-buildNetwork.wTX2094.RData")

##Network topology analysis

powers<-c(12)
pname<-paste("power=", powers, sep="")


# work with individual genome, then work with different soft threshold

for (set in 1:nSets )
{
    # Extract total read counts for each genome
    subDat    <-  multiExpr[[set]]$data
    genome <-  shortLabels[set]
    net1<-get(paste(genome,"net",powers[1],sep="") )
    gg<-net1$goodGenes    #get the good genes descision
    
    pdf(paste("s4.",genome, "_connectivity.pdf", sep=""))

    # The following computes the network connectivity (Connectivity)
    degree<-as.data.frame(matrix(NA,nrow=nGenes,ncol=length(powers)))
    names(degree)=pname
    for(j in 1:length(powers) ){
        k= softConnectivity(subDat,power=powers[j], type = "signed")
        # Creating Scale Free Topology Plots
        scaleFreePlot(k, truncated=T,main= paste("beta=",powers[j]))
        degree[,j]<-k
    }
    
    # Plot dendrograms all three clustering results
    plotDendroAndColors( net1$dendrograms[[1]], main = paste( pname[1], "dendrogram" ), labels2colors(net1$colors[gg]),pname,dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)

    dev.off()
}

#FOR MULTIBLOCK NETWORKS
for (set in 1:nSets )
{
    # Extract total read counts for each genome
    subDat    <-  multiExpr[[set]]$data
    genome <-  shortLabels[set]
    net1<-get(paste(genome,"net",powers[1],sep="") )
    gg<-net1$goodGenes    #get the good genes descision
    
    pdf(paste("s4.",genome, "_connectivity_wTX2094_partitioned.pdf", sep=""))

    # The following computes the network connectivity (Connectivity)
    degree<-as.data.frame(matrix(NA,nrow=nGenes,ncol=length(powers)))
    names(degree)=pname
    for(j in 1:length(powers) ){
        k= softConnectivity(subDat,power=powers[j], type = "signed")
        # Creating Scale Free Topology Plots
        scaleFreePlot(k, truncated=T,main= paste("beta=",powers[j]))
        degree[,j]<-k
    }
    
    # Plot dendrograms all three clustering results
    for(l in 1:length(net1$dendrograms))
	{
		plotDendroAndColors( net1$dendrograms[[l]], main = paste( "block",l, "dendrogram" ), labels2colors(net1$colors[net1$blockGenes[[l]]]),pname,dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
	}
    dev.off()
}

# plot all powers for all genomes
pdf("s4.Allsets_Connectivity.pdf")
for(j in 1:length(powers))
{
    degree<-as.data.frame(matrix(NA,nrow=nGenes,ncol=nSets))
    names(degree)<- shortLabels
    for(set in 1:nSets)
    {
        subDat <-  multiExpr[[set]]$data
        k= softConnectivity(subDat,power=powers[j], type = "signed")
        degree[,set]<-k
    }
    multiScaleFreePlot(degree, paste("Connectivity distribution, Power=",powers[j],sep=""))
}
dev.off()

samples = c("D05","D10","D15","D20","D05","D10","D15","D20","D05","D10","D15","D20","W05","W10","W15","W20","W05","W10","W15","W20","W05","W10","W15","W20")
sampleList = list(samples[1:12],samples[13:24],samples)

for (set in 1:3 ) # all
{
    # Extract total read counts for each genome
    subDat    <-  multiExprAll[[set]]$data
    subDat   <-  apply(subDat,2,as.numeric)  # important, otherwise report error
    genome <-  shortLabels[set]
    net<-get(paste(genome,"net12",sep="") )
    gg<-net$goodGenes    #get the good genes descision
    Nmodules= dim(net$MEs)[2]
    assign(paste("net", genome, sep=""),net)
    
    print(net$TOMFiles)
    print(paste("Number of modules in ",genome," network is ",Nmodules,sep=""))
    colorsa<- labels2colors(net$colors)
    
    pdf(paste("s4.",genome,"_modules.pdf",sep=""))
    
    # Eigengenes are the 1st principal component of modules in a given single dataset, which provide a summary profile for each module.
    # Displaying module heatmap and the eigengene
    MEs<-net$MEs
    plots <- list()  # new empty list
    ss = as.factor(sampleList[[set]])
    for(me in 0:(Nmodules-1)) {
       which.module=paste("ME",me,sep="")
       module.color=labels2colors(me)
       #heatmap
       par(mfrow=c(2,1), mar=c(0.3, 5.5, 4, 2))
       plotMat(t(scale(subDat[,net$colors==me ]) ),
               nrgcols=30,rlabels=T,rcols=module.color,
               main=paste(which.module, module.color, sep=": "), cex.main=2)
       #barplot
       par(mar=c(5, 4.2, 0, 0.7))
       barplot(MEs[,which.module], col=module.color, main="", cex.main=2,
       ylab="eigengene expression",xlab="fiber development (dpa)", names.arg=as.character(ss))
       #line, anova
       df<-data.frame(ME=MEs[,which.module], ss, module = which.module )
       fit<-aov(ME~ss,df)
       dfc<-summarySE(df, measurevar="ME", groupvars=c("ss", "module"))
       dfc$genome <- gsub("05|10|15|20","",dfc$ss)
       dfc$dpa <-gsub("D|W","",dfc$ss)
       plots[[me+1]]<- ggplot(dfc, aes(x=genome, y=ME, fill = dpa)) +
                       geom_bar(stat="identity",position=position_dodge(), color="black", size=0.3) +
       geom_errorbar(aes(ymin=ME-se, ymax=ME+se), width=.3,position=position_dodge(0.9)) +
       ggtitle(paste(which.module," ",module.color,", P=", round(anova(fit)$"Pr(>F)"[1], 4), sep="") )+
       theme_bw() +
       theme(plot.title=element_text( size=11),legend.position = "none")
    }
    for(page in 1:ceiling(Nmodules/9))
    {
        if(Nmodules>(9*page))
        {  multiplot(plotlist = plots[(9*page-8):(9*page)],  layout=matrix(1:9, nrow=3, byrow=TRUE) )  }
        else
        {  multiplot(plotlist = plots[(9*page-8):Nmodules],  layout=matrix(1:9, nrow=3, byrow=TRUE) )  }
    }

   plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 90)

   dev.off()
}

save(sampleList, netDom, netWild, netAll, file = "R-04-networkTopology.RData")

#########
##Topology interpretation
#########

library(WGCNA);
library(flashClust);
library(RColorBrewer);
options(stringsAsFactors = FALSE);

remove(list=ls())

lnames = load(file = "R-01-dataInput.RData")
lnames     #  "multiExpr"   "nSets"       "setLabels"   "shortLabels", "datTraits"
nSets      # 3
setLabels  #  "Five genomes"         "Additional synthetic" "Proteomic set"
shortLabels# "setT" "setF" "setP"
nGenes<-checkSets(multiExpr)$nGenes #36560
nGenes


load('R-04-networkTopology.RData')

net<-netAll
# get eigengenes
MEs<-net$MEs
# eigengene~sample, anova
samples <- sampleList[[3]]
ss<-as.factor(samples)
pval<-apply(MEs,2,function(x){round(anova(aov(x~ss) )$"Pr(>F)"[1],4)}) 
pval<-as.data.frame(pval)
pval$symbol<-ifelse(pval$pval<0.05,"*"," ")
pval$numeric<-as.numeric(substring(rownames(pval),3) )
pval<-pval[order(pval$numeric),]
pval$symbol[1]<-" "  # ME0 always meaningless
pval
#The significance measure is showing if module eigengene are controlled by dom status and timepoint

# For each module, we also define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profile. This allows us to quantify the similarity of all genes on the array to every module.
# calculate the module membership values
# (aka. module eigengene based connectivity kME):
MM=signedKME(multiExpr[[3]]$data, MEs)
rownames(MM)<-names(multiExpr[[3]]$data)
# Intramodular analysis: identifying genes with high MM
# Using the MM measures, we can identify genes that have a high significance for weight as well as high module membership in interesting modules.

####
##Annotation
####

aa<-load('D5annotation.Rdata')#from CottonGen
aa
aa<-annotation221[,c("transcript","tair10.defline", "pfam","panther","kog", "kegg.ec", "kegg.orthology", "go", "tair10", "tair10.symbol")]
dim(aa<-aa[grep("[.]1$",aa$transcript),])
aaA <- aa; aaD <- aa;
aa$gene = gsub("[.]1$","", aa$transcript)
aaA$gene = paste(gsub("[.]1$","", aaA$transcript),"_A",sep="")
aaD$gene = paste(gsub("[.]1$","", aaD$transcript),"_D",sep="")
aaAll<-rbind(aaA,aaD)
me<-data.frame(gene = colnames(multiExpr[[3]]$data), ME = net$colors)
dim(me<-merge(me, aa, all.x=TRUE, by="gene"))
write.table(me, file="s5.moduleAndannotation.partitioned.txt", row.names=FALSE,sep="\t")

## Functional enrichment analysis
# add Gorai id to module group
MM$ID<-colnames(multiExpr[[3]]$data)
for(me in module.sig)
{
    me<-paste("k",me,sep="")
    write.table(MM[,c("ID",me)], file=paste("GSEA/MM/rnk_",me,".rnk",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}
# Run GseaPreranked analysis
MM<-MM[,order(names(MM))]
write.table(MM,file="s5.moduleMembership.txt",sep="\t",row.names=FALSE,quote=FALSE)
save(net, colMat, pval, MM, file = "MM.RData")


me0<-me
names(MM)[1]<-"gene"
dim(me0<-merge(me0,MM,by="gene",all.x=TRUE, all.y=TRUE) )

write.table(me0,file="s5.moduleMembershipAndAnnotation.partitioned.txt",sep="\t",row.names=FALSE)
write.table(t(multiExpr[[3]]$data),file="s5.rldExpr_all.txt",sep="\t",quote = FALSE)

#Wild vs Dom conservation
checkSets(multiExpr)
nSets<-checkSets(multiExpr)$nSets  
nGenes<-checkSets(multiExpr)$nGenes
nSamples<-checkSets(multiExpr)$nSamples[1]

softPower = 12
shortLabels = c("Dom","Wild")
# put all net in one list in correct order
nets<-list(Dom = Domnet12, Wild = Wildnet12)
names(nets)==shortLabels  # TRUE make sure correct order
all.colors<-cbind(Domnet12$colors, Wildnet12$colors)
# put all anova p together
# 12 samples
dpa12=as.factor(rep(seq(5,20,by=5),3))
# make dpa data frame
dpas<-list(Dom = dpa12, Wild = dpa12)
anovaP<-list()
for(i in 1:nSets)
{
    MEs<-nets[[i]]$MEs
    dpa<-dpas[[i]]
    pval<-apply(MEs,2,function(x){round(anova(aov(x~dpa) )$"Pr(>F)"[1],4)})
    pval<-as.data.frame(pval)
    pval$symbol<-ifelse(pval$pval<0.05,"*"," ")
    pval$numeric<-as.numeric(substring(rownames(pval),3) )
    pval<-pval[order(pval$numeric),]
    pval$symbol[1]<-" "  # ME0 always meaningless
    anovaP[[i]]<-pval
}
names(anovaP)<-shortLabels

######### End of my consensus analysis unit ###########

###JOE's unique partitioned stuff###
#1) I am trying to get look at homoeolog bias within partitioned network
#2) As well as trying to get the differences in modules bt aggr and partitioned
modDat <- as.data.frame(cbind(colnames(multiExpr[[3]]$data),Allnet12$colors))
table(factor(modDat[grep(pattern = "*_A",modDat[,1]),2])) #gets A contribution to each module

table(factor(modDat[grep(pattern = "*_D",modDat[,1]),2]))

table(factor(modDat[,2]))

modPairCounts <- as.data.frame(rbind(table(factor(modDat[grep(pattern = "*_A",modDat[,1]),2])),table(factor(modDat[grep(pattern = "*_D",modDat[,1]),2]))),row.names = c("A","D"))
modPairCountsT <- t(modPairCounts)
modADat <- modDat[grep(pattern = "*_A",modDat[,1]),]
modADat <- modDat[grep(pattern = "*_D",modDat[,1]),]
modADat$V1 <- gsub(pattern = "_A", rep = "", x = modADat[,1])
modDDat$V1 <- gsub(pattern = "_D", rep = "", x = modDDat[,1])
modADDat <- merge(modADat,modDDat,all = TRUE)
modADDat$moduleIdentity <- modADDat$A.module == modADDat$D.module
table(modADDat$moduleIdentity)

#modulePreservation() analysis between the At network and Dt network, to see if the
#module structure of the aggregate network is preserved in each of the partitioned networks
#utilizes the module structure from aggregate network and the multiExpr of the At genes

splitExpr <- vector("list",3)
splitExpr[[1]]$data = aggrAllExpr$data #need gene expr from aggregate dataset
splitExpr[[2]]$data = multiExpr[[3]]$data[,A.geneSet]
splitExpr[[3]]$data = multiExpr[[3]]$data[,D.geneSet]
temp <- gsub(pattern = "_A", rep = "", x = colnames(splitExpr[[2]]$data))
colnames(splitExpr[[2]]$data) <- temp
temp <- gsub(pattern = "_D", rep = "", x = colnames(splitExpr[[3]]$data))
colnames(splitExpr[[3]]$data) <- temp

splitColor<-list(Aggr=aggrColorDat[,2])#needs to be just the numbers, but they need to be correlated

pwset<-combn(2,2)

nPermutations1=200

set.seed(1)
system.time({
    mp = modulePreservation(splitExpr, splitColor, networkType="signed", referenceNetworks = 1, testNetworks=c(2,3), nPermutations = nPermutations1,
    randomSeed = 1, quickCor = 0, verbose = 3)
})#This works bc the comparison only occurs between those that are in both the ref and test set

pdf("medianRankAndZsummaryForAggrVsPart.pdf")
par(mfrow=c(2,2),mar = c(4.5,4.5,2.5,1))
ref=1; test = c(2,3)
for(te in test){
	Obs.PreservationStats= mp$preservation$observed[[ref]][[te]]
	Z.PreservationStats=mp$preservation$Z[[ref]][[te]]
	# Let us now visualize the data.
	modColors = rownames(Obs.PreservationStats)
	moduleSize = Obs.PreservationStats$moduleSize
	# we will omit the grey module (background genes)
	# and the gold module (random sample of genes)
	selectModules = !(modColors %in% c("0", "gold"))
	# Text labels for points
	point.label = modColors[selectModules]
	#Composite preservation statistics
	medianRank=Obs.PreservationStats$medianRank.pres
	Zsummary=Z.PreservationStats$Zsummary.pres
	# plot medianRank versus module size
	plot(moduleSize[selectModules],medianRank[selectModules],col=1,
	bg=modColors[selectModules],pch = 21,main=paste(names(splitExpr)[te],"medianRank
	Preservation",sep=" "),
	cex = 2, ylab ="medianRank",xlab="Module size", log="x")
	labelPoints(moduleSize[selectModules],medianRank[selectModules],
	point.label,cex=1,offs=0.03)
	# plot Zsummary versus module size
	plot(moduleSize[selectModules],Zsummary[selectModules], col = 1,
	bg=modColors[selectModules],pch = 21,
	main=paste(names(splitExpr)[te],"Zsummary Preservation",sep=" "),
	cex=2,ylab ="Zsummary", xlab = "Module size", log = "x")
	labelPoints(moduleSize[selectModules],Zsummary[selectModules],
	point.label,cex=1,offs=0.03)
	# Add threshold lines for Zsummary
	abline(h=0); abline(h=2, col = "blue", lty = 2)
	abline(h=10, col = "red", lty = 2)
}


pdf("Allnet12_module15_heatmap_and_scatterplot.pdf")
whichmodule="15"
Eigengene15 <- Allnet12$MEs[,colnames(Allnet12$MEs) == "ME15"]
datExprModule=multiExpr[[3]]$data[,moduleColorsAggr=="midnightblue"]
# set the margins of the graphics window
# par(mar=c(0.3, 5.5, 3, 2))
# create a heatmap whose columns correspond to the libraries
# and whose rows correspond to genes
plot.mat(t(scale(datExprModule)),cex.axis=2,nrgcols=30,rlabels=F,rcols="midnightblue",main=paste("heatmap of,"whichmodule,"module"))
#scatter plot between eigengene and sample network connectivity
verboseScatterplot(Eigengene15,Z.k,xlab=paste("ME",whichmodule,sep=""))
abline(h=-2,col="red",lwd=2)
dev.off()

#compare domesticated and wild networks
splitExpr <- list(Dom = multiExpr[[1]],Wild = multiExpr[[2]])
splitColor<-list(Dom = Domnet12$colors, Wild = Wildnet12$colors)

pwset<-combn(2,2)

nPermutations1=200

set.seed(1)
system.time({
    mp = modulePreservation(splitExpr, splitColor, networkType="signed", referenceNetworks = c(1,2), testNetworks= list(2,1), nPermutations = nPermutations1,
    randomSeed = 1, quickCor = 0, verbose = 3)
})

pdf("modPresGraphsForWildVsDomPart.pdf")
par(mfrow=c(2,2),mar = c(4.5,4.5,2.5,1))
comp <- list(c(1,2),c(2,1))
for(co in comp){
	Obs.PreservationStats= mp$preservation$observed[[co[1]]][[co[2]]]
	Z.PreservationStats=mp$preservation$Z[[co[1]]][[co[2]]]
	# Let us now visualize the data.
	modColors = rownames(Obs.PreservationStats)
	moduleSize = Obs.PreservationStats$moduleSize
	# we will omit the grey module (background genes)
	# and the gold module (random sample of genes)
	selectModules = !(modColors %in% c("0", "0.1"))
	# Text labels for points
	point.label = modColors[selectModules]
	#Composite preservation statistics
	medianRank=Obs.PreservationStats$medianRank.pres
	Zsummary=Z.PreservationStats$Zsummary.pres
	# plot medianRank versus module size
	plot(moduleSize[selectModules],medianRank[selectModules],col=1,
	bg=modColors[selectModules],pch = 21,main=paste(names(splitExpr)[co[2]],"medianRank
	Preservation",sep=" "),
	cex = 2, ylab ="medianRank",xlab="Module size", log="x")
	labelPoints(moduleSize[selectModules],medianRank[selectModules],
	point.label,cex=1,offs=0.03)
	# plot Zsummary versus module size
	plot(moduleSize[selectModules],Zsummary[selectModules], col = 1,
	bg=modColors[selectModules],pch = 21,
	main=paste(names(splitExpr)[co[2]],"Zsummary Preservation",sep=" "),
	cex=2,ylab ="Zsummary", xlab = "Module size", log = "x")
	labelPoints(moduleSize[selectModules],Zsummary[selectModules],
	point.label,cex=1,offs=0.03)
	# Add threshold lines for Zsummary
	abline(h=0); abline(h=2, col = "blue", lty = 2)
	abline(h=10, col = "red", lty = 2)
}
dev.off()

#Looking at setting up the GSEA for DE in hub genes
#Get lists of DE genes for specific and general conditions (wild/dom, dpa)
DEgenes <-  list()
DEgenes_dev <- vector()
DEgenes_con <- vector()
for(comp in 1:length(DE.aggr.devByCon)){
	bool <- DE.aggr.devByCon[[comp]]$padj < 0.05
	bool <- ifelse(is.na(bool),FALSE,bool)
	geneIDs <- rownames(DE.aggr.devByCon[[comp]][bool,])
	DEgenes[[comp]] <- geneIDs
	if(comp < 7){
		temp <- c(DEgenes_dev,geneIDs)
		DEgenes_dev <- unique(temp)
	}
	if(comp >= 7){
		temp <- c(DEgenes_con,geneIDs)
		DEgenes_con <- unique(geneIDs)
	}
}
DEgenes[[comp+1]] <- DEgenes_dev
DEgenes[[comp+2]] <- DEgenes_con

#Set up gene list as a .gmt file
#"name of comparison"	"ID"	"Gene 1"	"Gene 2"	"Gene 3"	...	"Gene n"

labs <- c("w05vw10","w10vw15","w15vw20","d05vd10","d10vd15","d15vd20","w05vd05","w10vd10","w15vd15","w20vd20","all_dev","all_status")
names(DEgenes) <- labs
lapply(DEgenes, write, "DE_bycomparison.gmt", append = TRUE, ncolumns = 10000) #Still need to add row titles, but can do it in text editor

#for partitioned network
DEgenes <-  list()
DEgenes_dev <- vector()
DEgenes_con <- vector()
for(comp in 1:length(DE.homoeolog.devByCond)){
	bool <- DE.homoeolog.devByCond[[comp]]$padj < 0.05
	bool <- ifelse(is.na(bool),FALSE,bool)
	geneIDs <- rownames(DE.homoeolog.devByCond[[comp]][bool,])
	DEgenes[[comp]] <- geneIDs
	if(comp < 7){
		temp <- c(DEgenes_dev,geneIDs)
		DEgenes_dev <- unique(temp)
	}
	if(comp >= 7){
		temp <- c(DEgenes_con,geneIDs)
		DEgenes_con <- unique(geneIDs)
	}
}
DEgenes[[comp+1]] <- DEgenes_dev
DEgenes[[comp+2]] <- DEgenes_con

labs <- c("w05vw10","w10vw15","w15vw20","d05vd10","d10vd15","d15vd20","w05vd05","w10vd10","w15vd15","w20vd20","all_dev","all_status")
names(DEgenes) <- labs
lapply(DEgenes, write, "DE_bycomparison.homoeolog.gmt", append = TRUE, ncolumns = 14000) #Still need to add row titles, but can do it in text editor
