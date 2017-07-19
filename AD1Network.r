##AD1Network.R
##AD1 network analysis script

rldExpr <- read.table(file="s5.rldExpr_all.txt")
#remove column with Yuc 10dpa
#now remove reads with less than 1 read
dim(rldExprTrunc[rowSums(counts(rldExprTrunc))/23 >= 1,])

#Part 1: set up data sets
options(stringsAsFactors = FALSE)

#build 
row.A <- paste(rownames(counts.A.raw),"_A",sep = "")
row.D <- paste(rownames(counts.D.raw),"_D",sep = "")
rownames(counts.A.raw) <- row.A
rownames(counts.D.raw) <- row.D
counts.AD.raw <- rbind(counts.A.raw,counts.D.raw)
counts.AD.expr <- counts.AD.raw[!(rowSums(counts.AD.raw) %in% 0),]
dim(counts.AD.expr)
#convert genes to rlog
library(DESeq2)
AD1.dds <- DESeqDataSetFromMatrix(countData = counts.AD.expr, colData = colMat, design = ~ sample0)
AD.rld <- rlog(AD1.dds,blind=FALSE)
dim(assay(AD.rld))
library(WGCNA)
enableWGCNAThreads(nThreads=8)#for biocrunch
#need to run on biocrunch or other HPC
#Euclidean distance matrix calculation for SAMPLES! NOT GENES!
A <- adjacency((assay(AD1.rld)),type="distance")
#calculate connectivity
k=as.numeric(apply(A,2,sum))-1
#standardized connectivity
Z.k=scale(k)

thresholdZ.k=-5
outlierColor1=ifelse(Z.k<thresholdZ.k,"red","black")

thresholdZ.k= -2.5
outlierColor2=ifelse(Z.k<thresholdZ.k ,"red","black")

library(flashClust)
#create the clustering tree
sampleTree = flashClust(as.dist(1-A), method = "average")
datColors=data.frame(outlierC_5=outlierColor1,outlierC_2.5=outlierColor2)
# Plot the sample dendrogram and the colors underneath.
pdf("AD.part.dendrogram_and_trait_heatmap.rlog.pdf")
plotDendroAndColors(sampleTree,groupLabels=names(datColors),
colors=datColors,main="Sample dendrogram and outlier heatmap")
dev.off()

#let's run this in aggregate and partitioned
#Hooray, the samples are all in-liers! Can run with all samples, unless I wish to remove any
#let's run this similar to how Jing did it for later ease
#So our normalized expression is in AD.rld
# We can turn this into accession specific data sets
#Like Jing did
AD.rldDat <- as.data.frame(assay(AD.rld))
ADrldT <- t(AD.rldDat)
Tx2095data <- ADrldT[grep("TX2095",rownames(ADrldT)),]
Tx665data <- ADrldT[grep("TX665",rownames(ADrldT)),]
Yucdata <- ADrldT[c(grep("Yuc",rownames(ADrldT)),grep("TX2094",rownames(ADrldT))),]
Maxxadata <- ADrldT[grep("Maxxa",rownames(ADrldT)),]
CRBdata <- ADrldT[grep("CRB252",rownames(ADrldT)),]
TM1data <- ADrldT[grep("TM1",rownames(ADrldT)),]
#reminder, we have colMat which corresponds with the rows in our ADrldT

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

#make colMat - metadata file
accession <- c(rep("AD1.CRB252",4),rep("AD1.Maxxa",4),rep("AD1.TM1",4),rep("AD1.TX2095",4),
               rep("AD1.TX665",4),rep("AD1.Yuc",4))
timepoint <- rep(c(5,10,15,20),6)
cult <- c(rep("D",12),rep("W",12))

colMat <- data.frame(accession = accession, timepoint = factor(timepoint), cult = factor(cult))
colMat$sample0 <- paste(colMat$cult,colMat$timepoint,sep = ".")

#We have our data sets, now we can look at the quality of those data sets
gsg <- goodSamplesGenesMS(multiExpr, verbose =3)
#No bad samples, no bad genes
#Don't need to worry about removing genes or samples
#check on data sets and clustering again
pdf(file = "s1.SampleClustering.pdf", width = 12, height = 12);
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
#save the data so we have a reset point

##Part 2: Thresholding and Power
library(WGCNA)
library(RColorBrewer)
library(ggplot2)
library(Rmisc)
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads=8)

#load prior data set if not already loaded
nGenes<-checkSets(multiExpr)$nGenes

# Choose a set of soft-thresholding powers, consider three type of adjacnecy tables, although I am going to use "signed" network for this analysis
#types<-c("unsigned", "signed", "signed hybrid")
types<-c("signed")
for (type in types)
{
    powers = c(c(1:10), seq(from = 12, to=40, by=2))
    # Initialize a list to hold the results of scale-free analysis
    powerTables = vector(mode = "list", length = nSets);
    # Call the network topology analysis function for each set in turn
    for(set in 1:nSets){
        powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers, verbose = 2, networkType = type)[[2]])      }
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

#results of pickSoftThreshold for domesticated network
#   Power SFT.R.sq   slope truncated.R.sq mean.k. median.k. max.k.
#1      1  0.04350  5.5700          0.925 14900.0   14900.0  15700
#2      2  0.04380  1.7000          0.925  8910.0    8930.0  10400
#3      3  0.02760  0.6660          0.845  5900.0    5870.0   7950
#4      4  0.00946  0.2380          0.814  4190.0    4110.0   6420
#5      5  0.00177 -0.0708          0.805  3110.0    3010.0   5360
#6      6  0.05660 -0.3140          0.807  2400.0    2270.0   4600
#7      7  0.18200 -0.4900          0.824  1900.0    1760.0   4010
#8      8  0.33200 -0.6260          0.848  1540.0    1390.0   3540
#9      9  0.46200 -0.7410          0.865  1270.0    1110.0   3160
#10    10  0.56100 -0.8410          0.877  1060.0     906.0   2860
#11    12  0.67500 -0.9990          0.888   771.0     617.0   2380
#12    14  0.72700 -1.1100          0.894   581.0     437.0   2020
#13    16  0.75900 -1.1900          0.900   451.0     317.0   1750
#14    18  0.78000 -1.2500          0.907   358.0     236.0   1530
#15    20  0.78800 -1.3100          0.908   290.0     179.0   1350
#16    22  0.79700 -1.3500          0.912   238.0     139.0   1210
#17    24  0.81400 -1.3800          0.924   199.0     109.0   1090
#18    26  0.82000 -1.4200          0.929   168.0      86.6    985
#19    28  0.82800 -1.4400          0.935   143.0      69.6    899
#20    30  0.83500 -1.4600          0.940   123.0      56.2    824
#21    32  0.83600 -1.4900          0.942   106.0      45.9    759
#22    34  0.84400 -1.5000          0.948    92.8      37.8    701
#23    36  0.84700 -1.5100          0.951    81.4      31.4    650
#24    38  0.85500 -1.5200          0.958    71.9      26.3    605
#25    40  0.85800 -1.5300          0.962    63.8      22.1    565

#results of pickSoftThreshold for wild network (no Yuc 10 dpa)
#   Power SFT.R.sq    slope truncated.R.sq mean.k. median.k. max.k.
#1      1 2.30e-02  6.46000          0.923 14900.0   14900.0  15500
#2      2 1.36e-02  1.10000          0.721  8840.0    8800.0  10300
#3      3 4.05e-06  0.00779          0.572  5810.0    5700.0   7820
#4      4 3.20e-02 -0.40700          0.545  4100.0    3940.0   6340
#5      5 1.35e-01 -0.61000          0.571  3030.0    2850.0   5330
#6      6 2.76e-01 -0.74200          0.630  2320.0    2120.0   4580
#7      7 4.28e-01 -0.84900          0.693  1830.0    1620.0   4010
#8      8 5.40e-01 -0.93600          0.737  1480.0    1260.0   3560
#9      9 6.24e-01 -1.01000          0.773  1220.0     994.0   3190
#10    10 6.74e-01 -1.07000          0.791  1020.0     797.0   2890
#11    12 7.49e-01 -1.16000          0.825   738.0     530.0   2400
#12    14 7.77e-01 -1.22000          0.836   556.0     366.0   2050
#13    16 7.81e-01 -1.26000          0.832   432.0     259.0   1770
#14    18 7.88e-01 -1.30000          0.838   343.0     188.0   1550
#15    20 7.95e-01 -1.32000          0.843   279.0     141.0   1360
#16    22 7.96e-01 -1.34000          0.844   230.0     107.0   1220
#17    24 8.13e-01 -1.34000          0.861   192.0      82.9   1090
#18    26 8.19e-01 -1.35000          0.868   162.0      65.4    983
#19    28 8.27e-01 -1.35000          0.877   139.0      52.4    892
#20    30 8.38e-01 -1.35000          0.887   119.0      42.2    812
#21    32 8.47e-01 -1.35000          0.896   104.0      34.4    743
#22    34 8.55e-01 -1.35000          0.905    90.8      28.3    683
#23    36 8.64e-01 -1.34000          0.913    79.9      23.5    629
#24    38 8.68e-01 -1.34000          0.918    70.8      19.6    581
#25    40 8.62e-01 -1.36000          0.915    63.0      16.5    543

#results of pickSoftThreshold for all network (no Yuc 10 dpa)
#Power SFT.R.sq slope truncated.R.sq mean.k. median.k. max.k.
#1      1   0.0738 12.80          0.958 14900.0  14900.00  15500
#2      2   0.0138 -1.31          0.867  8480.0   8430.00   9880
#3      3   0.1400 -1.71          0.825  5270.0   5160.00   7210
#4      4   0.2250 -1.36          0.815  3490.0   3350.00   5610
#5      5   0.3240 -1.20          0.832  2430.0   2270.00   4530
#6      6   0.4410 -1.18          0.853  1760.0   1590.00   3760
#7      7   0.5470 -1.19          0.873  1320.0   1140.00   3190
#8      8   0.6310 -1.22          0.888  1010.0    841.00   2750
#9      9   0.6810 -1.28          0.890   791.0    629.00   2410
#10    10   0.7200 -1.33          0.895   632.0    478.00   2130
#11    12   0.7730 -1.41          0.909   421.0    287.00   1710
#12    14   0.7980 -1.47          0.914   295.0    180.00   1410
#13    16   0.8060 -1.53          0.915   214.0    117.00   1190
#14    18   0.8180 -1.55          0.921   160.0     79.10   1020
#15    20   0.8220 -1.57          0.926   123.0     54.50    879
#16    22   0.8320 -1.58          0.935    96.3     38.20    767
#17    24   0.8360 -1.59          0.938    76.7     27.30    676
#18    26   0.8460 -1.58          0.947    62.0     19.70    600
#19    28   0.8480 -1.59          0.950    50.7     14.40    536
#20    30   0.8570 -1.59          0.957    42.0     10.70    482
#21    32   0.8670 -1.58          0.964    35.1      8.03    435
#22    34   0.8760 -1.57          0.970    29.7      6.09    394
#23    36   0.8830 -1.57          0.976    25.2      4.66    361
#24    38   0.8870 -1.57          0.979    21.6      3.61    333
#25    40   0.8870 -1.58          0.980    18.6      2.79    307

##results for domesticated expression when MJ's TX2094 is substituted in
##hopefully not too different

#pickSoftThreshold: will use block size 1506.
#   Power SFT.R.sq   slope truncated.R.sq mean.k. median.k. max.k.
#1      1  0.03930  5.7100          0.924 14900.0   14900.0  15700
#2      2  0.05950  2.1600          0.923  8910.0    8930.0  10400
#3      3  0.03140  0.7130          0.856  5910.0    5870.0   7940
#4      4  0.00941  0.2320          0.827  4190.0    4120.0   6400
#5      5  0.00133 -0.0603          0.818  3120.0    3020.0   5350
#6      6  0.05470 -0.3030          0.821  2410.0    2280.0   4580
#7      7  0.18100 -0.4770          0.840  1910.0    1770.0   4000
#8      8  0.33500 -0.6210          0.860  1550.0    1400.0   3540
#9      9  0.46600 -0.7500          0.868  1270.0    1120.0   3170
#10    10  0.56400 -0.8460          0.883  1070.0     912.0   2860
#11    12  0.67700 -1.0100          0.893   775.0     624.0   2390
#12    14  0.73100 -1.1100          0.901   585.0     442.0   2030
#13    16  0.76400 -1.1900          0.908   455.0     322.0   1760
#14    18  0.78400 -1.2600          0.913   361.0     240.0   1550
#15    20  0.80100 -1.3100          0.921   293.0     183.0   1370
#16    22  0.81000 -1.3600          0.925   241.0     142.0   1230
#17    24  0.81900 -1.4000          0.931   201.0     111.0   1110
#18    26  0.83100 -1.4300          0.939   170.0      88.3   1010
#19    28  0.83600 -1.4500          0.943   145.0      71.0    922
#20    30  0.84200 -1.4800          0.947   125.0      57.6    847
#21    32  0.84400 -1.5000          0.949   108.0      47.2    781
#22    34  0.84900 -1.5100          0.953    94.5      38.9    723
#23    36  0.85600 -1.5300          0.959    83.1      32.4    672
#24    38  0.86200 -1.5400          0.964    73.4      27.1    627
#25    40  0.86600 -1.5500          0.968    65.3      22.8    586

##Here is for wild with MJ's TX2094 10 dpa

#pickSoftThreshold: will use block size 1506.
#   Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
#1      1   0.0110  4.360          0.972 14900.0   14900.0  15600
#2      2   0.0561  2.120          0.778  8770.0    8740.0  10100
#3      3   0.0011  0.124          0.575  5700.0    5580.0   7630
#4      4   0.0378 -0.429          0.533  3960.0    3800.0   6150
#5      5   0.1370 -0.613          0.576  2900.0    2700.0   5140
#6      6   0.2930 -0.761          0.643  2200.0    1980.0   4400
#7      7   0.4550 -0.890          0.709  1720.0    1490.0   3840
#8      8   0.5790 -0.990          0.764  1370.0    1140.0   3400
#9      9   0.6690 -1.070          0.805  1120.0     892.0   3040
#10    10   0.7170 -1.130          0.823   926.0     706.0   2740
#11    12   0.7620 -1.220          0.833   660.0     457.0   2280
#12    14   0.7900 -1.280          0.844   491.0     308.0   1930
#13    16   0.8030 -1.310          0.851   377.0     214.0   1660
#14    18   0.8030 -1.350          0.851   296.0     153.0   1450
#15    20   0.8130 -1.360          0.860   238.0     112.0   1280
#16    22   0.8200 -1.370          0.867   194.0      84.0   1130
#17    24   0.8210 -1.380          0.871   161.0      64.4   1010
#18    26   0.8310 -1.380          0.882   135.0      49.8    912
#19    28   0.8440 -1.380          0.895   115.0      39.1    826
#20    30   0.8610 -1.370          0.911    98.2      31.1    751
#21    32   0.8590 -1.380          0.912    84.8      24.9    689
#22    34   0.8610 -1.390          0.918    73.8      20.2    636
#23    36   0.8570 -1.400          0.919    64.6      16.5    589
#24    38   0.8610 -1.410          0.924    56.9      13.6    547
#25    40   0.8480 -1.420          0.916    50.4      11.3    509

##Here is power selection for all with TX2094 10 dpa from MJ
## Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
##1      1 0.033300  9.260          0.965 14900.0  14900.00  15600
##2      2 0.000378 -0.237          0.886  8450.0   8410.00   9790
##3      3 0.069400 -1.280          0.796  5220.0   5110.00   7100
##4      4 0.219000 -1.340          0.804  3440.0   3290.00   5490
##5      5 0.339000 -1.230          0.824  2380.0   2220.00   4420
##6      6 0.463000 -1.210          0.851  1710.0   1540.00   3670
##7      7 0.576000 -1.250          0.875  1270.0   1100.00   3110
##8      8 0.652000 -1.290          0.886   972.0    805.00   2680
##9      9 0.714000 -1.330          0.898   759.0    598.00   2340
##10    10 0.754000 -1.360          0.908   603.0    452.00   2070
##11    12 0.801000 -1.420          0.920   400.0    268.00   1660
##12    14 0.828000 -1.470          0.928   278.0    166.00   1370
##13    16 0.848000 -1.500          0.938   201.0    107.00   1150
##14    18 0.858000 -1.510          0.944   150.0     71.40    977
##15    20 0.867000 -1.520          0.951   115.0     48.70    843
##16    22 0.873000 -1.530          0.955    89.4     33.80    735
##17    24 0.882000 -1.530          0.962    71.0     24.00    647
##18    26 0.888000 -1.530          0.968    57.3     17.20    576
##19    28 0.890000 -1.540          0.971    46.8     12.50    517
##20    30 0.888000 -1.550          0.971    38.7      9.25    467
##21    32 0.892000 -1.550          0.975    32.3      6.90    424
##22    34 0.896000 -1.550          0.978    27.2      5.22    388
##23    36 0.897000 -1.550          0.980    23.1      3.97    356
##24    38 0.901000 -1.560          0.983    19.8      3.04    328
##25    40 0.904000 -1.560          0.985    17.1      2.35    303
##

# repeat above for powerTables.signed, and powerTables.hybrid

#Results for partitioned dom network with TX2094
#   Power SFT.R.sq    slope truncated.R.sq mean.k. median.k. max.k.
#1      1 5.40e-02  5.81000          0.938 25700.0   25700.0  27400
#2      2 1.92e-01  3.87000          0.885 15200.0   15300.0  17600
#3      3 9.80e-02  1.35000          0.828 10000.0   10000.0  13200
#4      4 4.00e-02  0.50800          0.807  7030.0    6940.0  10600
#5      5 1.83e-05  0.00737          0.796  5180.0    5030.0   8800
#6      6 4.17e-02 -0.27600          0.796  3960.0    3770.0   7490
#7      7 1.64e-01 -0.47000          0.822  3110.0    2890.0   6480
#8      8 3.26e-01 -0.61700          0.849  2500.0    2260.0   5680
#9      9 4.65e-01 -0.72500          0.872  2050.0    1800.0   5040
#10    10 5.78e-01 -0.81600          0.892  1700.0    1450.0   4510
#11    12 6.90e-01 -0.99600          0.895  1220.0     974.0   3740
#12    14 7.39e-01 -1.13000          0.895   909.0     678.0   3170
#13    16 7.60e-01 -1.22000          0.893   699.0     487.0   2730
#14    18 7.69e-01 -1.30000          0.891   550.0     358.0   2380
#15    20 7.79e-01 -1.36000          0.897   442.0     269.0   2110
#16    22 7.87e-01 -1.41000          0.903   361.0     205.0   1880
#17    24 7.91e-01 -1.45000          0.907   299.0     159.0   1700
#18    26 7.92e-01 -1.49000          0.911   251.0     125.0   1540
#19    28 8.04e-01 -1.51000          0.921   212.0      99.2   1400
#20    30 8.06e-01 -1.54000          0.926   182.0      79.5   1280
#21    32 8.10e-01 -1.57000          0.931   157.0      64.4   1180
#22    34 8.12e-01 -1.59000          0.935   136.0      52.6   1100
#23    36 8.12e-01 -1.61000          0.938   119.0      43.4   1020
#24    38 8.16e-01 -1.62000          0.942   105.0      36.0    950
#25    40 8.19e-01 -1.64000          0.946    92.5      30.1    888

#Softhresholding for partitioned wild data set with TX2094
#   Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
#1      1  0.05260  7.510          0.973 25600.0   25600.0  27200
#2      2  0.11800  3.510          0.867 15000.0   15000.0  17300
#3      3  0.00116  0.141          0.709  9650.0    9500.0  12900
#4      4  0.03870 -0.473          0.641  6640.0    6390.0  10200
#5      5  0.11900 -0.625          0.637  4800.0    4500.0   8440
#6      6  0.27400 -0.813          0.688  3600.0    3270.0   7170
#7      7  0.44100 -0.956          0.750  2780.0    2440.0   6210
#8      8  0.57900 -1.070          0.799  2200.0    1850.0   5450
#9      9  0.66500 -1.150          0.827  1780.0    1430.0   4840
#10    10  0.72400 -1.210          0.848  1460.0    1120.0   4340
#11    12  0.78200 -1.290          0.867  1020.0     718.0   3560
#12    14  0.81600 -1.330          0.880   749.0     478.0   2990
#13    16  0.83600 -1.360          0.890   567.0     329.0   2550
#14    18  0.84700 -1.380          0.897   441.0     234.0   2200
#15    20  0.85000 -1.400          0.898   350.0     170.0   1930
#16    22  0.86300 -1.400          0.908   284.0     127.0   1700
#17    24  0.87300 -1.400          0.916   233.0      95.7   1510
#18    26  0.88700 -1.390          0.929   194.0      73.5   1350
#19    28  0.87800 -1.420          0.925   163.0      57.5   1230
#20    30  0.86600 -1.440          0.919   139.0      45.5   1120
#21    32  0.86200 -1.460          0.919   119.0      36.4   1030
#22    34  0.84800 -1.480          0.913   103.0      29.5    956
#23    36  0.82500 -1.520          0.898    89.6      24.0    887
#24    38  0.81100 -1.540          0.890    78.5      19.7    826
#25    40  0.80500 -1.560          0.890    69.2      16.3    771

##Softthresholding for partitioned all with TX2094
#   Power SFT.R.sq slope truncated.R.sq mean.k. median.k. max.k.
#1      1   0.0441  7.68          0.968 25600.0  25600.00  27100
#2      2   0.0123  1.32          0.914 14500.0  14400.00  16700
#3      3   0.0528 -1.19          0.827  8840.0   8700.00  11900
#4      4   0.1730 -1.31          0.810  5750.0   5540.00   9120
#5      5   0.3220 -1.31          0.829  3930.0   3670.00   7260
#6      6   0.4690 -1.31          0.857  2790.0   2520.00   5940
#7      7   0.5820 -1.34          0.877  2040.0   1770.00   4980
#8      8   0.6600 -1.39          0.887  1540.0   1280.00   4260
#9      9   0.7210 -1.42          0.899  1180.0    934.00   3690
#10    10   0.7580 -1.45          0.906   931.0    696.00   3230
#11    12   0.8100 -1.48          0.921   602.0    403.00   2540
#12    14   0.8440 -1.50          0.934   410.0    245.00   2050
#13    16   0.8520 -1.53          0.936   291.0    153.00   1710
#14    18   0.8600 -1.55          0.942   213.0     99.20   1450
#15    20   0.8570 -1.57          0.941   160.0     65.90   1240
#16    22   0.8690 -1.57          0.952   123.0     44.70   1080
#17    24   0.8770 -1.57          0.960    96.6     31.00    947
#18    26   0.8810 -1.57          0.964    77.0     22.00    837
#19    28   0.8880 -1.56          0.970    62.2     15.80    744
#20    30   0.8880 -1.57          0.973    50.9     11.50    672
#21    32   0.8880 -1.58          0.975    42.1      8.45    610
#22    34   0.8900 -1.59          0.978    35.1      6.27    557
#23    36   0.8910 -1.59          0.980    29.6      4.71    511
#24    38   0.8940 -1.59          0.983    25.1      3.56    470
#25    40   0.8980 -1.60          0.986    21.5      2.71    435


save(powerTables.signed , file = "R-02-choosePower.wTX2094.RData")

#This finished: now we have the power tables and the graphs
#A power of 12, as Jing suggests for signed networks (and Horvath does too), seems good

##NETWORK CONSTRUCTION
#We are just continuing so we don't need to reload our saved data
#reset powers with our determined value

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
    
    pdf(paste("s4.",genome, "_connectivity_wTX2094.pdf", sep=""))

    # Network topology and other concepts need to be explored more!!!!!!!!!!!!!!!!!!
    # The following computes the network connectivity (Connectivity)
    degree<-as.data.frame(matrix(NA,nrow=nGenes,ncol=length(powers)))
    names(degree)=pname
    for(j in 1:length(powers) ){
        k= softConnectivity(subDat,power=powers[j], type = "signed")
        # Creating Scale Free Topology Plots (SUPPLEMENTARY FIGURE S1 in our article)
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

    # Network topology and other concepts need to be explored more!!!!!!!!!!!!!!!!!!
    # The following computes the network connectivity (Connectivity)
    degree<-as.data.frame(matrix(NA,nrow=nGenes,ncol=length(powers)))
    names(degree)=pname
    for(j in 1:length(powers) ){
        k= softConnectivity(subDat,power=powers[j], type = "signed")
        # Creating Scale Free Topology Plots (SUPPLEMENTARY FIGURE S1 in our article)
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
##BACK TO NORMAL

# plot all powers for all genomes
pdf("s4.Allsets_Connectivity_wTX2094.pdf")
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
#change this 

for (set in 1:3 ) # all and trio
{
    # Extract total read counts for each genome
    subDat    <-  multiExprAll[[set]]$data
    subDat   <-  apply(subDat,2,as.numeric)  # important, otherwise report error
    genome <-  shortLabels[set]
    net<-get(paste(genome,"net12",sep="") )
    gg<-net$goodGenes    #get the good genes descision
#    adjacency = adjacency(subDat, power = 12, type = "signed")
    Nmodules= dim(net$MEs)[2]
    assign(paste("net", genome, sep=""),net)
    
#   load(net$TOMFiles)
    print(net$TOMFiles)
    print(paste("Number of modules in ",genome," network is ",Nmodules,sep=""))
    colorsa<- labels2colors(net$colors)
    
    pdf(paste("s4.",genome,"_modules.withTX2094.redo.pdf",sep=""))
    
    # Eigengenes are the 1st principal component of modules in a given single dataset, which provide a summary profile for each module.
    # Displaying module heatmap and the eigengene
    # sizeGrWindow(8,7);
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

   dev.off()# end the big plot
}


#Hooray, Modules!
#[1] "Dom_power12_TOM-block.1.RData"
#[1] "Number of modules in Dom network is 26"
#[1] "Wild_power12_TOM-block.1.RData"
#[1] "Number of modules in Wild network is 43"
#[1] "All_power12_TOM-block.1.RData"
#[1] "Number of modules in All network is 43"

#Modules in Aggregate with TX2094
#[1] "Dom_power12_TOM-block.1.RData"
#[1] "Number of modules in Dom network is 29"
#[1] "Wild_power12_TOM-block.1.RData"
#[1] "Number of modules in Wild network is 33"
#[1] "All_power12_TOM-block.1.RData"
#[1] "Number of modules in All network is 27"

#Module in Partitoned with TX2094
#[1] "Dom_power12_TOM-block.1.RData" "Dom_power12_TOM-block.2.RData"
#[3] "Dom_power12_TOM-block.3.RData"
#[1] "Number of modules in Dom network is 48"
#[1] "Wild_power12_TOM-block.1.RData" "Wild_power12_TOM-block.2.RData"
#[3] "Wild_power12_TOM-block.3.RData"
#[1] "Number of modules in Wild network is 108"
#[1] "All_power12_TOM-block.1.RData" "All_power12_TOM-block.2.RData"
#[3] "All_power12_TOM-block.3.RData"
#[1] "Number of modules in All network is 53"


save(sampleList, netDom, netWild, netAll, file = "R-04-networkTopology.RData")
#########
##Topology interpretation
#########

library(WGCNA);
library(flashClust);
library(RColorBrewer);
options(stringsAsFactors = FALSE);

remove(list=ls())#clean everything up

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
#Helps eliminate noise



# define a gene significance variable. We quantify associations of individual genes with our trait of interest (dpa) by defining Gene Significance GS as (the absolute value of) the correlation between the gene and the trait.
GS.dpa = cor(multiExpr[[3]]$data,colMat[,2],use="p")
GS.cult = cor(multiExpr[[3]]$data,as.numeric(colMat[,3]),use="p")#probably not kosher
#For some reason, it won't quantify the GS variable by cult
# This translates the numeric values into colors
GS.Color=numbers2colors(GS.dpa,signed=T)

# conduct DE analysis
## source('deseq2.r', chdir = TRUE)
#This is done, but need to 
# import differential expression results
load('DEid.Rdata')

pdf("Allnet_DendroWithAggrAndPartMods.pdf")
sigDEs2colors<-function(x,sigC) { ifelse( colnames(multiExpr[[3]]$data) %in% x, sigC, "white")}
DE_dev.Color   <- as.data.frame( sapply(DEid_dev,    function(x) sigDEs2colors(x, sigC="purple")) )
DE_con1.Color<- as.data.frame( sapply(DEid_con[1:3], function(x) sigDEs2colors(x, sigC="darkgreen")) )
DE_con2.Color<- as.data.frame( sapply(DEid_con[4:10], function(x) sigDEs2colors(x, sigC="brown")) )
#DE.homoeolog.devAll and *.devByCon from DESeq2
DE.05v10.genes <- rownames(DE.homoeolog.devAll[[1]])[DE.homoeolog.devAll[[1]]$padj < 0.05]
DE.05v10.genes <- unique(DE.05v10.genes)
DE.05v10.A.genes <- DE.05v10.genes[grep("_A",DE.05v10.genes)]
DE.05v10.A.genes <- gsub(pattern = "_A",replacement = "",DE.05v10.A.genes)




#  Hierarchical cluster tree (average linkage, dissTOM) of the 37222 proteins. The first color bands provide a simple visual comparison of module assignments (branch cuttings) based on the dynamic hybrid branch cutting method. Other color bands visualize the gene significance measure: "red" indicates a high positive correlation with seed weight, oil content and dpa. Note that the brown, blue and red module contain many genes that have high positive correlations with traits.
#plotColor <- cbind( labels2colors(net$colors), GS.Color, DE_dev.Color, DE_con1.Color, DE_con2.Color)
#plotColor <- cbind( labels2colors(net$colors), GS.Color)
plotColor <- d3[,2:4]
#plotLabel <- c("Aggregate Module","Partitioned A Modules","Partitioned D Modules",paste("GS.dpa",colnames(GS.dpa),sep=""))
plotLabel <- c("Aggregate Module","Partitioned A Modules","Partitioned D Modules")

for(l in 1:length(net$dendrograms))
{
	plotDendroAndColors(net$dendrograms[[1]],colors=plotColor,dendroLabels=FALSE, hang=0.03,addGuide=TRUE,guideHang=0.05, main="Cluster dendrogram with modules")
}

dev.off() 
#I realize this may not make as much sense as I hope (comparing gene expression to dpa)
#GS now signifies expression correlation with dpa, so genes that increases progressively will be highly correlated, that is all
#Probably makes more sense to use the DE than anything else



# Relate eigengenes to external traits or sample conditions
# Add the weight to existing module eigengenes
MET=orderMEs(cbind(MEs,as.numeric(colMat$timepoint)))
#Visualization of the eigengene network representing the relationships among the modules and sample traits. The top panel shows a hierarchical clustering dendrogram of the eigengenes based on the dissimilarity diss(q_1,q_2)=1-cor(E^{(q_1)},E^{(q_2)}). The bottom panel shows the shows the eigengene adjacency A_{q1,q2}=0.5+0.5 cor(E^{(q_1)},E^{(q_2)}).
pdf("s5.relateTraits_eigengeneNetworks.withTX2094.pdf")
plotEigengeneNetworks(MET,"",marDendro=c(0,4,1,2), marHeatmap=c(3,4,1,2),cex.lab=0.8,xLabelsAngle=90)
dev.off()
# plot eigengene network for only the significant modules
pdf("s5.relateTraits_eigengeneNetworksSig.pdf")
module.sig<-rownames(pval[pval$symbol=="*",])
MET.sig<-MET[,module.sig]
plotEigengeneNetworks(MET.sig,"",marDendro=c(0,4,1,2), marHeatmap=c(3,4,1,2),cex.lab=0.8,xLabelsAngle=90)
# plot it again with color names
names(MET.sig)<-paste("ME",labels2colors(as.numeric(substring(names(MET.sig),3) ) ),sep="")
plotEigengeneNetworks(MET.sig,"",marDendro=c(0,4,1,2), marHeatmap=c(3,4,1,2),cex.lab=0.8,xLabelsAngle=90)
dev.off()
#I think the ordering of the MEs is going a little funny, need to check on that

######### I want to plot a color bar for the sig only dendrogram ##
sigModule<-as.numeric(gsub("ME","",rownames(pval[pval$symbol=="*",]) ) )
ch.col <- labels2colors(sigModule)
pdf("s5.sigOnlyDendro.pdf")
plot(1:length(sigModule), 1:length(sigModule), type="n", yaxt="n", ylab="")
for (k in 1:length(sigModule)) {   rect(k, 1, k+1,  2, col = ch.col[k], border=NA) }
dev.off()
#########

pdf("s5.moduleTraitsCorrelation.partitioned.withTX2094.pdf")
# graphical representation for corelation with modules
moduleTraitCor = cor(MEs, colMat[,2], use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples=24);
MEcolors<-paste("ME",labels2colors(as.numeric(gsub("ME","",names(MEs)) ) ), sep="")
# sizeGrWindow(5,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), " (", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mfrow=c(1,1), mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
# Table of module-trait correlations and p-values. Each cell reports the correlation (and p-value) resulting from  correlating module eigengenes (rows) to traits (columns). The table is color-coded by correlation according to the color legend.
labeledHeatmap(Matrix = moduleTraitCor, xLabels = colnames(GS.dpa), yLabels = MEcolors, ySymbols = names(MEs), colorLabels = TRUE, colors = blueWhiteRed(50), textMatrix = as.matrix(textMatrix), setStdMargins = FALSE, cex.text = 0.7,zlim = c(-1,1), main = paste("Module-trait relationships"))
# plot it for only significant modules
where.sig<-sort(match(module.sig, rownames(moduleTraitCor)) )
moduleTraitCor.sig <- moduleTraitCor[where.sig,]
textMatrix.sig <- textMatrix[where.sig,]
labeledHeatmap(Matrix = moduleTraitCor.sig, xLabels = colnames(GS.dpa), yLabels = MEcolors[where.sig], ySymbols = rownames(moduleTraitCor.sig), colorLabels = TRUE, colors = blueWhiteRed(50), textMatrix = as.matrix(textMatrix.sig), setStdMargins = FALSE, cex.text = 0.7,zlim = c(-1,1), main = paste("Module-trait relationships: sig only"))
dev.off()

# For each module, we also define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profile. This allows us to quantify the similarity of all genes on the array to every module.
# calculate the module membership values
# (aka. module eigengene based connectivity kME):
pdf("s5.moduleTraitAndMembershipCor.aggr.withTX2094.pdf")
MM=signedKME(multiExpr[[3]]$data, MEs)
rownames(MM)<-names(multiExpr[[3]]$data)
# Intramodular analysis: identifying genes with high GS and MM
# Using the GS and MM measures, we can identify genes that have a high significance for weight as well as high module membership in interesting modules. As an example, we look at the brown module that a high correlation with body weight. We plot a scatterplot of Gene Significance vs. Module Membership in select modules...
colorOfColumn=substring(names(MM),4)[-1]
par(mfrow = c(2,1))
for (module in sigModule) {
    column = match(module,colorOfColumn)
    restModule=net$colors==module
    verboseScatterplot(MM[restModule,column], GS.dpa[restModule], xlab=paste("Module Membership of ME",module),ylab=paste("GS.dpa"),    main=paste("kME.",module,"vs. GS.dpa"), col=labels2colors(module))
}
dev.off()



####
##Annotation
####

aa<-load('D5annotation.Rdata')#from Phytozome
aa  # "annotation221" "annot"
aa<-annotation221[,c("transcript","tair10.defline", "pfam","panther","kog", "kegg.ec", "kegg.orthology", "go", "tair10", "tair10.symbol")]
dim(aa<-aa[grep("[.]1$",aa$transcript),])  #37505
aaA <- aa; aaD <- aa;
aa$gene = gsub("[.]1$","", aa$transcript)
aaA$gene = paste(gsub("[.]1$","", aaA$transcript),"_A",sep="")
aaD$gene = paste(gsub("[.]1$","", aaD$transcript),"_D",sep="")
aaAll<-rbind(aaA,aaD)
me<-data.frame(gene = colnames(multiExpr[[3]]$data), ME = net$colors)
#I have D and A genes, merge will not work, need something a bit sloppier
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
write.table(MM,file="s5.moduleMembership.partitioned.txt",sep="\t",row.names=FALSE,quote=FALSE)
save(net, colMat, pval, GS.dpa, MM, file = "R-05-all12_GSandMM.RData")


me0<-me
names(MM)[1]<-"gene"
dim(me0<-merge(me0,MM,by="gene",all.x=TRUE, all.y=TRUE) )
#Build table of just hub genes
mefields <- names(me0)[grep(pattern = "kME",x= names(me0))]
bareMe <- gsub(names(me0)[grep(pattern = "kME",x= names(me0))],pattern = "kME",rep = "")
hubgenes <- data.frame()
for(j in 1:length(mefields)){
	temp <- me0[me0[,mefields[j]]>0.9 & me0$ME==bareMe[j], c("gene","ME","tair10.defline","go")]
	hubgenes <- rbind(hubgenes,temp)
}

#get hub genes and order them by MM for their respective modules (prep for GSEA)
for(j in 1:length(mefields)){
	temp <- me0[me0[,mefields[j]]>0.9 & me0$ME==bareMe[j], c("gene",mefields[j])]
	temp <- temp[order(temp[,2],decreasing = TRUE),]
	write.table(temp, file=paste("GSEA/HG/rnk_homoeolog_",mefields[j],".rnk",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}

write.table(me0,file="s5.moduleMembershipAndAnnotation.partitioned.txt",sep="\t",row.names=FALSE)
write.table(hubgenes,file="s5.HubGenesAndAnnotation.partitioned.txt",sep="\t",row.names=FALSE)
write.table(t(multiExpr[[3]]$data),file="s5.rldExpr_all.txt",sep="\t",quote = FALSE)

##Need to write code to count genes per module, average module membership in each module
table(Allnet12$colors)[[1]]
#aggregate
#   0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
#1221 7659 7600 2192 1154 1113  920  916  806  737  533  468  446  412  410  408
#  16   17   18   19   20   21   22   23   24   25   26
# 402  382  375  349  308  270  186  125  124  118   72

#partitioned
#   0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
# 895 9314 8711 3103 1942 1923 1891 1673 1623 1614 1577 1436 1131 1083 1052  971
#  16   17   18   19   20   21   22   23   24   25   26   27   28   29   30   31
# 952  922  757  673  654  606  570  551  491  485  482  382  348  233  231  221
#  32   33   34   35   36   37   38   39   40   41   42   43   44   45   46   47
# 174  173  167  165  160  155  155  153  144  141  136  123  112  106  105   81
#  48   49   50   51   52
#  57   57   49   48   38

#get MM averages
MMavg <- vector()
for(d in 1:length(Allnet12$MEs)){
	dtrue = d - 1
	index <- grep(paste("\\b",dtrue,"\\b",sep = ""),Allnet12$colors)#the \\b marks delimit the word so that 1 does not match 10, 11, etc.
	mod.MM <- MM[index,grep(pattern = paste("\\bkME",dtrue,"\\b",sep=""),x= names(MM))]
	MMavg[d] <- mean(mod.MM)
	names(MMavg)[d] <- paste("kME",dtrue,sep="")
}

#for aggregate
#[1] 0.004739606 #ME0
#[1] 0.6518444
#[1] 0.6449256
#[1] 0.696829
#[1] 0.694955
#[1] 0.6915861
#[1] 0.6636334
#[1] 0.6766454
#[1] 0.7092471
#[1] 0.6924697
#[1] 0.7542956
#[1] 0.6838829
#[1] 0.7158009
#[1] 0.7246418
#[1] 0.7015744
#[1] 0.7740259
#[1] 0.6767776
#[1] 0.7493389
#[1] 0.7556119
#[1] 0.7439153
#[1] 0.7387778
#[1] 0.8213382
#[1] 0.749697
#[1] 0.7432625
#[1] 0.7683519
#[1] 0.7206928
#[1] 0.7979632 #ME26

#for partitioned
#[1] 0.05592105 #ME0
#[1] 0.7159461
#[1] 0.7169094
#[1] 0.6691636
#[1] 0.6346418
#[1] 0.6956112
#[1] 0.6596426
#[1] 0.6854211
#[1] 0.6327866
#[1] 0.6351114
#[1] 0.6270643
#[1] 0.6421879
#[1] 0.6395067
#[1] 0.6309056
#[1] 0.660738
#[1] 0.6643247
#[1] 0.6597614
#[1] 0.6835923
#[1] 0.6911691
#[1] 0.6704279
#[1] 0.6941389
#[1] 0.6763429
#[1] 0.68851
#[1] 0.8365576
#[1] 0.6530131
#[1] 0.6978336
#[1] 0.6091138
#[1] 0.7017622
#[1] 0.7032193
#[1] 0.7458701
#[1] 0.7113125
#[1] 0.6289287
#[1] 0.741264
#[1] 0.7136483
#[1] 0.7243897
#[1] 0.7343503
#[1] 0.6971909
#[1] 0.8781597
#[1] 0.6812239
#[1] 0.6934196
#[1] 0.6619941
#[1] 0.7190801
#[1] 0.8043843
#[1] 0.7025075
#[1] 0.7467695
#[1] 0.7039559
#[1] 0.7334005
#[1] 0.7285506
#[1] 0.746544
#[1] 0.7029991
#[1] 0.8207001
#[1] 0.8123384
#[1] 0.7835598 #ME52

####
##Module conservation between dom and wild
####

#library(WGCNA);
#library(flashClust);
#library(RColorBrewer);
#library(ggplot2);
#options(stringsAsFactors = FALSE);

#remove(list=ls())
#load("R-07-prep.RData") # "multiExpr"          "lanova"             "powerTables.signed"
#load("R-07-individualNetworks.RData") # "A2net12"  "AD3net12" "D5net12"  "Synnet12" "TM1net12" "Yucnet12"
#source('multiscalefreeplot.r', chdir = TRUE)
#source('summarySE.r', chdir = TRUE)
#source('multiplot.r', chdir = TRUE)

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
#dpa24=as.factor(rep(seq(5,20,by=5),6))
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


# for any consensus module analysis, define new CmultiExpr and TOMs
comparisons <- list( all =1:2)          ## all five consensus
)

#Need to edit this block for my set, don't worry about the loop, just need to run through once
print(paste("Consensus module analysis for ",sep=""))
print(shortLabels[1:2])
    
    
######### Start of Jing's consensus analysis unit ########### CANNOT DO THIS FOR WILD AND DOM IF TOO BIG FOR ONE TOM
CmultiExpr = vector(mode = "list", length = nSets)
for(i in 1:nSets){ CmultiExpr[[i]] = multiExpr[[i]]}
fileName<-paste(shortLabels[1:2],sep="",collapse="_")

print("Start to read in TOMs.")
# Initialize an appropriate array to hold the TOMs
TOMs = array(0, dim = c(nSets, nGenes, nGenes));
# load and store TOMs from each individual data set
for (set in 1:nSets)
{
    load(paste(shortLabels[set], "_power12_TOM-block.1.RData", sep="") )   #biocrunch
    TOMs[set, , ]<-as.matrix(TOM)
}


# Step-by-step consensus network construction 
# consensus is defined as the component-wise minimum of the multiple TOMs.

# first, TOMs need to be scaled to be comparable, using the 95th percentile
# BUT BE CAREFUL, for each individual network, original TOMs and scaled TOMs might generate different numbers of modules!!!!
print("Scale TOMs using 95th percentile: values, scaling powers")
scaleP = 0.95
# Set RNG seed for reproducibility of sampling
set.seed(12345)
# Sample sufficiently large number of TOM entries
nSamples = as.integer(1/(1-scaleP) * 1000)
# Choose the sampled TOM entries
scaleSample = sample(nGenes*(nGenes-1)/2, size = nSamples)
TOMScalingSamples = list();
# These are TOM values at reference percentile
scaleQuant = rep(1, nSets)
# Scaling powers to equalize reference TOM values
scalePowers = rep(1, nSets)
# prepare the TOM scaled
TOMscaled <-TOMs
# Loop over sets
for (set in 1:nSets)
{
    # Select the sampled TOM entries
    TOMScalingSamples[[set]] = as.dist(TOMs[set, , ])[scaleSample]
    # Calculate the 95th percentile
    scaleQuant[set] = quantile(TOMScalingSamples[[set]], probs = scaleP, type = 8)
    # Scale the other TOMs
    if (set>1)
    {
        scalePowers[set] = log(scaleQuant[1])/log(scaleQuant[set]);
        TOMscaled[set, ,] = TOMs[set, ,]^scalePowers[set];
    }
}
# check the scaling achieved using a quantile-quantile plot
scaledTOMSamples = list();
for (set in 1:nSets){scaledTOMSamples[[set]] = TOMScalingSamples[[set]]^scalePowers[set]}
# Open a suitably sized graphics window
pdf(paste("s8.checkTOMscaling-",fileName,".pdf",sep=""))
pwset<-combn(nSets,2)
for(i in 1:choose(nSets,2))
{
    # qq plot of the unscaled samples
    qqUnscaled = qqplot(TOMScalingSamples[[pwset[1,i]]], TOMScalingSamples[[pwset[2,i]]], plot.it = TRUE, cex = 0.6, xlab = paste("TOM in", setLabels[pwset[1,i]]), ylab = paste("TOM in", setLabels[pwset[2,i]]), main = "Q-Q plot of TOM", pch = 20)
    # qq plot of the scaled samples
    qqScaled = qqplot(scaledTOMSamples[[pwset[1,i]]], scaledTOMSamples[[pwset[2,i]]], plot.it = FALSE)
    points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20);
    abline(a=0, b=1, col = "blue")
    legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black", "red"))
}
dev.off()
# the red scaled quantile plot should be closer than black unscaled plot to the theoretic blue line
print(scaleQuant)   # [1] 0.2126422 0.1969166
print(scalePowers)  # [1] 1.000000 0.952719

# Second, calculate the consensus Topological Overlap by taking the component-wise ("parallel") minimum of the TOMs in individual sets.
# Between diploids - diploid consensus network vs A2 and D5 each
if(nSets==2) consensusTOM = pmin(TOMscaled[1, , ], TOMscaled[2, , ]) #I know that I have nSets = 2, so this one
#if(nSets==3) consensusTOM = pmin(TOMscaled[1, , ], TOMscaled[2, , ], TOMscaled[3, , ])
#if(nSets==4) consensusTOM = pmin(TOMscaled[1, , ], TOMscaled[2, , ], TOMscaled[3, , ], TOMscaled[4, , ])
#if(nSets==5) consensusTOM = pmin(TOMscaled[1, , ], TOMscaled[2, , ], TOMscaled[3, , ], TOMscaled[4, , ], TOMscaled[5, , ])
#if(nSets==6) consensusTOM = pmin(TOMscaled[1, , ], TOMscaled[2, , ], TOMscaled[3, , ], TOMscaled[4, , ], TOMscaled[5, , ], TOMscaled[6, , ])
# Clustering
consTree = flashClust(as.dist(1-consensusTOM), method = "average")
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
unmergedLabels = cutreeDynamic(dendro = consTree, distM = as.matrix(1-consensusTOM), deepSplit = 2, cutHeight = 0.995, minClusterSize = minModuleSize, pamRespectsDendro = TRUE )
unmergedColors = labels2colors(unmergedLabels)
# a quick summary of the module detection,
print(paste("Before merging, ", length(unique(unmergedLabels)), " modules were found.",sep=""))
#[1] "Before merging, 94 modules were found."
# Calculate module eigengenes
unmergedMEs = multiSetMEs(CmultiExpr, colors = NULL, universalColors = unmergedColors)
# Calculate consensus dissimilarity of consensus module eigengenes
consMEDiss = consensusMEDissimilarity(unmergedMEs);
# Cluster consensus modules
consMETree = flashClust(as.dist(consMEDiss), method = "average");
# merge modules
merge = mergeCloseModules(CmultiExpr, unmergedLabels, cutHeight = 0.25, verbose = 3)
# Numeric module labels
moduleLabels = merge$colors;
# the corresponding colors for large module ID need to be adjusted, otherwise they cannot be plotted
sortedModules = sort(unique(moduleLabels))
moduleLabels.adjusted = match(moduleLabels,sortedModules)-1
# Convert labels to colors
moduleColors = labels2colors(moduleLabels.adjusted)
# Eigengenes of the new merged modules:
consMEs = merge$newMEs;
# Calculate new module eigengenes
mergedMEs = multiSetMEs(CmultiExpr, colors = NULL, universalColors = moduleLabels)
# Calculate consensus dissimilarity of consensus module eigengenes, Cluster consensus modules
consMETree.merged = flashClust(as.dist(consensusMEDissimilarity(mergedMEs) ), method = "average");

# Save useful info for the consensus modules
save(consMEs, unmergedLabels, moduleColors, moduleLabels, moduleLabels.adjusted, consMETree, consMETree.merged, file = paste("s8.moduleConsensus.", fileName, ".RData", sep="") )
print("Save consensus modules information and start to draw figures.")


# Plot the consensus module results
pdf(paste("s8.consensus-",fileName,".pdf",sep=""))
par(mfrow = c(1,1))
plot(consMETree, main = "Consensus clustering of consensus module eigengenes before merging",xlab = "", sub = "")
abline(h=0.25, col = "red")
plot(consMETree.merged, main = "Consensus clustering of consensus module eigengenes after merging",xlab = "", sub = "")
plotDendroAndColors(consTree, cbind(unmergedColors, moduleColors), c("Unmerged", "Merged"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
plotDendroAndColors(consTree, cbind(moduleColors, labels2colors(all.colors) ), c("Consensus",shortLabels[1:2]), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Consensus gene dendrogram and module colors")
# Eigengenes are the 1st principal component of modules in a given single dataset, which provide a summary profile for each module.
# Displaying module heatmap and the eigengene
# sizeGrWindow(8,7);
for(i in 1:2)
{   MEs<-consMEs[[i]]$data
    Nmodules<-dim(MEs)[2]
    module.names<-names(MEs)
    dpa<-dpas[[i]]
    plots <- list()  # new empty list
    for(me in 1:Nmodules)
    {
        which.module=module.names[me]
        module.color=labels2colors(match(as.numeric(substring(which.module,3)),sortedModules)-1)
        #heatmap
        par(mfrow=c(2,1), mar=c(0.3, 5.5, 4, 2))
        plotMat(t(scale(CmultiExpr[[i]]$data[,moduleColors==module.color ]) ), nrgcols=30,rlabels=T,rcols=module.color,
        main=paste(shortLabels[i],which.module, module.color, sep=": "), cex.main=2)
        #barplot
        par(mar=c(5, 4.2, 0, 0.7))
        barplot(MEs[,which.module], col=module.color, main="", cex.main=2,
        ylab="eigengene expression",xlab="seed development (dpa)", names.arg=as.character(dpa) )
        #line, anova
        df<-data.frame(ME=MEs[,me], dpa, module = which.module )
        fit<-aov(ME~dpa,df)
        dfc<-summarySE(df, measurevar="ME", groupvars=c("dpa", "module"))
        plots[[me]]<- ggplot(dfc, aes(x=dpa, y=ME, group=module)) +
        geom_errorbar(aes(ymin=ME-se, ymax=ME+se), width=.1) +
        geom_line(colour=module.color) + geom_point( ) +
        ggtitle(paste(which.module," ",module.color,", anova P=", round(anova(fit)$"Pr(>F)"[1], 4), sep="") )+
        theme(plot.title=element_text( size=8))
    }
    for(page in 1:ceiling(Nmodules/9))
    {
        if(Nmodules>(9*page))
        {  multiplot(plotlist = plots[(9*page-8):(9*page)],  layout=matrix(1:9, nrow=3, byrow=TRUE) )  }
        else
        {  multiplot(plotlist = plots[(9*page-8):Nmodules],  layout=matrix(1:9, nrow=3, byrow=TRUE) )  }
    }
}
dev.off()


# Calcultae and Plot consensus module preservation results
pdf(paste("s8.modulePreservation-",fileName,".pdf",sep=""),width=8, height=10)
# Recalculate consMEs to give them color names, or
# consMEs = multiSetMEs(multiExpr, universalColors = merge$colors);
#sizeGrWindow(8,10);
par(cex = 0.4)  #very import for cex
plotEigengeneNetworks(consMEs, setLabels[1:2], marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1), zlimPreservation = c(0.5, 1), xLabelsAngle = 90)
# Characterizing consensus modules by differential expression of their corresponding eigengenes in the various time points. Red means over-expression, green under-expression; numbers in each cell give the corresponding t-test p-value. Each column corresponds to an eigengene and each row corresponds to a time point.
# setCorrelationPreservation(consMEs, setLabels[comp], excludeGrey = TRUE, greyLabel = "grey")
# pairwise plot
#pwset<-combn(nSets,2)
#par(cex=0.4);
#plotEigengeneNetworks(list(consMEs[[pwset[1,1]]],consMEs[[pwset[2,1]]]), setLabels[c(pwset[1,1],pwset[2,1])], marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1), zlimPreservation = c(0.5, 1), xLabelsAngle = 90)

# with GS and MM
par(mfrow=c(4,1) )
#Need to redo again for differential expression data
for(i in 1:1 )
{
    # the pair for comparison
    MEs1<-consMEs[[pwset[1,i]]]$data
    MEs2<-consMEs[[pwset[2,i]]]$data
    dpa1<-dpas[[pwset[1,i]]]
    dpa2<-dpas[[pwset[2,i]]]
    # diffTable <- apply(MEs1-MEs2,2, function(x){unlist(lapply(split(x,dpa),mean ))} )
    diffTable <- apply(MEs1,2, function(x){unlist(lapply(split(x,dpa1),mean ))} ) - apply(MEs2,2, function(x){unlist(lapply(split(x,dpa2),mean ))} )
    pvalTable <- matrix(0, nrow = 4, ncol = dim(MEs1)[2]);
    for(cc in 1:dim(MEs1)[2])
    {
        for(rr in 1:4)
        { pvalTable[rr,cc] = t.test(split(MEs1[,cc],dpa1)[[rr]],split(MEs2[,cc],dpa2)[[rr]])$p.value }
    }
    labeledHeatmap( Matrix = diffTable, xLabels = names(MEs1), yLabels = c(5,10,15,20),
    colorLabels = TRUE, colors = blueWhiteRed(50),
    textMatrix = round(as.matrix(pvalTable),2),
    cex.text = 0.7,  zlim = c(-0.8,0.8),setStdMargins = FALSE,
    main = paste("Consensus MEs differential expression: ", paste(shortLabels[pwset[,i]],collapse=" vs " ), sep=""  ) )
}
dev.off()


# Note that preservation measures can also be generated through (setCorrelationPreservation)
# do a permutation test for significance
resultP<-data.frame(matrix(ncol=2))
names(resultP)<-c()
pwset<-combn(nSets,2)
for(i in 1:1 )
    {
    names(resultP)[i]<-paste(shortLabels[pwset[1,i]],"vs",shortLabels[pwset[2,i]],sep="")
    }
# make 1000 times permutation
nP = 1000
for(i in 1:nP) {
# permutation
colorsP<-sample(moduleLabels, size=length(moduleLabels), replace=FALSE)
# permutate module color assignment
consMEsP = multiSetMEs(multiExpr, universalColors = colorsP);
# given module colors, calculate conMEs
dd<-setCorrelationPreservation(consMEsP, setLabels, excludeGrey = TRUE, greyLabel = "grey")
resultP[i,]<-as.matrix(dd)[2]
}
# observed top 2.5% values should be larger than
apply(resultP, 2, function(x){sort(x, decreasing=TRUE)[25]})
#DomvsWild
#0.7780247
apply(resultP, 2, function(x){sort(x, decreasing=FALSE)[25]})
#DomvsWild
#0.5389602
#0.5389602

############### actual values are
D<-setCorrelationPreservation(consMEs, setLabels[1:2], excludeGrey = TRUE, greyLabel = "grey")
#aggregate network
#                         Domesticated G. hirsutum Wild G. hirsutum
#Domesticated G. hirsutum                0.0000000        0.7730014
#Wild G. hirsutum                        0.7730014        0.0000000

#Partitioned network
#                         Domesticated G. hirsutum Wild G. hirsutum
#Domesticated G. hirsutum                0.0000000        0.7633359
#Wild G. hirsutum                        0.7633359        0.0000000


save(D, resultP, file = "s8.preservationD.RData" )

# Plot marginal analysis between consensus modules and each individual dataset
pdf(paste("s8.ConsensusAndMarginal-",fileName,".pdf",sep=""),width=10,height=7)
par(mfrow=c(1,1));
par(cex = 1.0);
par(mar=c(8, 10.4, 2.7, 1)+0.3);
consMEs.no<-length(unique(moduleLabels))
consMEs.module<-labels2colors(match(as.numeric(substring(names(consMEs[[1]]$data ),3) ), sortedModules)-1 )  #colors in order
# loop pairwise comparison
for(i in 1:2)  
{
    coln<-consMEs.no                  #put consensus modules in columns
    rown<-ncol(nets[[i]]$MEs )  # put each individual set modules in rows
    # Initialize tables of p-values and of the corresponding counts
    pTable = matrix(0, nrow = rown, ncol = coln);
    CountTbl = matrix(0, nrow = rown, ncol = coln);
    # color list of MEs in the color of decreasing numbers of memebers
    colModule  <- consMEs.module
    rowModule  <- labels2colors(as.numeric(names(table(nets[[i]]$colors)) ))
    # colors for each gene
    colColors  <- moduleColors
    rowColors  <- labels2colors(nets[[i]]$colors )
    # anova significance sybol
    rowP  <- anovaP[[i]]$symbol
    # Initialize tables of p-values and of the corresponding counts
    pTable = matrix(0, nrow = rown, ncol = coln);
    CountTbl = matrix(0, nrow = rown, ncol = coln);
    # Execute all pairwaise comparisons
    for (rmod in 1:rown) {
		for (cmod in 1:coln){
			rMembers = (rowColors == rowModule[rmod] );
			cMembers = (colColors == colModule[cmod] );
			pTable[rmod, cmod] = -log10(fisher.test(rMembers, cMembers, alternative = "greater")$p.value);
			CountTbl[rmod, cmod] = sum(rowColors == rowModule[rmod]  & colColors == colModule[cmod]  )
		}
	}
    
    # display the p-value and counts in a color-coded table. The colors will indicate the p-value signicance
    # Truncate p values smaller than 10^{-50} to 10^{-50}
    pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
    pTable[pTable>50 ] = 50 ;
    # Marginal counts (really module sizes)
    rModTotals = apply(CountTbl, 1, sum)
    cModTotals = apply(CountTbl, 2, sum)
    # Use function labeledHeatmap to produce the color-coded table with all the trimmings
    labeledHeatmap( Matrix = pTable, colorLabels = TRUE,
    xLabels = paste(" ", colModule), yLabels = paste(" ", rowModule),
    xSymbols = paste(names(consMEs[[1]]$data ),"-", colModule, ": ", cModTotals, " ", sep=""),
    ySymbols = paste(names(nets)[i], rowModule, ": ", rModTotals, rowP, sep=""),
    textMatrix = CountTbl, colors = blueWhiteRed(100)[50:100],
    main = "Correspondence of dataset-specific to consensus modules ",
    cex.text = 1, cex.lab = 1, setStdMargins = FALSE      )
}
dev.off()

# Plot marginal analysis between consensus modules and each individual dataset
# comp<-comparisons[[li]]
# nSets<-length(comp)
# fileName<-paste(shortLabels[comp],sep="",collapse="_")
# sortedModules = sort(unique(moduleLabels))
pdf(paste("s8.ConsensusAndMarginal-",fileName,".sig.pdf",sep=""),width=10,height=7)
par(mfrow=c(1,1));
par(cex = 1.0);
par(mar=c(8, 10.4, 2.7, 1)+0.3);
consMEs.no<-length(unique(moduleLabels))
consMEs.module<-labels2colors(match(as.numeric(substring(names(consMEs[[1]]$data ),3) ), sortedModules)-1 )
cMod_colors <- cbind(match(as.numeric(substring(names(consMEs[[1]]$data ),3) ), sortedModules)-1 ,consMEs.module)  #colors in order
# loop pairwise comparison
for(i in 1:2)  
{
    coln<-consMEs.no                  #put consensus modules in columns
    rown<-ncol(nets[[i]]$MEs )  # put each individual set modules in rows
    # Initialize tables of p-values and of the corresponding counts
    pTable = matrix(0, nrow = rown, ncol = coln);
    CountTbl = matrix(0, nrow = rown, ncol = coln);
    # color list of MEs in the color of decreasing numbers of memebers
    colModule  <- consMEs.module
    rowModule  <- labels2colors(as.numeric(names(table(nets[[i]]$colors)) ))
    # colors for each gene
    colColors  <- moduleColors
    rowColors  <- labels2colors(nets[[i]]$colors )
    # anova significance sybol
    rowP  <- anovaP[[i]]$symbol
    # Initialize tables of p-values and of the corresponding counts
    pTable = matrix(0, nrow = rown, ncol = coln);
    CountTbl = matrix(0, nrow = rown, ncol = coln);
    # Execute all pairwaise comparisons
    for (rmod in 1:rown){
    for (cmod in 1:coln)
    {
        rMembers = (rowColors == rowModule[rmod] );
        cMembers = (colColors == colModule[cmod] );
        pTable[rmod, cmod] = -log10(fisher.test(rMembers, cMembers, alternative = "greater")$p.value);
        CountTbl[rmod, cmod] = sum(rowColors == rowModule[rmod]  & colColors == colModule[cmod]  )
    }
    }
    # display the p-value and counts in a color-coded table. The colors will indicate the p-value signicance
    # Truncate p values smaller than 10^{-50} to 10^{-50}
    pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
    pTable[pTable>50 ] = 50 ;
    # Marginal counts (really module sizes)
    rModTotals = apply(CountTbl, 1, sum)
    cModTotals = apply(CountTbl, 2, sum)
    # Use function labeledHeatmap to produce the color-coded table with all the trimmings
    labeledHeatmap( Matrix = pTable[rowP=="*",], colorLabels = TRUE,
    xLabels = paste(" ", colModule), yLabels = paste(" ", rowModule[rowP=="*"]),
    xSymbols = paste(names(consMEs[[1]]$data ),"-", colModule, ": ", cModTotals, " ", sep=""),
    ySymbols = paste(names(nets)[i], rowModule[rowP=="*"], ": ", rModTotals[rowP=="*"], rowP[rowP=="*"], sep=""),
    textMatrix = CountTbl[rowP=="*",], colors = blueWhiteRed(100)[50:100],
    main = "Correspondence of dataset-specific to consensus modules ",
    cex.text = 0.5, cex.lab = 0.5, setStdMargins = FALSE      )
}
dev.off()


######### End of my consensus analysis unit ###########

###JOE's unique partitioned stuff###
#1) I am trying to get look at homoeolog bias within partitioned network
#2) As well as trying to get the differences in modules bt aggr and partitioned
modDat <- as.data.frame(cbind(colnames(multiExpr[[3]]$data),Allnet12$colors))
table(factor(modDat[grep(pattern = "*_A",modDat[,1]),2])) #gets A contribution to each module
#   0    1   10   11   12   13   14   15   16   17   18   19    2   20   21   22
# 421 4645  782  722  555  551  522  435  478  439  360  337 4347  325  299  285
#  23   24   25   26   27   28   29    3   30   31   32   33   34   35   36   37
# 274  274  242  231  276  170  108 1582  129  114   79   85   81   71   83   86
#  38   39    4   40   41   42   43   44   45   46   47   48   49    5   50   51
#  70   85 1059   77   66   69   60   57   50   45   36   31   22  935   17   27
#  52    6    7    8    9
#  21  932  829  788  810
table(factor(modDat[grep(pattern = "*_D",modDat[,1]),2]))
#   0    1   10   11   12   13   14   15   16   17   18   19    2   20   21   22
# 474 4669  795  714  576  532  530  536  474  483  397  336 4364  329  307  285
#  23   24   25   26   27   28   29    3   30   31   32   33   34   35   36   37
# 277  217  243  251  106  178  125 1521  102  107   95   88   86   94   77   69
#  38   39    4   40   41   42   43   44   45   46   47   48   49    5   50   51
#  85   68  883   67   75   67   63   55   56   60   45   26   35  988   32   21
#  52    6    7    8    9
#  17  959  844  835  804
#To confirm the totals
table(factor(modDat[,2]))
#
#   0    1   10   11   12   13   14   15   16   17   18   19    2   20   21   22
# 895 9314 1577 1436 1131 1083 1052  971  952  922  757  673 8711  654  606  570
#  23   24   25   26   27   28   29    3   30   31   32   33   34   35   36   37
# 551  491  485  482  382  348  233 3103  231  221  174  173  167  165  160  155
#  38   39    4   40   41   42   43   44   45   46   47   48   49    5   50   51
# 155  153 1942  144  141  136  123  112  106  105   81   57   57 1923   49   48
#  52    6    7    8    9
#  38 1891 1673 1623 1614

modPairCounts <- as.data.frame(rbind(table(factor(modDat[grep(pattern = "*_A",modDat[,1]),2])),table(factor(modDat[grep(pattern = "*_D",modDat[,1]),2]))),row.names = c("A","D"))
modPairCountsT <- t(modPairCounts)
modADat <- modDat[grep(pattern = "*_A",modDat[,1]),]
modADat <- modDat[grep(pattern = "*_D",modDat[,1]),]
modADat$V1 <- gsub(pattern = "_A", rep = "", x = modADat[,1])
modDDat$V1 <- gsub(pattern = "_D", rep = "", x = modDDat[,1])
modADDat <- merge(modADat,modDDat,all = TRUE)
modADDat$moduleIdentity <- modADDat$A.module == modADDat$D.module
table(modADDat$moduleIdentity)
#
#FALSE  TRUE
#15891  7561
#length(modADDat$moduleIdentity)
#[1] 27544 #lots of NAs too where just one homoeolog is used
#uniqGenePair <- unique(gsub(pattern = "_[AD]",replacement = "",x=colnames(multiExpr[[3]]$data),perl = T))
#> table(is.na(modADDat$moduleIdentity))
#
#FALSE  TRUE
#23452  4092

#Jing suggested doing a modulePreservation() analysis between the At network and Dt network, to see if the
#module structure of the aggregate network is preserved in each of the partitioned networks
#utilizes the module structure from aggregate network and the multiExpr of the At genes

#So I run this in the partitioned section after porting over the colors list
#set up expr and color sets
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

pdf("FIGURE3_ROUGH_medianRankAndZsummaryForAggrVsPart.pdf")
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

pdf("Allnet12_module26_heatmap_and_scatterplot.pdf")
whichmodule="26"
Eigengene26 <- Allnet12$MEs[,colnames(Allnet12$MEs) == "ME26"]
datExprModule=multiExpr[[3]]$data[,moduleColors=="darkorange"]
# set the margins of the graphics window
# par(mar=c(0.3, 5.5, 3, 2))
# create a heatmap whose columns correspond to the libraries
# and whose rows correspond to genes
plot.mat(t(scale(datExprModule)),cex.axis=2,nrgcols=30,rlabels=F,rcols="darkorange",main=paste("heatmap of,"whichmodule,"module"))
#scatter plot between eigengene and sample network connectivity
verboseScatterplot(Eigengene26,Z.k,xlab=paste("ME",whichmodule,sep=""))
abline(h=-2,col="red",lwd=2)
dev.off()

#looking at the scatterplot, this shows that all the TX665 samples have a very different eigengene value for module 26 ("dark orange") in the aggregate Allnet12 network
#need to look at GO terms, TAIR defline or something to find any interesting genes

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
#I forgot to do this so setting it up now
splitExpr <- list(Dom = multiExpr[[1]],Wild = multiExpr[[2]])
splitColor<-list(Dom = Domnet12$colors, Wild = Wildnet12$colors)

pwset<-combn(2,2)

nPermutations1=200

set.seed(1)
system.time({
    mp = modulePreservation(splitExpr, splitColor, networkType="signed", referenceNetworks = c(1,2), testNetworks= list(2,1), nPermutations = nPermutations1,
    randomSeed = 1, quickCor = 0, verbose = 3)
})

pdf("FIGURE4_modPresGraphsForWildVsDomPart.pdf")
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

#Do the dom to wild comparison again with consensus moduleSize
splitExpr <- list(Dom = multiExpr[[1]],Wild = multiExpr[[2]])
splitColor2<-list(Dom = cModColors, Wild = cModColors)

nPermutations1=200

set.seed(1)
system.time({
    mp_cMod = modulePreservation(splitExpr, splitColor2, networkType="signed", referenceNetworks = list(1,2), testNetworks= list(2,1), nPermutations = nPermutations1,
    randomSeed = 1, quickCor = 0, verbose = 3)
})

pdf("medianRankAndZsummaryForWildVsDomAggrwCModules.pdf")
par(mfrow=c(2,2),mar = c(4.5,4.5,2.5,1))
comp <- list(c(1,2),c(2,1))
for(co in comp){
	Obs.PreservationStats= mp_cMod$preservation$observed[[co[1]]][[co[2]]]
	Z.PreservationStats=mp_cMod$preservation$Z[[co[1]]][[co[2]]]
	# Let us now visualize the data.
	modColors = rownames(Obs.PreservationStats)
	moduleSize = Obs.PreservationStats$moduleSize
	# we will omit the grey module (background genes)
	# and the gold module (random sample of genes)
	selectModules = !(modColors %in% c("grey", "gold"))
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

#I would like to look at differential expression by module and see if there is correlation between module membership and differential expression
#1. Assign significant DE to all genes, and log fold change
#2. Look at simple presence absence of DE in modules (both aggr and part) .. CHECK
#3. Look at correlation between log fold change and module membership

#DE.aggr.devByCon contains all comparisons between timepoints within condition and between conditions within a timepoint
#So we want to get those that have a significant value (padj < 0.05) and look at which modules they are a part of
#In partitioned network, replace DE.aggr.devByCon with DE.homoeolog.devByCond
gnames <- list()
for(i in 1:length(DE.aggr.devByCon)){
	index <- DE.aggr.devByCon[[i]]$padj < 0.05; index[is.na(index)] <- FALSE;
	compnames <- rownames(DE.aggr.devByCon[[i]][index,])
	gnames[[i]] <- compnames
}
gnames <- unique(sort(unlist(gnames)))
length(gnames)
#[1] 10279 in aggr
#[1] 14669 in partitioned

DE.aggrByMod <- Allnet12Modules[Allnet12Modules$Gene %in% gnames,]
#In partitioned, DE.homoeoByMod <- Allnet12.homoeo.modules[Allnet12.homoeo.modules$Gene %in% gnames,]
table(DE.aggrByMod$Modules)
#compare these numbers totals
#   0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
#  39 3332 2635  943  485  116  204  191  291  154  206  125  167  287  206  124
#  16   17   18   19   20   21   22   23   24   25   26
#  71  117  199   78  150   77   14    7    1   57    3
#In partitioned
#  0    1   10   11   12   13   14   15   16   17   18   19    2   20   21   22
#  21 4080  428  273  378  260  368   57  161  486  117   42 3201   82  226   90
#  23   24   25   26   27   28   29    3   30   31   32   33   34   35   36   37
# 119   11  183   28   75   91    8  734   12   21   36   82   26    6   13   53
#  38   39    4   40   41   42   43   44   45   46   47   48   49    5   51   52
#   2   12  196    3   39    2   11   30    3   39    3    1    1  418   10   21
#   6    7    8    9
# 528  213  903  466

table(Allnet12Modules$Modules)
#
#   0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
#1221 7659 7600 2192 1154 1113  920  916  806  737  533  468  446  412  410  408
#  16   17   18   19   20   21   22   23   24   25   26
# 402  382  375  349  308  270  186  125  124  118   72
#In partitioned
#   0    1   10   11   12   13   14   15   16   17   18   19    2   20   21   22
# 895 9314 1577 1436 1131 1083 1052  971  952  922  757  673 8711  654  606  570
#  23   24   25   26   27   28   29    3   30   31   32   33   34   35   36   37
# 551  491  485  482  382  348  233 3103  231  221  174  173  167  165  160  155
#  38   39    4   40   41   42   43   44   45   46   47   48   49    5   50   51
# 155  153 1942  144  141  136  123  112  106  105   81   57   57 1923   49   48
#  52    6    7    8    9
#  38 1891 1673 1623 1614


table(DE.aggrByMod$Modules)/table(Allnet12Modules$Modules)
#ratio of genes in each module that are DE in at least one comparison
#          0           1           2           3           4           5
#0.031941032 0.435043739 0.346710526 0.430200730 0.420277296 0.104222821
#          6           7           8           9          10          11
#0.221739130 0.208515284 0.361042184 0.208955224 0.386491557 0.267094017
#         12          13          14          15          16          17
#0.374439462 0.696601942 0.502439024 0.303921569 0.176616915 0.306282723
#         18          19          20          21          22          23
#0.530666667 0.223495702 0.487012987 0.285185185 0.075268817 0.056000000
#         24          25          26
#0.008064516 0.483050847 0.041666667

#Because one module doesn't have any DEGs in it, the arrays don't line up
#Will need some extra code
#FOR PARTITIONED NETWORK START
temp1 <- data.frame(table(Allnet12.homoeo.modules$Modules))
temp2 <- data.frame(table(DE.homoeoByMod$Modules))
names(temp1) <- c("Modules","TotalFreq")
names(temp2) <- c("Modules","DEFreq")
ModuleNumbersAndDE <- merge(x=temp1,y=temp2, all.x = TRUE)
#END
#Make sure to turn NAs into 0 in the DE column
#Then create new column
ModuleNumbersAndDE$Ratio <- ModuleNumbersAndDE$DEFreq/ModuleNumbersAndDE$TotalFreq
#        [,1]         [,2]         [,3]         [,4]         [,5]
#Modules "0"          "1"          "10"         "11"         "12"
#Ratio   "0.02346369" "0.43805025" "0.27140140" "0.19011142" "0.33421751"
#        [,6]         [,7]         [,8]         [,9]         [,10]
#Modules "13"         "14"         "15"         "16"         "17"
#Ratio   "0.24007387" "0.34980989" "0.05870237" "0.16911765" "0.52711497"
#        [,11]        [,12]        [,13]        [,14]        [,15]
#Modules "18"         "19"         "2"          "20"         "21"
#Ratio   "0.15455746" "0.06240713" "0.36746642" "0.12538226" "0.37293729"
#        [,16]        [,17]        [,18]        [,19]        [,20]
#Modules "22"         "23"         "24"         "25"         "26"
#Ratio   "0.15789474" "0.21597096" "0.02240326" "0.37731959" "0.05809129"
#        [,21]        [,22]        [,23]        [,24]        [,25]
#Modules "27"         "28"         "29"         "3"          "30"
#Ratio   "0.19633508" "0.26149425" "0.03433476" "0.23654528" "0.05194805"
#        [,26]        [,27]        [,28]        [,29]        [,30]
#Modules "31"         "32"         "33"         "34"         "35"
#Ratio   "0.09502262" "0.20689655" "0.47398844" "0.15568862" "0.03636364"
#        [,31]        [,32]        [,33]        [,34]        [,35]
#Modules "36"         "37"         "38"         "39"         "4"
#Ratio   "0.08125000" "0.34193548" "0.01290323" "0.07843137" "0.10092688"
#        [,36]        [,37]        [,38]        [,39]        [,40]
#Modules "40"         "41"         "42"         "43"         "44"
#Ratio   "0.02083333" "0.27659574" "0.01470588" "0.08943089" "0.26785714"
#        [,41]        [,42]        [,43]        [,44]        [,45]
#Modules "45"         "46"         "47"         "48"         "49"
#Ratio   "0.02830189" "0.37142857" "0.03703704" "0.01754386" "0.01754386"
#        [,46]        [,47]        [,48]        [,49]        [,50]
#Modules "5"          "50"         "51"         "52"         "6"
#Ratio   "0.21736869" "0.00000000" "0.20833333" "0.55263158" "0.27921735"
#        [,51]        [,52]        [,53]
#Modules "7"          "8"          "9"
#Ratio   "0.12731620" "0.55637708" "0.28872367"

#Get these split up by comparison
gnames <- list()
for(i in 1:length(DE.aggr.devByCon)){
	index <- DE.aggr.devByCon[[i]]$padj < 0.05; index[is.na(index)] <- FALSE;
	compnames <- rownames(DE.aggr.devByCon[[i]][index,])
	gnames[[i]] <- compnames
}
DE.byConByMod <- list()
for(j in 1:length(gnames)){
	temp <- Allnet12Modules[Allnet12Modules$Gene %in% gnames[[j]],]
	DE.byConByMod[[j]] <- data.frame(table(temp$Modules))
	names(DE.byConByMod[[j]]) <- c("Gene","Freq")
}

for(j in 1:length(gnames)){
	temp <- Allnet12.homoeo.modules[Allnet12.homoeo.modules$Gene %in% gnames[[j]],]
	DE.byConByMod[[j]] <- data.frame(table(temp$Modules))
	names(DE.byConByMod[[j]]) <- c("Gene","Freq")
}

#Need to build into a table

##Network construction and consensuns modules part 2 - multiple TOM blocks
#From https://github.com/huguanjing/AD1_RNA-seq/blob/master/AD1_Networks120116.r
############### Step 4.  Network Construction  ###############
## nohup R CMD BATCH s4.networkConstruction.R &
################

library(WGCNA)
options(stringsAsFactors = FALSE);
enableWGCNAThreads(nThreads=4)
# Allowing parallel execution with up to 4 working processes.
library(RColorBrewer)
library(ggplot2);
sessionInfo()

ptm <- proc.time()

load("R-03-chooseSoftThreshold.Rdata")
Powers

rdatafiles<-c("R-03-dataInput.RData")

# rdatafiles<-c("R-02-dataInput.rsem_rpkm.RData", "R-02-dataInput.polycat_rpkm.RData")

for(file in rdatafiles)
{
    print(file)
    ll<-load(file)
    # compr<-gsub("R-02-dataInput.|.RData","",file)
    # powerEach = Powers[compr]
    # print(paste0(file,", construct network with b = ",powerEach))
    print(shortLabels)
    print(checkSets(multiExpr)$nGenes)
    compr <- c("WildvDom")
	#in this case, multiExpr only contained [[1]] and [[2]], dom and wild
    
    ###### calculate individual TOMs
    print("###### Calculate individual TOMs:")
    iTOMs = blockwiseIndividualTOMs(
    # Input data
    multiExpr,
    # Data checking options
    checkMissingData = TRUE,
    # Options for splitting data into blocks
    blocks = NULL,
    randomSeed = 12345,
    maxBlockSize =  20000,  # as individual network # 5000 for 4G memory, 20000 for 16G, 30000 for 32 G
    # Network construction arguments: correlation options, use bicor instead of default pearson
    corType = "pearson",
    # Adjacency and topology overlap function options
    #power = powerEach,
    power=max(powers),
    networkType = "signed", TOMType = "signed",
    # Save individual TOMs?
    saveTOMs = TRUE,
    individualTOMFileNames = paste0(compr,".iTOM-%N-block.%b.RData")  )
    
    ###### calculate consensus modules
    print("###### Construct consensus networks:")
    cnet = blockwiseConsensusModules(
    # Input data
    multiExpr,
    # Data checking options
    checkMissingData = TRUE,
    # Options for splitting data into blocks
    blocks = NULL,
    randomSeed = 12345,
    maxBlockSize =  20000,  # as individual network # 5000 for 4G memory, 20000 for 16G, 30000 for 32 G
    # Network construction arguments: correlation options, use bicor instead of default pearson
    corType = "pearson",
    # Adjacency and topology overlap function options
    #power = powerEach,
    power=max(powers),
    networkType = "signed", TOMType = "signed",
    # load previous TOMs
    individualTOMInfo = iTOMs,
    # Saving the consensus TOM
    saveConsensusTOMs = TRUE,
    consensusTOMFileNames = paste0(compr,".cTOM-Block%b.RData"),
    # Basic tree cut options
    deepSplit = 2,  #default, known to reasonable
    minModuleSize = 100, #default 20, use 30 for transcriptome, or ncol(subDat)/2
    pamStage = TRUE, pamRespectsDendro = TRUE, #default, known to reasonable
    # Thredhold to merge modules: a height cut of 0.25 corresponding to correlation of 0.75
    mergeCutHeight = 0.25,
    # others
    reassignThreshold = 0,
    numericLabels = TRUE,
    verbose = 3)
    
    ###### calculate individual modules
    print("###### Construct individual networks:")
    # stupid blockwiseModules only load TOM rdata file with "TOM", not "tomDS"
    tomFiles<-grep(paste0(compr,".iTOM"),list.files(), value=TRUE)
    for(fl in tomFiles)
    {
        load(fl)
        TOM<-tomDS
        save(TOM, file=fl)
    }
    rm(TOM)
    collectGarbage()
    for(i in 1:nSets)
    {
        inet = blockwiseModules(
        # Input data
        multiExpr[[i]]$data,
        # Data checking options
        checkMissingData = TRUE,
        # Options for splitting data into blocks
        blocks =  iTOMs$blocks,
        #randomSeed = 12345,
        #maxBlockSize =  500,  # as individual network # 5000 for 4G memory, 20000 for 16G, 30000 for 32 G
        # Network construction arguments: correlation options, use bicor instead of default pearson
        corType = "pearson",
        # Adjacency and topology overlap function options
        #power = powerEach,
        power=max(powers),
        networkType = "signed", TOMType = "signed",
        # load previous TOMs
        loadTOM = TRUE,
        saveTOMFileBase = paste0(compr,".iTOM-",i),
        # Basic tree cut options
        deepSplit = 2,  #default, known to reasonable
        minModuleSize = 100, #default 20, use 30 for transcriptome, or ncol(subDat)/2
        pamStage = TRUE, pamRespectsDendro = TRUE, #default, known to reasonable
        # Thredhold to merge modules: a height cut of 0.25 corresponding to correlation of 0.75
        mergeCutHeight = 0.25,
        # others
        reassignThreshold = 0,
        numericLabels = TRUE,
        verbose = 3)
        
        assign(paste0("inet",i), inet)
        
    }
    
    save(list=c( "iTOMs", "cnet", grep("inet.",ls(), value=TRUE)), file = paste0("R-04-buildNetwork.", compr,".RData"))
    collectGarbage()

}

consMEs = cnet$multiMEs;
moduleLabels = cnet$colors;
# Convert the numeric labels to color labels
moduleColors = labels2colors(moduleLabels)
consTree = cnet$dendrograms[[1]];

proc.time() - ptm

# With maxBlockSize = 20000, too two days total for one comparison
#Module preservation between partitoned wild and dom networks
splitExpr <- list(Dom = multiExpr[[1]],Wild = multiExpr[[2]])
splitColor<-list(Dom = cnet$colors, Wild = cnet$colors)

nPermutations1=200

set.seed(1)
system.time({
    mp_cMod = modulePreservation(splitExpr, splitColor, networkType="signed", referenceNetworks = list(1,2), testNetworks= list(2,1), nPermutations = nPermutations1,
    randomSeed = 1, quickCor = 0, verbose = 3)
})

pdf("medianRankAndZsummaryForWildVsDomPartwCModules.pdf")
par(mfrow=c(2,2),mar = c(4.5,4.5,2.5,1))
comp <- list(c(1,2),c(2,1))
for(co in comp){
	Obs.PreservationStats= mp_cMod$preservation$observed[[co[1]]][[co[2]]]
	Z.PreservationStats=mp_cMod$preservation$Z[[co[1]]][[co[2]]]
	# Let us now visualize the data.
	modColors = rownames(Obs.PreservationStats)
	moduleSize = Obs.PreservationStats$moduleSize
	# we will omit the grey module (background genes)
	# and the gold module (random sample of genes)
	selectModules = !(modColors %in% c("grey", "gold"))
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

#number of DE genes in wild development
length(unique(c(DEgenes[[1]],DEgenes[[2]],DEgenes[[3]])))
#number of DE genes in dom development
length(unique(c(DEgenes[[4]],DEgenes[[5]],DEgenes[[6]])))
#number of DE genes between wild and dom for all timepoints
length(unique(c(DEgenes[[7]],DEgenes[[8]],DEgenes[[9]],DEgenes[[10]])))

#Get hub genes for overall network
All.scaleK <- degree[,3]/max(degree[,3])#degree comes from s4
All.scaleK.df <- data.frame(Gene = colnames(multiExpr[[3]]$data),Scaled.K = All.scaleK)
All.scaleK.df <- All.scaleK.df[order(All.scaleK.df$Scaled.K,decreasing = TRUE),]
write.table(All.scaleK.df,row.names = FALSE,col.names = FALSE, quote = FALSE, sep = "\t", file = "Allnet12.part.hubgenes.rnk")

############### Step 6. export networks for cytoscape  ###############
#nohup R CMD BATCH s6......R &
################
library(WGCNA);
library(flashClust);
library(RColorBrewer);
library(ggplot2);
options(stringsAsFactors = FALSE);


# Extract only FA genes
CytoskelCellWall<-read.table("yang_D5_from_daojun.sorted.txt",header=FALSE)
CCWdiffCor <-read.table("CCWgenes.diffCor.txt",sep = "\t")
dim(CytoskelCellWall) #703
CSCWGenes <-unique(sort(CytoskelCellWall[,1])) # or CCWdiffCor[,2] Gorais
probes = colnames(multiExpr[[3]]$data)   #29706
# locate FA genes
asCSCWs = is.finite(match(probes,CSCWGenes))
asCSCWs2 = is.finite(match(CSCWGenes, probes))
table(asCSCWs)
# FALSE  TRUE
# 29443   263
#No Gorai.N... genes
dim(CCWdiffCor.noAT <- CCWdiffCor.noAT[!duplicated(CCWdiffCor.noAT$Gr),])

# locate CSCWs in modules
moduleCSCWs<-as.data.frame(table((moduleLabels.adjusted[asCSCWs]))
#   Var1 Freq
#1     0   10
#2     1  137
#3     2   71
#4     3   33
#5     4    8
#6     5    5
#7     6    8
#8     7    4
#9     8    8
#10    9    6
#11   10    9
#12   11    2
#13   12    4
#14   13    4
#15   14    7
#16   16    4
#17   17    4
#18   18    9
#19   19    4
#20   20    2
#21   21    3
#22   22    2
#23   23    2
#24   24    1
#25   25    2
#26   26    1

names(moduleCSCWs)<-c("moduleLabels","CSCWs")
moduleAll <-as.data.frame(table(moduleLabels.adjusted))
names(moduleAll)<-c("moduleLabels","All")
moduleCSCWs<-merge(moduleCSCWs, moduleAll, by="moduleLabels" ,all.y=TRUE)
moduleCSCWs$CSCWs[is.na(moduleCSCWs$CSCWs)]<-0
moduleCSCWs <- moduleCSCWs[order(moduleCSCWs$moduleLabels),]
## calculate enrichment
tt<-colSums(moduleCSCWs[,2:3])
moduleCSCWs$fisherP<-round( apply(moduleCSCWs[,2:3],1,
function(x)fisher.test(matrix(as.numeric(c(x,tt-x)), nrow = 2, dimnames = list( c( "CSCWs","all"),c("inModule", "out"))) ,alternative="greater" )$p.value) ,3)
# FAs enriched in below modules
moduleCSCWs[moduleCSCWs$fisherP<0.05,]
#   moduleLabels CSCWs  All fisherP
#1             1   105 6729   0.000
#18           23    15  721   0.002
#21           27     5  144   0.010


# list modules for each CSCW gene
modules<-data.frame(nodeName=probes, module=moduleLabels.adjusted)
CSCWs<-merge(data.frame(nodeName = CSCWGenes), modules, by = "nodeName")
dim(CSCWs<-CSCWs[!is.na(CSCWs$module),])
#[1] 263 2

####
##FROM HERE ON IS ABOUT GETTING FISHER'S TEST RESULTS, GO FURTHER DOWN TO GET BACK TO THE CYTOSCAPE EXPORT
#
#unique(CSCWs$Biological.process)#Do not have Biological.process at the ready
#
#xtabs(~Biological.process+module,data=FAs)
## now I want to do fishers test to show which module is each category enriched at
#pdf("s6.FAs2modules.pdf",width=10,height=7)
#par(mfrow=c(1,1));
#par(cex = 1.0);
#par(mar=c(8, 10.4, 2.7, 1)+0.3);
#coln<-dim(pval)[1]                             # put modules in columns
#rown<-length(unique(FAs$Biological.process) )  # put FAs categories in rows
## Initialize tables of p-values and of the corresponding counts
#pTable = matrix(NA, nrow = rown, ncol = coln);
#CountTbl = matrix(0, nrow = rown, ncol = coln);
## color list of MEs in the color of decreasing numbers of memebers
#colModule  <- as.numeric(gsub("ME","",rownames(pval)))
## category
#rowCategory <- unique(FAs$Biological.process)
## colors for each gene
#colColors  <- labels2colors(1:coln)
## anova significance sybol
#colP  <- pval$symbol
## Execute all pairwaise comparisons
#for (rmod in 1:rown)
#    for (cmod in 1:coln)
#    {
#        rMembers = (FAs$Biological.process == rowCategory[rmod] );
#        cMembers = (FAs$module == colModule[cmod] );
#        CountTbl[rmod, cmod] = sum(FAs$Biological.process == rowCategory[rmod]  & FAs$module == colModule[cmod]   )
#        if(CountTbl[rmod, cmod]>0) pTable[rmod, cmod] = -log10(fisher.test(rMembers, cMembers, alternative = "greater")$p.value);
#    }
## display the p-value and counts in a color-coded table. The colors will indicate the p-value signicance
#    # Truncate p values smaller than 10^{-50} to 10^{-50}
#pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
#pTable[pTable>50 ] = 50 ;
## Marginal counts (really module sizes)
#rModTotals = apply(CountTbl, 1, sum)
#cModTotals = apply(CountTbl, 2, sum)
#select<-(cModTotals>0)
## Use function labeledHeatmap to produce the color-coded table with all the trimmings
#labeledHeatmap( Matrix = pTable[,select], colorLabels = TRUE,
#    xLabels = paste(colP[select],rownames(pval)[select]), yLabels = paste(" ", rowCategory),
#    textMatrix = CountTbl[,select], colors = blueWhiteRed(100)[50:100],
#    main = "Correspondence of FA-related gene categories to modules ",
#    cex.text = 0.5, cex.lab = 0.5, setStdMargins = FALSE      )
#dev.off()
#
## check some FA families
#xtabs(~Protein.Gene.Abbreviation+module, FAs[FAs$Biological.process=="Lipid Transfer Proteins",])
#
#xtabs(~Protein.Gene.Abbreviation+module, FAs[FAs$Biological.process=="Transcription Factors Associated with Lipid Synthesis",])
#
####
##TUNE BACK IN FOR CSCWs
# Export CSCW network to cytoscape
# Cytoscape [2] allows the user to input an edge file and a node file, allowing the user to specify for example the link weights and the node colors.
## R
# Get topological overlap, "TOM"
load("Dom_power12_TOM-block.1.RData")   #or wild
# Select the corresponding Topological Overlap
subTOM = as.matrix(TOM)[asCSCWs, asCSCWs];
str(subTOM)
subProbes = probes[asCSCWs];
dimnames(subTOM) = list(subProbes, subProbes)
subColors<-moduleLabels.adjusted[asCSCWs]
quantile(TOM, 1:10/10)
#wild TOM         10%          20%          30%          40%          50%          60%
#			8.691218e-05 4.596893e-04 1.421986e-03 3.474652e-03 7.505313e-03 1.514111e-02
#				  70%          80%          90%         100%
#			2.971458e-02 5.906640e-02 1.249964e-01 4.683024e-01
#Dom TOM
#         10%          20%          30%          40%          50%          60%
#6.152193e-05 3.541641e-04 1.271477e-03 3.554854e-03 8.510618e-03 1.837867e-02
#         70%          80%          90%         100%
#3.712910e-02 7.251904e-02 1.422749e-01 4.806390e-01

aa<-read.table("s5.moduleAndAnnotation.txt", header=TRUE,sep="\t")
rownames(aa)<-aa$gene
aa_CSCWs<-aa[subProbes,]
aa_CSCWs<-aa_CSCWs[,-1]
table(rownames(aa_CSCWs) %in% CSCWs$Gr )#All of them match
nn<-cbind(CSCWs,aa_CSCWs)
nn$module.color<-labels2colors(nn$module)
library(gplots)
nn$module.hex<-col2hex(nn$module.color)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(subTOM,
    edgeFile = "CytoscapeInput-edges-CSCWs.dom.txt",#or wild if that's what you use
    nodeFile = "CytoscapeInput-nodes-CSCWs.dom.txt",
    weighted = TRUE,
    threshold = 0.02,   #
    nodeNames = rownames(nn),
    nodeAttr = nn )
	
#Already added the CSCW genes, so now go get cModules
# plot eigengene network
# plotEigengeneNetworks(consMEs[[1]]$data,setLabels[1])
for(i in 1:2){
    A = adjacency(consMEs[[i]]$data, power = 1, type="signed")
    # below code does the same to plot the eigengene network dendrogram
    plot(flashClust(as.dist(1-A),method="average") )
    #abline(h=0.2, col = "red")
    # export A to cytoscape for netwrok visualization
    cyt = exportNetworkToCytoscape(A,
    edgeFile = paste("cytoscape_consensus-consMEs",shortLabels[i],"-edge.txt", sep=""),
    nodeFile = paste("cytoscape_consensus-consMEs",shortLabels[i],"-nodes.txt", sep=""),
    weighted = TRUE,
    threshold = 0,  #if higher than minValue, losing nodes
    nodeNames = colnames(consMEs[[i]]$data),
    altNodeNames = NA,
    nodeAttr = NA)
}
##Get interesting stats for these genes and significant modules 
##for dom and wild
A = adjacency(multiExpr[[1]]$data[,asCSCWs],type = "signed")
fNC_dom = fundamentalNetworkConcepts(A)
A = adjacency(multiExpr[[2]]$data[,asCSCWs],type = "signed")
fNC_wild = fundamentalNetworkConcepts(A)

fNC_Dom.sigMods = list()
fNC_Wild.sigMods = list()
sigMods = c("1","23","27")
for(i in 1:2){
	for(mod in 1:length(sigMods)){
		sigModuleGenes <- is.finite(match(moduleLabels.adjusted,sigMods[mod]))
		A = adjacency(multiExpr[[i]]$data[,sigModuleGenes],type = "signed")
		ifelse(i == 1,fNC_Dom.sigMods[[mod]] <- fundamentalNetworkConcepts(A),fNC_Wild.sigMods[[mod]] <- fundamentalNetworkConcepts(A))
	}
}

save(fNC_dom,fNC_wild,fNC_Dom.sigMods,fNC_Wild.sigMods,file = "CSCWnetworkConcepts.RData")


i=1
x<-read.table(file=paste("cytoscape_consensus-consMEs",shortLabels[i],"-edge.txt", sep=""),header=TRUE)
x$name<-paste(x$fromNode,x$toNode,sep=".")
x<-x[,c(7,3)]
names(x)[2]<-shortLabels[i]
edge<-x
for(i in 2)
{
    x<-read.table(file=paste("cytoscape_consensus-consMEs",shortLabels[i],"-edge.txt", sep=""),header=TRUE)
    x$name<-paste(x$fromNode,x$toNode,sep=".")
    x<-x[,c(7,3)]
    names(x)[2]<-shortLabels[i]
    edge<-merge(edge,x, by="name",all.x=TRUE,all.y=TRUE)
}
edge$fromNode <-gsub("[.].*","",edge$name)
edge$toNode <-gsub(".*[.]","",edge$name)
edge$direction<-"undirected"
edge<-edge[,c("fromNode","direction","toNode","Dom","Wild")]
write.table(edge, "cytoscape_consensus-consMEs-edges.txt",sep="\t",quote=FALSE,row.names=FALSE)

#Checking on overlap between homoeologs and aggregate expression
#in partitioned
DE.all.nodiff <- list()
DE.all.nodiff[[1]] <- unique(sort(gsub(pattern = "_[AD]", rep="", unique(c(DEgenes[[1]],DEgenes[[2]],DEgenes[[3]])))))
DE.all.nodiff[[2]] <- unique(sort(gsub(pattern = "_[AD]", rep="", unique(c(DEgenes[[4]],DEgenes[[5]],DEgenes[[6]])))))
DE.all.nodiff[[3]] <- unique(sort(gsub(pattern = "_[AD]", rep="", unique(c(DEgenes[[7]],DEgenes[[8]],DEgenes[[9]],DEgenes[[10]])))))
DE.homoeo.nodiff <- DE.all.nodiff
save(DE.homoeo.nodiff,file = "DE.homoeo.nodiff.RData")
#bring this over to Aggregate
load("DE.homoeo.nodiff")
#Wild share
table(!is.na(match(unique(c(DEgenes[[1]],DEgenes[[2]],DEgenes[[3]])),DE.homoeo.nodiff[[1]])))
#FALSE  TRUE
# 1168  4707

table(!is.na(match(DE.homoeo.nodiff[[1]],unique(c(DEgenes[[1]],DEgenes[[2]],DEgenes[[3]])))))
#FALSE  TRUE
# 1415  4707

#Dom share
table(!is.na(match(unique(c(DEgenes[[4]],DEgenes[[5]],DEgenes[[6]])),DE.homoeo.nodiff[[2]])))
#FALSE  TRUE
# 1211  5978
table(!is.na(match(DE.homoeo.nodiff[[2]],unique(c(DEgenes[[4]],DEgenes[[5]],DEgenes[[6]])))))
#FALSE  TRUE
# 1731  5978
#Between Wild and Dom
table(!is.na(match(unique(c(DEgenes[[7]],DEgenes[[8]],DEgenes[[9]],DEgenes[[10]])),DE.homoeo.nodiff[[3]])))
#FALSE  TRUE
#  632  1558
table(!is.na(match(DE.homoeo.nodiff[[3]],unique(c(DEgenes[[7]],DEgenes[[8]],DEgenes[[9]],DEgenes[[10]])))))
#FALSE  TRUE
# 1363  1558


AllnetColorDat <- aggrColorDat
partColorDat <- data.frame(Gene = colnames(multiExprAll[[3]]$data),Part.Mod = Allnet12$colors)
partColorDat.A <- partColorDat[grep("_A",partColorDat$Gene),]
partColorDat.D <- partColorDat[grep("_D",partColorDat$Gene),]
colnames(partColorDat.A)[2] <- "Part.Mod.A"
colnames(partColorDat.D)[2] <- "Part.Mod.D"
partColorDat.A$Gene <- gsub("_A","",partColorDat.A$Gene)
partColorDat.D$Gene <- gsub("_D","",partColorDat.D$Gene)
partColorDat.AD <- merge(partColorDat.A, partColorDat.D, by = "Gene", all = TRUE)
AllnetColorDat <- merge(AllnetColorDat,partColorDat.AD, all = TRUE)

rown = length(unique(AllnetColorDat$Aggr))
coln = length(unique(c(AllnetColorDat$Part.Mod.A,AllnetColorDat$Part.Mod.D)))
rMods <- sort(as.numeric(unique(AllnetColorDat$Aggr)))
cMods <- sort(as.numeric(unique(c(AllnetColorDat$Part.Mod.A,AllnetColorDat$Part.Mod.D))))
AllnetOverlapColorMat.A <- matrix(nrow = rown, ncol = coln)
for (rmod in 1:rown) {
		for (cmod in 1:coln){
			AllnetOverlapColorMat.A[rmod, cmod] = sum(AllnetColorDat$Aggr[!is.na(AllnetColorDat$Part.Mod.A)] == rMods[rmod]  & AllnetColorDat$Part.Mod.A[!is.na(AllnetColorDat$Part.Mod.A)] == cMods[cmod] )
		}
}

AllnetOverlapColorMat.D <- matrix(nrow = rown, ncol = coln)
for (rmod in 1:rown) {
		for (cmod in 1:coln){
			AllnetOverlapColorMat.D[rmod, cmod] = sum(AllnetColorDat$Aggr[!is.na(AllnetColorDat$Part.Mod.D)] == rMods[rmod]  & AllnetColorDat$Part.Mod.D[!is.na(AllnetColorDat$Part.Mod.D)] == cMods[cmod] )
		}
}

write.table(AllnetOverlapColorMat.A,file = "AggrVsPartAModulesCrossTable.txt",sep = "\t", quote = FALSE, row.names = FALSE)
write.table(AllnetOverlapColorMat.D,file = "AggrVsPartDModulesCrossTable.txt",sep = "\t", quote = FALSE, row.names = FALSE)


pdf("crosstabAgenomeAndAggr.pdf")
#make heatmap of A vs Aggr
rModTotals = apply(AllnetOverlapColorMat.A, 1, sum)
cModTotals = apply(AllnetOverlapColorMat.A, 2, sum)

#Doesn't produce a great looking result bc of large diff in numbers
#labeledHeatmap( Matrix = AllnetOverlapColorMat.A, colorLabels = TRUE,
#xLabels = paste0("ME", cMods), yLabels = paste0("ME", rMods),
#textMatrix = AllnetOverlapColorMat.A, colors = blueWhiteRed(100)[50:100],
#main = "Cross-tabulation of A-subgenome modules and aggregate modules",
#cex.text = 1, cex.lab = 1, setStdMargins = FALSE)
#dev.off()

#Dom network
WildnetColorDat <- as.data.frame(WildColorDat)
colnames(WildnetColorDat) <- c("Gene","Aggr")
WildpartColorDat <- data.frame(Gene = colnames(multiExprAll[[3]]$data),Part.Mod = Wildnet12$colors)
WildpartColorDat.A <- WildpartColorDat[grep("_A",WildpartColorDat$Gene),]
WildpartColorDat.D <- WildpartColorDat[grep("_D",WildpartColorDat$Gene),]
colnames(WildpartColorDat.A)[2] <- "Part.Mod.A"
colnames(WildpartColorDat.D)[2] <- "Part.Mod.D"
WildpartColorDat.A$Gene <- gsub("_A","",WildpartColorDat.A$Gene)
WildpartColorDat.D$Gene <- gsub("_D","",WildpartColorDat.D$Gene)
WildpartColorDat.AD <- merge(WildpartColorDat.A, WildpartColorDat.D, by = "Gene", all = TRUE)
WildnetColorDat <- merge(WildnetColorDat,WildpartColorDat.AD, all = TRUE)

rown = length(unique(WildnetColorDat$Aggr))
coln = length(unique(c(WildnetColorDat$Part.Mod.A,WildnetColorDat$Part.Mod.D)))
rMods <- sort(as.numeric(unique(WildnetColorDat$Aggr)))
cMods <- sort(as.numeric(unique(c(WildnetColorDat$Part.Mod.A,WildnetColorDat$Part.Mod.D))))
WildnetOverlapColorMat.A <- matrix(nrow = rown, ncol = coln)
for (rmod in 1:rown) {
		for (cmod in 1:coln){
			WildnetOverlapColorMat.A[rmod, cmod] = sum(WildnetColorDat$Aggr[!is.na(WildnetColorDat$Part.Mod.A)] == rMods[rmod]  & WildnetColorDat$Part.Mod.A[!is.na(WildnetColorDat$Part.Mod.A)] == cMods[cmod] )
		}
}

WildnetOverlapColorMat.D <- matrix(nrow = rown, ncol = coln)
for (rmod in 1:rown) {
		for (cmod in 1:coln){
			WildnetOverlapColorMat.D[rmod, cmod] = sum(WildnetColorDat$Aggr[!is.na(WildnetColorDat$Part.Mod.D)] == rMods[rmod]  & WildnetColorDat$Part.Mod.D[!is.na(WildnetColorDat$Part.Mod.D)] == cMods[cmod] )
		}
}

write.table(WildnetOverlapColorMat.A,file = "WildAggrVsPartAModulesCrossTable.txt",sep = "\t", quote = FALSE, row.names = FALSE)
write.table(WildnetOverlapColorMat.D,file = "WildAggrVsPartDModulesCrossTable.txt",sep = "\t", quote = FALSE, row.names = FALSE)

