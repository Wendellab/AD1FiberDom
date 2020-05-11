#Correspondence table code
wildPartLabels = Wildnet12$colors
domPartLabels = Domnet12$colors

wildPartGeneTree = Wildnet12$dendrograms[[1]]
domPartGeneTree = Domnet12$dendrograms[[1]]

wildPartColors = labels2colors(wildPartLabels)
domPartColors = labels2colors(domPartLabels)

wildPartMEs = Wildnet12$MEs
domPartMEs = Domnet12$MEs

wildPartOrderedMEs = orderMEs(wildPartMEs, greyName= "ME0")
domPartOrderedMEs = orderMEs(domPartMEs, greyName= "ME0")

wildPartModuleLabels = substring(names(wildPartOrderedMEs), 3)
domPartModuleLabels = substring(names(domPartOrderedMEs), 3)

wildPartModules = labels2colors(as.numeric(wildPartModuleLabels))
domPartModules = labels2colors(as.numeric(domPartModuleLabels))

nWildPartMods = length(wildPartModules)
nDomPartMods = length(domPartModules)

pTable = matrix(0, nrow = nWildPartMods, ncol = nDomPartMods);
CountTbl = matrix(0, nrow = nWildPartMods, ncol = nDomPartMods);

for (wmod in 1:nWildPartMods)
for (dmod in 1:nDomPartMods)
{
wildMembers = (wildPartColors == wildPartModules[wmod]);
domMembers = (domPartColors == domPartModules[dmod]);
pTable[wmod, dmod] = -log10(fisher.test(wildMembers, domMembers, alternative = "greater")$p.value);
CountTbl[wmod, dmod] = sum(wildPartColors == wildPartModules[wmod] & domPartColors ==
domPartModules[dmod])
}

pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
pTable[pTable>50 ] = 50 ;

wildModTotals = apply(CountTbl, 1, sum)
domModTotals = apply(CountTbl, 2, sum)

sizeGrWindow(10,7 );
par(mfrow=c(1,1));
par(cex = 1.0);
par(mar=c(5, 30, 4, 3) + 0.3);

pdf(file = "WildVsDomPartModules.pdf", wi =8.5, he = 11);
labeledHeatmap(Matrix = pTable,
xLabels = paste(" ", domPartModules),
yLabels = paste(" ", wildPartModules),
colorLabels = TRUE,
xSymbols = paste("Dom ", domPartModuleLabels, ": ", domModTotals, sep=""),
ySymbols = paste("Wild ", wildPartModuleLabels, ": ", wildModTotals, sep=""),
textMatrix = CountTbl,
colors = blueWhiteRed(100)[50:100],
main = "Correspondence of wild and domesticated modules in partitioned data set",
cex.text = 0.3, cex.lab = 0.5, setStdMargins = FALSE);
dev.off();
