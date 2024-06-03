library(WGCNA)
datExpr<-t(datExpr)

library(readr)
dis2<-read_csv("dis2.csv") #reorganizing the dataframe
dis2<-dis2[ !dis2$ID %in% remove, ]
asp<-subset(datExpr,rownames(datExpr) %in% dis2$ID[(1:23)])
Autism<-subset(datExpr,rownames(datExpr) %in% dis2$ID[(24:91)])
#Checking for any samples that might be outliers.
sampleTree = hclust(dist(Autism), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
#removing samples outliers
outlier<-c("GSM650514")
Autism<-Autism[!(rownames(Autism) %in% outlier),]
#Rechecking
sampleTree = hclust(dist(Autism), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

#start builiding WGCNA network
powers = c(c(1:30), seq(from = 1, to=30, by=5))
sft1 = pickSoftThreshold(asp, powerVector = powers, verbose = 1, networkType = "signed") #soft-thresholding power for network construction is 15 
sft2 = pickSoftThreshold(Autism, powerVector = powers, verbose = 1, networkType = "signed") #soft-thresholding power for network construction is 18 
net<- blockwiseModules(asp, power = 15,networkType="signed",TOMType = "signed",minModuleSize = 50,deepSplit = 2,verbose = 3)
net2<- blockwiseModules(Autism, power = 18,networkType="signed",TOMType = "signed",minModuleSize = 50,deepSplit = 2,verbose = 3)
#extracting modules
colasp<-labels2colors(net$colors)
colautism<-labels2colors(net2$colors)

#data preparation
setLabels = c("Autism","Asperger")
multiExpr = list(Autism = list(data = Autism),Asperger = list(data = asp))
multiColor=list(Autism=colautism,Asperger=colasp)

#preservation analysis where Autism is the reference network. The analysis check if the network in the reference is preserved into the test network
system.time( {
  mpAut = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = 1,
                          nPermutations = 1000,
                          randomSeed = 1,
                          quickCor = 0,
                          verbose = 3,networkType = "signed",maxModuleSize = 4000)
} );


#ASperger's syndrome is the reference
system.time( {
  mpAsp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = 2,
                          nPermutations = 1000,
                          randomSeed = 1,
                          quickCor = 0,
                          verbose = 3,networkType = "signed",maxModuleSize = 4000)
} );

#Dendrogram Figure and overlap table

nSets = 2
ref=1
test=2
library(flashClust)
dendrograms = list();
for (set in 1:nSets)
{
  adj = abs(cor(multiExpr[[set]]$data, use = "p"))^9;
  dtom = TOMdist(adj);
  dendrograms[[set]] = flashClust(as.dist(dtom), method = "a");
}
# Get eigengenes
mes = list()
for (set in 1:nSets)
{
  mes[[set]] = moduleEigengenes(multiExpr[[set]]$data, multiColor[[ref]])$eigengenes

overlap = overlapTable(colasp, colautism)
numMat = -log10(overlap$pTable);
numMat[numMat >50] = 50;
# Prepare for generating a color-coded plot of the overlap table. The text of the table will consist of
# counts and corresponding p-values.
textMat = paste(overlap$countTable, "\n", signif(overlap$pTable, 2));
dim(textMat) = dim(numMat)
# Additional information for the plot. These will be used shortly.
xLabels = paste("M", sort(unique(colautism)));
yLabels = paste("M", sort(unique(colasp)));
ySymbols = paste(sort(unique(colasp)), ": ", table(colasp), sep = "")
xSymbols = paste(sort(unique(colautism)), ": ", table(colautism), sep = "")



# Open a graphical window. If plotting into a file, skip the next line.
sizeGrWindow(7, 7); fp = FALSE;
#pdf(fi = spaste("Plots/Asperger's syndromeAutism-motivationFigure-dendrosAndTable.pdf"), w = 7, h = 7.0); fp = TRUE
layout(matrix(c(1,2,5, 3,4,5), 3, 2),
       heights = c(3, 1, 5.5), widths = c(1, 1));
#layout.show(5);
par(mgp = c(3, 2, 0));
plotDendroAndColors(dendrograms[[1]],
                    cbind(colasp, colautism),
                    c("Asperger's syndrome modules", "Autism modules"),
                    setLayout = FALSE,
                    marAll = c(1, 6, 2.7, 0.2),
                    addGuide = FALSE,
                    main = "A. Asperger's syndrome gene dendrogram\nand module colors", cex.main = 1.2,
                    dendroLabels = FALSE, hang = 0.03, cex.colorLabels = 0.7, abHeight = 0.95);
par(mgp = c(2, 1, 0));
plotDendroAndColors(dendrograms[[2]],
                    cbind(colasp, colautism),
                    c("Asperger's syndrome modules", "Autism modules"),
                    setLayout = FALSE,
                    marAll = c(1, 6, 2.7, 0.2),
                    addGuide = FALSE,
                    main = "B. Autism gene dendrogram\nand module colors", cex.main = 1.2,
                    dendroLabels = FALSE, hang = 0.03, cex.colorLabels = 0.7, abHeight = 0.95);
# Plot the overlap table
fcex = 0.3;
pcex = 0.3;
fcexl = 0.3;
pcexl = 0.3;
par(mar = c(6, 7, 2, 1.0));
labeledHeatmap(Matrix = numMat,
               xLabels = xLabels, xSymbols = xSymbols,
               yLabels = yLabels, ySymbols = ySymbols,
               colorLabels = TRUE,
               colors = greenWhiteRed(100)[50:100],
               textMatrix = textMat, cex.text = if (fp) fcex else pcex, setStdMargins = FALSE,
               cex.lab = if (fp) fcexl else pcexl,
               xColorWidth = 0.08,
               main = "C. Asperger's syndrome modules (rows) vs. Autism modules (columns)", cex.main = 1.2)
}

#extracting Z summary statistics
write.csv(mpAut[["preservation"]][["Z"]][["ref.Autism"]][["inColumnsAlsoPresentIn.Asperger"]], "autz.csv")
write.csv(mpAut[["quality"]][["observed"]][["ref.Autism"]][["inColumnsAlsoPresentIn.Asperger"]],"autmedianrank.csv")

write.csv(mpAsp[["quality"]][["observed"]][["ref.Asperger"]][["inColumnsAlsoPresentIn.Autism"]], "aspmedian.csv")
write.csv(mpAsp[["preservation"]][["Z"]][["ref.Autism"]][["inColumnsAlsoPresentIn.Asperger"]],"aspz.csv")

#enrichment
  enableWGCNAThreads(n=3)
  TOM = TOMsimilarityFromExpr(asp, power = 15);
  modules = c("lightyellow");
  # Select module probes
  probes = colnames(datExpr)
  moduleColors=colasp
  inModule = is.finite(match(colasp, modules));
  modProbes = probes[inModule];
  modGenes = probes[inModule]
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)
  # Export the network into edge and node list files Cytoscape can read
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("lightyellowCytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("lightyellowCytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.02,
                                 nodeNames = modProbes,
                                 altNodeNames = modGenes,
                                 nodeAttr = moduleColors[inModule]);

  #Enrichment analysis
  library(magrittr)
  library(cowplot)
  library(clusterProfiler)
  library(DOSE)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(brainImageR)
  library(ReactomePA)
  library(enrichplot)
  library(ggplot2)
  library(enrichR)
  modules = c("lightyellow");
  # Select module probes
  probes = colnames(datExpr)
  inModule = is.finite(match(colasp, modules));
  modProbes = probes[inModule];
  modGenes = probes[inModule]
  ly<-as.data.frame(modGenes)
  ly$gene<-mapIds(org.Hs.eg.db,
              keys = modGenes,
              column = "ENTREZID",
              keytype = "SYMBOL",
              multiVals = "first")
  
  #identify top 80 hub genes
  nTop = 80;
  IMConn = softConnectivity(asp[, modProbes]);
  top = (rank(-IMConn) <= nTop)
  g<-as.data.frame(top)
  inputt<-ly[!(g$top == FALSE),] #exported the genelist to gene mania for gene network visualization
#GO enrichment CC
  lyCC <- enrichGO(gene          = ly$gene,OrgDb         = org.Hs.eg.db,
                  keyType       = 'ENTREZID',
                  ont           = "CC",
                 
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)
  
  lyBB <- enrichGO(gene          = ly$gene,OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                
                 pvalueCutoff  = 0.015,
                 qvalueCutoff  = 0.05
  )
  
  lyMF <- enrichGO(gene          = ly$gene,OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "MF",
                 
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.015,
                 qvalueCutoff  = 0.05)
  
  a<-dotplot(ego,showCategory=15)
  b<-dotplot(bp,showCategory=15)
  c<-dotplot(MF,showCategory=15)
 
  d <- enrichPathway(gene=g$y, pvalueCutoff = 0.05, readable=TRUE)
  dd<-dotplot(d,showCategory=10)
  theme_set(theme_cowplot(font_size=10))
  plot_grid(a, b, c,dd,
            labels = c("A", "B", "C","D"),
            ncol = 2, nrow = 2)
  #brain enrichment analysis
  brainImageR:::loadworkspace()
  composite <- SpatialEnrichment(genes = ly$modGenes, reps = 2, refset = "developing")
  res <- testEnrich(composite, method = "fisher")
  composite <- CreateBrain(composite, res, slice = 5, pcut =1)
PlotBrain(composite)
  #disease enrichment analysis
  websiteLive <- TRUE
  dbs <- c("DisGeNET")
  if (websiteLive) {
    enriched <- enrichr(ly$modGenes, dbs)
  }
  enriched$DisGeNET
  if (websiteLive) plotEnrich(enriched[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
  
  
  #salmon module enrichment analysis
  
  modules = c("salmon");
  # Select module probes
  inModule = is.finite(match(colasp, modules));
  modProbes = probes[inModule];
  modGenes = probes[inModule]
  salmon<-as.data.frame(modGenes)
  salmon$entrezID<-mapIds(org.Hs.eg.db,
              keys = modGenes,
              column = "ENTREZID",
              keytype = "SYMBOL",
              multiVals = "first")
  
  #identify top 80 hub genes
  nTop = 80;
  IMConn = softConnectivity(asp[, modProbes]);
  top = (rank(-IMConn) <= nTop)
  g<-as.data.frame(top)
  inputt<-salmon[!(g$top == FALSE),] #exported the genelist to gene mania for gene network visualization
  #GO enrichment CC
salmonCC <- enrichGO(gene          = salmon$entrezID,OrgDb         = org.Hs.eg.db,
                  keyType       = 'ENTREZID',
                  ont           = "CC",
                  
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)
  
 salmonBP <- enrichGO(gene          = salmon$entrezID,OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05
  )
  
  salmonMF <- enrichGO(gene          = salmon$entrezID,OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "MF",
                 
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
  
  a<-dotplot(ego,showCategory=15)
  b<-dotplot(bp,showCategory=15)
  c<-dotplot(MF,showCategory=15)
  d <- enrichPathway(gene=salmon$entrezID, pvalueCutoff = 0.05, readable=TRUE)
  dd<-dotplot(d,showCategory=10)
  theme_set(theme_cowplot(font_size=10))
  plot_grid(a, b, c,dd,
            labels = c("A", "B", "C","D"),
            ncol = 2, nrow = 2)
  
#Cytoscape network files for network visualization
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)
  # Export the network into edge and node list files Cytoscape can read
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("salmonCytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("salmonCytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.02,
                                 nodeNames = modProbes,
                                 altNodeNames = modGenes,
                                 nodeAttr = moduleColors[inModule]);
  
  
  #brain enrichment analysis
  composite <- SpatialEnrichment(genes = salmon$modGenes, reps = 2, refset = "developing")
  res <- testEnrich(composite, method = "fisher")
  composite <- CreateBrain(composite, res, slice = 5, pcut =1)
  salmonbrai<-PlotBrain(composite)
  #disease enrichment analysis
  websiteLive <- TRUE
  dbs <- c("DisGeNET")
  if (websiteLive) {
    enriched <- enrichr(salmon$modGenes, dbs)
  }
  enriched$DisGeNET
  if (websiteLive) plotEnrich(enriched[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
  

#green yellow
  modules = c("darkturquoise");
  # Select module probes
  inModule = is.finite(match(colasp, modules));
  modProbes = probes[inModule];
  modGenes = probes[inModule]
  DT<-as.data.frame(modGenes)
  DT$entrezID<-mapIds(org.Hs.eg.db,
              keys = modGenes,
              column = "ENTREZID",
              keytype = "SYMBOL",
              multiVals = "first")
 
  #identify top 80 hub genes  
   nTop = 80;
  IMConn = softConnectivity(asp[, modProbes]);
  top = (rank(-IMConn) <= nTop)
  g<-as.data.frame(top)
  inputt<-DT[!(g$top == FALSE),] #exported the genelist to gene mania for gene network visualization
  #GO enrichment CC
  DTCC <- enrichGO(gene          = DT$entrezID,OrgDb         = org.Hs.eg.db,
                  keyType       = 'ENTREZID',
                  ont           = "CC",
                  
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)
  
  DTbp <- enrichGO(gene          = DT$entrezID,OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05
  )
  
  DTMF <- enrichGO(gene          = DT$entrezID,OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "MF",
                 
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
  
  a<-dotplot(ego,showCategory=12)
  b<-dotplot(bp,showCategory=15)
  c<-dotplot(MF,showCategory=15)
  d <- enrichPathway(gene=DT$entrezID, pvalueCutoff = 0.05, readable=TRUE)
  dd<-dotplot(d,showCategory=10)
  theme_set(theme_cowplot(font_size=10))
  plot_grid(a, b, c,dd,
            labels = c("A", "B", "C","D"),
            ncol = 2, nrow = 2)
  
  #Cytoscape network files for network visualization
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)
  # Export the network into edge and node list files Cytoscape can read
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("darkturCytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("darkturCytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.02,
                                 nodeNames = modProbes,
                                 altNodeNames = modGenes,
                                 nodeAttr = moduleColors[inModule]);
  
  
  #brain enrichment analysis
  composite <- SpatialEnrichment(genes = DT$entrezID, reps = 2, refset = "developing")
  res <- testEnrich(composite, method = "fisher")
  composite <- CreateBrain(composite, res, slice = 5, pcut =1)
PlotBrain(composite)
  #disease enrichment analysis
  websiteLive <- TRUE
  dbs <- c("DisGeNET")
  if (websiteLive) {
    enriched <- enrichr(DT$entrezID, dbs)
  }
  enriched$DisGeNET

  
