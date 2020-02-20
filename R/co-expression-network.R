require(matrixStats)
require(data.table)
require(igraph)
require(Hmisc)
require(compiler)
require(reshape2 )
### main functions ###
source("/Users/ida.larsson/Desktop/Exjobb/Uppmax/bin/analysis.function.coexp.net.r")
filterExprMat <- function(exprMat, meanCut, mode="mean") {
  if (mode == "mean") rMeans = rowMeans(exprMat)
  if (mode == "median") rMeans = rowMedians(exprMat)
  return(exprMat[rMeans>meanCut,])
}
makeCorTable <- function(exprMat, cutOff=0.99, mode="pearson", self=F, debug=F) {
  minValue = -2
  if (mode == "pearson") {
    corMat = cor(t(exprMat))
  }
  if (mode == "incomplete") {
    corMat = cor(t(exprMat), use="pairwise.complete.obs")
  }
  if (mode == "spearman") {
    corMat = cor(t(exprMat), method = "spearman")
  }
  if (mode == "r2") {
    corMat = cor(t(exprMat))
    corMat = corMat^2
  }
  ### checking NA
  if (debug) print("na amount:")
  if (debug) print(table(is.na(corMat)))
  corMat[is.na(corMat)]=minValue
  cutOff = as.numeric( quantile(corMat[lower.tri(corMat)], cutOff, na.rm=T) )
  if(debug) print("quantile: ")
  if(debug) print(cutOff)
  resTable = melt(corMat)[,1:3]
  resTableNew = as.data.table(resTable)
  setkey(resTableNew, value)
  resTableNew = resTableNew[order(value, decreasing = T)]
  resTable = as.data.frame(resTableNew)
  if (!self) {
    resTable = resTable[resTable[,1]!=resTable[,2],]
  }
  if(debug) print("cut-off before")
  if(debug) print(dim(resTable))
  resTable = resTable[resTable[,3]>=cutOff,]
  if(debug) print("cut-off on-the-way")
  if(debug) print(dim(resTable))
  resTable = resTable[!is.na(resTable[,3]),]
  resTable = resTable[resTable[,3]!=minValue,]
  if(debug) print("cut-off after")
  if(debug) print(dim(resTable))
  return(resTable)
}

makeCorNet <- function(corTable) {
  corNet = igraph::graph.data.frame(corTable[,1:2], directed=F)
  corNet = igraph::simplify(corNet, remove.multiple=T, remove.loops=T)
  return(corNet)
}

makeModuleList <- function(corNet, debug=F) {
  fc = cluster_walktrap(corNet)
  cluster = fc$membership
  geneCluster = data.frame(gene=V(corNet)$name, cluster=cluster, stringsAsFactors=F)
  geneClusterList = split.data.frame(geneCluster, geneCluster$cluster)
  geneClusterList = lapply(geneClusterList, "[[", "gene")
  geneClusterSizes = do.call(c, lapply(geneClusterList, length))
  if (debug) {
    print("... done")
    print(paste("... modularity:", as.character(modularity(fc))))
    print(paste("... no. clusters:", as.character(length(geneClusterList))))
    print(paste("... no. genes of max cluster:", as.character(sort(geneClusterSizes,T)[1])))
  }
  return(geneClusterList)
}

annotateModulesByCC <-function(corNet, moduleList, cutCluster=5, cutCC = 0.5, debug=F) {
  modulePrunedList = getModulePruned(moduleList, cutCluster)
  moduleNames = names(modulePrunedList)
  # module matrix
  moduleMat = getModuleMatrix(corNet, modulePrunedList)
  print(moduleMat)
  # first axis
  moduleCC = getClusterScores(corNet, modulePrunedList)
  highCCModuleNames = names(moduleCC[moduleCC >= cutCC])
  # summary of two axes
  summaryTypes = rep("None", length(modulePrunedList))
  names(summaryTypes) = moduleNames
  summaryTypes[moduleNames %in% highCCModuleNames ] = "HighCC"
  summaryTypes[!moduleNames %in% highCCModuleNames ] = "LowCC"
  moduleCut=1
  moduleMelt = melt(moduleMat)
  moduleMelt = moduleMelt[moduleMelt[,3]>moduleCut,]
  moduleMelt[] = lapply(moduleMelt, as.character)
  modulesShown = unique.default(c(moduleMelt[,1], moduleMelt[,2]))
  moduleSizes = unlist(lapply(modulePrunedList, length))
  moduleAttr = data.frame(node=moduleNames, size=moduleSizes, summary = summaryTypes, cc = moduleCC, stringsAsFactors = F)
  modulesNotShown = moduleNames[!moduleNames %in% modulesShown]
  return(list(nodeTable = moduleAttr, edgeTable=moduleMelt, nodesNotShown = modulesNotShown ))
}

annotateModulesByThreeAxes <-function(corNet, moduleList, clockGenes, tissueSet, cutCluster=5, cutCC = 0.5, gsPval = 0.01, checkMode = F, debug=F) {
  modulePrunedList = getModulePruned(moduleList, cutCluster)
  moduleNames = names(modulePrunedList)
  # module matrix
  moduleMat = getModuleMatrix(corNet, modulePrunedList)
  # first axis
  moduleCC = getClusterScores(corNet, modulePrunedList)
  highCCModuleNames = names(moduleCC[moduleCC >= cutCC])
  # second axis
  if (checkMode) {
    clockIncluded = unlist(lapply(modulePrunedList, function(x) any(x %in% clockGenes)))
    clockModuleNames = moduleNames[clockIncluded]
  }
  else {
    clockSig = getEnrichMentGsByList(modulePrunedList, clockGenes)
    clockModuleNames = moduleNames[clockSig <= gsPval]
  }
  #clockIncluded = unlist(lapply(modulePrunedList, function(x) any(x %in% clockGenes)))
  #clockModuleNames = moduleNames[clockIncluded]
  enrichMat = getEnrichMentGsByGs(modulePrunedList, tissueSet)
  enrichMatSig = (enrichMat <= gsPval)*1
  # third axis
  tissueInfo = sapply(1:length(modulePrunedList), function(ind) {
  currGs = enrichMatSig[,ind]
  names(currGs) = names(tissueSet)
  sigGs = names(currGs[currGs==1])
   if (length(sigGs)==0) return("None")
    return(getStrFromGenes_cmp(unlist(sigGs)))
  })
  # summary of two axes
  summaryTypes = rep("None", length(modulePrunedList))
  names(summaryTypes) = moduleNames
  summaryTypes[moduleNames %in% highCCModuleNames & moduleNames %in% clockModuleNames] = "HighCC.Circadian"
  summaryTypes[moduleNames %in% highCCModuleNames & !moduleNames %in% clockModuleNames] = "HighCC"
  summaryTypes[!moduleNames %in% highCCModuleNames & moduleNames %in% clockModuleNames] = "Circadian"
  moduleCut=1
  moduleMelt = melt(moduleMat)
  moduleMelt = moduleMelt[moduleMelt[,3]>moduleCut,]
  moduleMelt[] = lapply(moduleMelt, as.character)
  modulesShown = unique.default(c(moduleMelt[,1], moduleMelt[,2]))
  moduleSizes = unlist(lapply(modulePrunedList, length))
  moduleAttr = data.frame(node=moduleNames, size=moduleSizes, summary = summaryTypes, tissueInfo = tissueInfo, cc = moduleCC, stringsAsFactors = F)
  modulesNotShown = moduleNames[!moduleNames %in% modulesShown]
  return(list(nodeTable = moduleAttr, edgeTable=moduleMelt, nodesNotShown = modulesNotShown ))
}

write.edge.cytoscape <- function(edgeTable, nodesNotShown, outFile) {
  lenTable=dim(edgeTable)[1]
  lenNodesNotShown = length
  print(head(edgeTable))
  sink(outFile, append=F)
  cat("source")
  cat("\t")
  cat("target")
  cat("\n")
  for (i in 1:lenTable) {
    # ints = unlist(edgeTable[i,])
    cat(edgeTable[i,1])
    cat("\t")
    cat(edgeTable[i,2])
    cat("\n")
  }
  for (node in nodesNotShown) {
    cat(node)
    cat("\n")
  }
  sink()
}

getGeneSetAsList <- function(geneSetFile ) {
  geneSetData = read.table(geneSetFile, header=F, stringsAsFactors = F, sep="\t")
  geneSetNames = geneSetData[,1]
  #if (change) geneSetNames = fixGeneSetNames(geneSetNames)
  geneSetGeneStrs = geneSetData[,2]
  geneSetGenes = sapply(geneSetGeneStrs, getGenesFromStr_cmp)
  names(geneSetGenes) = geneSetNames
  return(geneSetGenes)
}

getInteractionsFromModules <- function( coexpTable, moduleList, targets, debug=T ) {
  moduleNames = names(moduleList)
  if (!all(targets %in% moduleNames)) {
     print("... some targets were not shown in modulePruned")
    print("... please check it again")
    return(NULL)
  }
  targetModules = moduleList[targets]
  targetModuleGenes = unique.default(do.call(c, targetModules))
  targetModuleMaps = rep(NA, length(targetModuleGenes))
  for (i in 1:length(targets)) {
    #currTarget = targets[1]
    currTarget = targets[i]
    currTargetGenes = unlist(moduleList[[currTarget]])
    targetModuleMaps[targetModuleGenes %in% currTargetGenes] = currTarget
  }
  targetInd = coexpTable[,1]%in%targetModuleGenes & coexpTable[,2]%in%targetModuleGenes
  targetEdgeTable = coexpTable[targetInd, 1:3]
  colnames(targetEdgeTable) = c("Source", "Target","Corr")
  targetNodeTable = data.frame(index=targetInd,
  node=targetModuleGenes,
  symbol=getSymbols(targetModuleGenes),
  moduleType=targetModuleMaps,
  stringsAsFactors=F)
  return(list(nodeTable = targetNodeTable, edgeTable = targetEdgeTable))
}

### execution ###
#data pre-processing
gbmGenes = read.csv("/Users/ida.larsson/Documents/Exjobb/R/SurvivalInputForGBM.txt", stringsAsFactors = F,check.names = F, sep='\t')
gbmGenes = gbmGenes[7:19577,]
rownames(gbmGenes) = gbmGenes[,1]
gbmGenes = gbmGenes[,-1]
for (i in 1:153){
  gbmGenes[,i] = as.numeric(gbmGenes[,i])
}

#filter out lowly expressed genes
exprMat = filterExprMat(gbmGenes,meanCut = 1,mode = "mean")

#correlation test
corTable = makeCorTable(exprMat,cutOff = 0.99,mode = "pearson",self = F, debug = F)

#igraph
corNet = makeCorNet(corTable)

#cluster using random walk method
moduleList = makeModuleList(corNet, debug = F)

#prepape for cytoscape
cytoscapeMaterial = annotateModulesByCC(corNet,moduleList,cutCluster=5,cutCC=0.5,debug=F)

#cytoscape
write.edge.cytoscape(cytoscapeMaterial$edgeTable,cytoscapeMaterial$nodesNotShown, 'cytoEdgeFile.txt')
write.node.cytoscape(cytoscapeMaterial$nodeTable, 'cytoNodeFile.txt')
