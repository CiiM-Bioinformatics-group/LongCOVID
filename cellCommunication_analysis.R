library(CellChat)
library(Seurat)
library(Signac)
library(future)

celltypes = c("cd8T.C0", "cd8T.C1","cd8T.C2","cd8T.C3","cd8T.C4","MC1","MC2","MC3","MC4","NK.C0","NK.C1","NK.C2","NK.C3","NK.C4")

timePoints = c("R_Recov","nonICU_T2","nonICU_T3","nonICU_T4")

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search =c("Secreted Signaling","Cell-Cell Contact"), key = "annotation")
showDatabaseCategory(CellChatDB.use)
data.input <- GetAssayData(bigObj, assay = "RNA", layer = "data") # normalized data matrix
metaAll <- bigObj@meta.data

bigCellChatList_perTP = list()
for (j in timePoints) {
  print(j)
  cell.use = rownames(metaAll)[metaAll$AcuteCovStatus2 == j & metaAll$CellCluster %in% celltypes]
  dataSubset = data.input[, cell.use]
  meta = metaAll[cell.use, ]
  meta$samples = meta$poolID
  cellchat <- createCellChat(object = dataSubset, meta = meta, group.by = "CellCluster")
  cellchat <- addMeta(cellchat, meta = meta)
  cellchat <- setIdent(cellchat, ident.use = "CellCluster") # set "labels" as default cell identity
  levels(cellchat@idents) # show factor levels of the cell labels
  groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database # nolint: line_length_linter.
  future::plan("multicore", workers = 4) 
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat,type =  "truncatedMean", trim = 0.15)
  cellchat <- filterCommunication(cellchat, min.cells = 30)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  bigCellChatList_perTP[[j]] = cellchat
  #networkDFList[[j]] = df.net
}

nonICUOnlyList = subset(bigCellChatList_perTP,names(bigCellChatList_perTP) %in% c("nonICU_T2","nonICU_T3","nonICU_T4","R_Recov"))


weight.max <- getMaxWeight(nonICUOnlyList, attribute = c("idents","count"))

gg <- list()
num.link <- sapply(nonICUOnlyList, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
for (i in 1:length(nonICUOnlyList)) {
  nonICUOnlyList[[i]] = netAnalysis_computeCentrality(nonICUOnlyList[[i]],slot.name = "netP")
  gg[[i]] <- netAnalysis_signalingRole_scatter(nonICUOnlyList[[i]], title = names(nonICUOnlyList)[i], weight.MinMax = weight.MinMax)
  #pp[[i]] <- netAnalysis_signalingChanges_scatter(nonICUOnlyList[[i]],idents.use = "CD8 T cells", comparison = c(4,1))
}

library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- Reduce(union,list(nonICUOnlyList[[i]]@netP$pathways,
                                   nonICUOnlyList[[i+1]]@netP$pathways,
                                   nonICUOnlyList[[i+2]]@netP$pathways))
ht1 = netAnalysis_signalingRole_heatmap(nonICUOnlyList[[i]], pattern = "outgoing", signaling = pathway.union,color.heatmap = "Reds", title = names(nonICUOnlyList)[i], width = 7, height = 8)
ht2 = netAnalysis_signalingRole_heatmap(nonICUOnlyList[[i+1]], pattern = "outgoing", signaling = pathway.union,color.heatmap = "Reds", title = names(nonICUOnlyList)[i+1], width = 7, height = 8)
ht3 = netAnalysis_signalingRole_heatmap(nonICUOnlyList[[i+2]], pattern = "outgoing", signaling = pathway.union,color.heatmap = "Reds", title = names(nonICUOnlyList)[i+2], width = 7, height = 8)
#ht4 = netAnalysis_signalingRole_heatmap(nonICUOnlyList[[i+3]], pattern = "outgoing", signaling = pathway.union,color.heatmap = "Reds", title = names(nonICUOnlyList)[i+3], width = 7, height = 8)

ht6 = netAnalysis_signalingRole_heatmap(nonICUOnlyList[[i]], pattern = "incoming", signaling = pathway.union, color.heatmap = "Blues",title = names(nonICUOnlyList)[i], width = 7, height = 8)
ht7 = netAnalysis_signalingRole_heatmap(nonICUOnlyList[[i+1]], pattern = "incoming", signaling = pathway.union,color.heatmap = "Blues", title = names(nonICUOnlyList)[i+1], width = 7, height = 8)
ht8 = netAnalysis_signalingRole_heatmap(nonICUOnlyList[[i+2]], pattern = "incoming", signaling = pathway.union,color.heatmap = "Blues", title = names(nonICUOnlyList)[i+2], width = 7, height = 8)





