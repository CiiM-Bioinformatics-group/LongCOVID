library(Seurat)

myeloidCells[["RNA"]] <- split(myeloidCells[["RNA"]], f = myeloidCells$SampleID)
options(future.globals.maxSize = 30 * 1024 ^ 3) # for 30 Gb RAM
myeloidCells <- NormalizeData(myeloidCells)
myeloidCells <- FindVariableFeatures(myeloidCells)
myeloidCells <- ScaleData(myeloidCells)
myeloidCells <- RunPCA(myeloidCells)
myeloidCells <- IntegrateLayers(object = myeloidCells, method = RPCAIntegration,orig.reduction = "pca", new.reduction = "integrated.rpca",verbose = FALSE, k.weight = 60)

# re-join layers after integration
myeloidCells[["RNA"]] <- JoinLayers(myeloidCells[["RNA"]])

myeloidCells <- FindNeighbors(myeloidCells, reduction = "integrated.rpca", dims = 1:10)
myeloidCells <- FindClusters(myeloidCells, resolution = 0.2)
myeloidCells <- RunUMAP(myeloidCells, dims = 1:10, reduction = "integrated.rpca")
myeloidCellsMarkers = FindAllMarkers(myeloidCells,only.pos = T)
myeloidCellsMarkers_signi = subset(myeloidCellsMarkers, myeloidCellsMarkers$p_val_adj < 0.05)

########subsetting only cd14 cells ######
cd14Mono = subset(myeloidCells, subset = RNA_snn_res.0.2 %in% c(0,1,3))
cd14Mono[["RNA"]] <- split(cd14Mono[["RNA"]], f = cd14Mono$SampleID)
options(future.globals.maxSize = 30 * 1024 ^ 3) # for 30 Gb RAM
cd14Mono <- NormalizeData(cd14Mono)
cd14Mono <- FindVariableFeatures(cd14Mono)
cd14Mono <- ScaleData(cd14Mono)
cd14Mono <- RunPCA(cd14Mono)
cd14Mono <- IntegrateLayers(object = cd14Mono, method = RPCAIntegration,orig.reduction = "pca", new.reduction = "integrated.rpca",verbose = FALSE, k.weight = 60)

# re-join layers after integration
cd14Mono[["RNA"]] <- JoinLayers(cd14Mono[["RNA"]])

cd14Mono <- FindNeighbors(cd14Mono, reduction = "integrated.rpca", dims = 1:10)
cd14Mono <- FindClusters(cd14Mono, resolution = 0.5)
cd14Mono <- RunUMAP(cd14Mono, dims = 1:10, reduction = "integrated.rpca")
cd14MonoMarkers = FindAllMarkers(cd14Mono,only.pos = T)
cd14MonoMarkers_signi = subset(cd14MonoMarkers, cd14MonoMarkers$p_val_adj < 0.05)

#####balData analysis ###########
balMyeloid[["RNA"]] <- split(balMyeloid[["RNA"]], f = balMyeloid$SampleID)
options(future.globals.maxSize = 30 * 1024 ^ 3) # for 30 Gb RAM
balMyeloid <- NormalizeData(balMyeloid)
balMyeloid <- FindVariableFeatures(balMyeloid)
balMyeloid <- ScaleData(balMyeloid)
balMyeloid <- RunPCA(balMyeloid)
balMyeloid <- IntegrateLayers(object = balMyeloid, method = RPCAIntegration,orig.reduction = "pca", new.reduction = "integrated.rpca",verbose = FALSE, k.weight = 60)

# re-join layers after integration
balMyeloid[["RNA"]] <- JoinLayers(balMyeloid[["RNA"]])

balMyeloid <- FindNeighbors(balMyeloid, reduction = "integrated.rpca", dims = 1:10)
balMyeloid <- FindClusters(balMyeloid, resolution = 0.3)
balMyeloid <- RunUMAP(balMyeloid, dims = 1:10, reduction = "integrated.rpca")


balMyeloidMarkers = FindAllMarkers(balMyeloid,only.pos = T)
balMyeloidMarkers_signi = subset(balMyeloidMarkers, balMyeloidMarkers$p_val_adj < 0.05)


####merging data from both tissues #####
bothTissues_myeloidOnly = merge(cd14Mono,y = balMyeloid, add.cell.ids = c("pbmc","bal"),project = "bothTissue",merge.data = TRUE)
bothTissues_myeloidOnly[["RNA"]] = JoinLayers(bothTissues_myeloidOnly[["RNA"]])

metdata = bothTissues_myeloidOnly@meta.data
export = CreateSeuratObject(counts = LayerData(bothTissues_myeloidOnly,assay = "RNA",layer = "counts"),project = "Export",meta.data = metdata)
export[["RNA"]] <- split(export[["RNA"]], f = export$SampleID)
options(future.globals.maxSize = 30 * 1024 ^ 3) # for 30 Gb RAM
export <- NormalizeData(export)
export <- FindVariableFeatures(export)
export <- ScaleData(export)
export <- RunPCA(export)
export <- IntegrateLayers(object = export, method = RPCAIntegration,orig.reduction = "pca", new.reduction = "integrated.rpca",verbose = FALSE, k.weight = 60)

# re-join layers after integration
export[["RNA"]] <- JoinLayers(export[["RNA"]])

export <- FindNeighbors(export, reduction = "integrated.rpca", dims = 1:15)
export <- FindClusters(export, resolution = c(0.2,0.3,0.4,0.5))
export <- RunUMAP(export, dims = 1:15, reduction = "integrated.rpca")

table(export$RNA_snn_res.0.3)

pro-fibroticMarkers = c("CD14","FCGR3A","S100A8","S100A12","FCN1","VCAN","FPR1","CD93","SELL","CCR2","DDX58","NFKBIA","SERPINB9","CXCL8",
                       "NAMPT","SLC25A37","APOBEC3A","TGFB1","ITGAM","CLEC5A","MYC","DDX21","DPYD","TYMP","GNG5","CCT5","LSM7","CALM1",
                       "CCL2","TGFBI","OLFML2B","CCL13","CD163","CD84","MERTK","LGMN","SDC3","SPP1","PLA2G7","NRP1","NPL","MARCKS","MAF","A2M",
                       "ADAP2","IDH1","FMN1","CALR","PMP22","HSPA5","SLAMF8","SPRED1","HS3ST1","GPNMB","LIPA","MATK","HTRA4","COL8A2","ABCC3","MRC1",
                       "TREM2","APOE","ACP5","GCHFR","CD9","EVL","GSN","CD52","APOC1","MARCO","FBP1","HLA-DRA","ALOX5AP","LPL","ITGB8","CES1","INHBA",
                       "GPD1","MME","VMO1","NUPR1","FN1","PPIC","S100A13","PLA2G16","FABP4","MKI67","TOP2A","CENPA","BIRC5")

myeloidcounts = GetAssayData(cd14Mono,slot = "counts")
myeloid_cell_rankings <- AUCell_buildRankings(myeloidcounts)
load("MonocyteCS_dgeList.RData")

cs4_dge = Monocytes_CS_dges_subset$CS4Markers

cs4_dge = subset(cs4_dge, cs4_dge$avg_log2FC > 0.7)

myeloidCells_cs4AUC = AUCell_calcAUC(cs4_dge$rowname,myeloid_cell_rankings)

myeloidCells_cs4AUCscore = t(getAUC(myeloidCells_cs4AUC))
cd14Mono$CS4AUC =myeloidCells_cs4AUCscore[,1]

exportcounts = GetAssayData(export,slot = "counts")
export_cell_rankings <- AUCell_buildRankings(exportcounts)

export_cs4AUC = AUCell_calcAUC(cs4_dge$rowname,export_cell_rankings)

export_cs4AUCscore = t(getAUC(export_cs4AUC))
export$CS4AUC =export_cs4AUCscore[,1]


