library(Seurat)
library(Signac)
bigObj = readRDS("/vol/projects/skumar/LongCovid/bothCohortsMergedAnalysis/bothCohorts_Annotated.rds")

celltypes = c("B cells","NK cells","CD14 Monocytes","CD16 Monocytes", "CD4 T cells", "CD8 T cells")
dge_vsR_PerTP_perCT_perCat = list()
dge_vsT1_PerTP_perCT_perCat = list()
dge_vsHC_PerTP_perCT_perCat = list()
dge_T1vsHC_perCT = list()
timeP = c("T2","T3","T4","T5")

for (i in celltypes) {
  for (j in timeP){
    ident1 = paste0(i,"_nonsevere_",j)
    ident2 = paste0(i,"_T1")
    print(paste0("Comparison ",ident1," vs ",ident2))
    LCvsT1_cats = FindMarkers(bigObj,ident.1 = ident1,ident.2 = ident2,group.by = "CellType_AcuteCov_TP")
    LCvsT1_cats = subset(LCvsT1_cats, LCvsT1_cats$p_val_adj < 0.05)
    if(nrow(LCvsT1_cats) > 0){
      LCvsT1_cats$Comparison = "LCVsAcute_nonsevere"
      dge_vsT1_PerTP_perCT_perCat[[paste0(i,"_",j,"_nonsevere")]] = LCvsT1_cats
    }
    
    ident2 = paste0(i,"_R_Recov")
    print(paste0("Comparison ",ident1," vs ",ident2))
    LCvsR_cats = FindMarkers(bigObj,ident.1 = ident1,ident.2 = ident2,group.by = "CellType_AcuteCov_TP")
    LCvsR_cats = subset(LCvsR_cats, LCvsR_cats$p_val_adj < 0.05)
    if(nrow(LCvsR_cats) > 0){
      LCvsR_cats$Comparison = "LCVsRecov_nonsevere"
      dge_vsR_PerTP_perCT_perCat[[paste0(i,"_",j,"_nonsevere")]] = LCvsR_cats
    }
    
    ident2 = paste0(i,"_HC_HC")
    print(paste0("Comparison ",ident1," vs ",ident2))
    LCvsHC_cats = FindMarkers(bigObj,ident.1 = ident1,ident.2 = ident2,group.by = "CellType_AcuteCov_TP")
    LCvsHC_cats = subset(LCvsHC_cats, LCvsHC_cats$p_val_adj < 0.05)
    if(nrow(LCvsHC_cats) > 0){
      LCvsHC_cats$Comparison = "LCVsHC_nonsevere"
      dge_vsHC_PerTP_perCT_perCat[[paste0(i,"_",j,"_nonsevere")]] = LCvsHC_cats
    }
    
    
    ident1 = paste0(i,"_severe_",j)
    ident2 = paste0(i,"_T1")
    print(paste0("Comparison ",ident1," vs ",ident2))
    LCvsT1_cats = FindMarkers(bigObj,ident.1 = ident1,ident.2 = ident2,group.by = "CellType_AcuteCov_TP")
    LCvsT1_cats = subset(LCvsT1_cats, LCvsT1_cats$p_val_adj < 0.05)
    if(nrow(LCvsT1_cats) > 0){
      LCvsT1_cats$Comparison = "LCVsAcute_severe"
      dge_vsT1_PerTP_perCT_perCat[[paste0(i,"_",j,"_severe")]] = LCvsT1_cats
    }
    
    ident2 = paste0(i,"_R_Recov")
    print(paste0("Comparison ",ident1," vs ", ident2))
    LCvsR_cats = FindMarkers(bigObj,ident.1 = ident1,ident.2 = ident2,group.by = "CellType_AcuteCov_TP")
    LCvsR_cats = subset(LCvsR_cats, LCvsR_cats$p_val_adj < 0.05)
    if(nrow(LCvsR_cats) > 0){
      LCvsR_cats$Comparison = "LCVsRecov_severe"
      dge_vsR_PerTP_perCT_perCat[[paste0(i,"_",j,"_severe")]] = LCvsR_cats
    }
    
    
    ident2 = paste0(i,"_HC_HC")
    print(paste0("Comparison ",ident1," vs ",ident2))
    LCvsHC_cats = FindMarkers(bigObj,ident.1 = ident1,ident.2 = ident2,group.by = "CellType_AcuteCov_TP")
    LCvsHC_cats = subset(LCvsHC_cats, LCvsHC_cats$p_val_adj < 0.05)
    if(nrow(LCvsHC_cats) > 0){
      LCvsHC_cats$Comparison = "LCVsHC_severe"
      dge_vsHC_PerTP_perCT_perCat[[paste0(i,"_",j,"_severe")]] = LCvsHC_cats
    
    }}}

for (i in celltypes){
  ident.1 = paste0(i,"_T1")
  ident.2 = paste0(i,"_HC_HC")
  T1vsHC = FindMarkers(bigObj,ident.1 = ident1,ident.2 = ident2,group.by = "CellType_AcuteCov_TP")
  T1vsHC = subset(T1vsHC, T1vsHC$p_val_adj < 0.05)
  if(nrow(T1vsHC) > 0){
    T1vsHC$Comparison = "T1vsHC"
    dge_T1vsHC_perCT[[paste0(i,"_T1vsHC")]] = T1vsHC
  }
}


dge_vsR_PerTP_perCT_perCat = lapply(dge_vsR_PerTP_perCT_perCat,tibble::rownames_to_column)
dge_vsT1_PerTP_perCT_perCat = lapply(dge_vsT1_PerTP_perCT_perCat,tibble::rownames_to_column)
dge_vsHC_PerTP_perCT_perCat = lapply(dge_vsHC_PerTP_perCT_perCat,tibble::rownames_to_column)
dge_T1vsHC_perCT = lapply(dge_T1vsHC_perCT,tibble::rownames_to_column)

dge_vsR_PerTP_perCT_perCat_up <- lapply(dge_vsR_PerTP_perCT_perCat, subset, avg_log2FC > 0.8)
dge_vsHC_PerTP_perCT_perCat_up <- lapply(dge_vsHC_PerTP_perCT_perCat, subset, avg_log2FC > 0.8)
dge_vsT1_PerTP_perCT_perCat_up <- lapply(dge_vsT1_PerTP_perCT_perCat, subset, avg_log2FC > 0.8)

dge_T1vsHC_perCT_subset = lapply(dge_T1vsHC_perCT, subset, abs(avg_log2FC) > 0.8)
dge_vsR_PerTP_perCT_perCat_subset = lapply(dge_vsR_PerTP_perCT_perCat,subset, abs(avg_log2FC) > 0.8)
dge_vsT1_PerTP_perCT_perCat_subset = lapply(dge_vsT1_PerTP_perCT_perCat,subset, abs(avg_log2FC) > 0.8)

dge_vsR_PerTP_perCT_perCat_up <-lapply(dge_vsR_PerTP_perCT_perCat_up,function(x) x[!grepl("^MT-",x$rowname),])
dge_vsHC_PerTP_perCT_perCat_up <-lapply(dge_vsHC_PerTP_perCT_perCat_up,function(x) x[!grepl("^MT-",x$rowname),])
dge_vsT1_PerTP_perCT_perCat_up <-lapply(dge_vsT1_PerTP_perCT_perCat_up,function(x) x[!grepl("^MT-",x$rowname),])
dge_T1vsHC_perCT_subset = lapply(dge_T1vsHC_perCT_subset,function(x) x[!grepl("^MT-",x$rowname),])

dge_vsR_PerTP_perCT_perCat_down <- lapply(dge_vsR_PerTP_perCT_perCat, subset, avg_log2FC < -0.8)
dge_vsHC_PerTP_perCT_perCat_down <- lapply(dge_vsHC_PerTP_perCT_perCat, subset, avg_log2FC < -0.8)
dge_vsT1_PerTP_perCT_perCat_down <- lapply(dge_vsT1_PerTP_perCT_perCat, subset, avg_log2FC < -0.8)

dge_vsR_PerTP_perCT_perCat_down <-lapply(dge_vsR_PerTP_perCT_perCat_down,function(x) x[!grepl("^MT-",x$rowname),])
dge_vsHC_PerTP_perCT_perCat_down <-lapply(dge_vsHC_PerTP_perCT_perCat_down,function(x) x[!grepl("^MT-",x$rowname),])
dge_vsT1_PerTP_perCT_perCat_down <-lapply(dge_vsT1_PerTP_perCT_perCat_down,function(x) x[!grepl("^MT-",x$rowname),])


cd14Mono_nonsevereGenes_LC = Reduce(union,list(dge_vsR_PerTP_perCT_perCat_up$`CD14 Monocytes_T2_nonsevere`$rowname,
                                            dge_vsR_PerTP_perCT_perCat_up$`CD14 Monocytes_T3_nonsevere`$rowname,
                                            dge_vsR_PerTP_perCT_perCat_up$`CD14 Monocytes_T4_nonsevere`$rowname,
                                            dge_vsT1_PerTP_perCT_perCat_up$`CD14 Monocytes_T2_nonsevere`$rowname,
                                            dge_vsT1_PerTP_perCT_perCat_up$`CD14 Monocytes_T3_nonsevere`$rowname,
                                            dge_vsT1_PerTP_perCT_perCat_up$`CD14 Monocytes_T4_nonsevere`$rowname,
                                            dge_vsT1_PerTP_perCT_perCat_down$`CD14 Monocytes_T2_nonsevere`$rowname,
                                            dge_vsT1_PerTP_perCT_perCat_down$`CD14 Monocytes_T3_nonsevere`$rowname,
                                            dge_vsT1_PerTP_perCT_perCat_down$`CD14 Monocytes_T4_nonsevere`$rowname,
                                            dge_vsT1_PerTP_perCT_perCat_down$`CD14 Monocytes_T5_nonsevere`$rowname,
                                            dge_T1vsHC_perCT_subset$`CD14 Monocytes_T1vsHC`$rowname))


cd14Mono_nonsevereGenes_LC = cd14Mono_nonsevereGenes_LC[!grepl("MT-",cd14Mono_nonsevereGenes_LC)]
cd14Mono_nonsevereGenes_LC = cd14Mono_nonsevereGenes_LC[!grepl("IGL",cd14Mono_nonsevereGenes_LC)]
cd14Mono_nonsevereGenes_LC = cd14Mono_nonsevereGenes_LC[!grepl("IGH",cd14Mono_nonsevereGenes_LC)]
cd14Mono_nonsevereGenes_LC = cd14Mono_nonsevereGenes_LC[!grepl("^AC[0-9].+",cd14Mono_nonsevereGenes_LC)]
cd14Mono_nonsevereGenes_LC = cd14Mono_nonsevereGenes_LC[!grepl("^AL[0-9].+",cd14Mono_nonsevereGenes_LC)]
cd14Mono_nonsevereGenes_LC = cd14Mono_nonsevereGenes_LC[!grepl("^AP[0-9].+",cd14Mono_nonsevereGenes_LC)]
cd14Mono_nonsevereGenes_LC = cd14Mono_nonsevereGenes_LC[!grepl("^RPL",cd14Mono_nonsevereGenes_LC)]
cd14Mono_nonsevereGenes_LC = cd14Mono_nonsevereGenes_LC[!grepl("^RPS",cd14Mono_nonsevereGenes_LC)]
cd14Mono_nonsevereGenes_LC = cd14Mono_nonsevereGenes_LC[!grepl("^IGK",cd14Mono_nonsevereGenes_LC)]
cd14Mono_nonsevereGenes_LC = cd14Mono_nonsevereGenes_LC[!grepl("-AS1",cd14Mono_nonsevereGenes_LC)]



####### 
cd14Mono = subset(bigObj,subset = CellType == "CD14 Monocytes")
cd14Mono$AcuteCovidStatus_TimePoint2 = paste0(cd14Mono$AcuteCovidStatus,"_",cd14Mono$TimePointAll)
cd14Mono$AcuteCovidStatus_TimePoint2 = ifelse(cd14Mono$AcuteCovidStatus_TimePoint2 %in% c("R_T3","R_T4","R_T5"), "R_Recov", as.character(cd14Mono$AcuteCovidStatus_TimePoint2))
MonogeneMat = AverageExpression(cd14Mono,assay = "RNA",features = cd14Mono_nonsevereGenes_LC,group.by = "AcuteCovidStatus_TimePoint2")
MonogeneMat = as.matrix(MonogeneMat$RNA)

MonogeneMat = MonogeneMat[,c(1,12,7,2,8,9,10,11,3,4,5,6)]

scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

MonogeneMat_scaled = scale_rows(MonogeneMat)

#########
cd16Mono_nonsevereGenes_LC = Reduce(union,list(dge_vsR_PerTP_perCT_perCat_up$`CD16 Monocytes_T2_nonsevere`$rowname,
                                            dge_vsR_PerTP_perCT_perCat_up$`CD16 Monocytes_T3_nonsevere`$rowname,
                                            dge_vsR_PerTP_perCT_perCat_up$`CD16 Monocytes_T4_nonsevere`$rowname,
                                            dge_vsT1_PerTP_perCT_perCat_up$`CD16 Monocytes_T2_nonsevere`$rowname,
                                            dge_vsT1_PerTP_perCT_perCat_up$`CD16 Monocytes_T3_nonsevere`$rowname,
                                            dge_vsT1_PerTP_perCT_perCat_up$`CD16 Monocytes_T4_nonsevere`$rowname,
                                            dge_vsT1_PerTP_perCT_perCat_down$`CD16 Monocytes_T2_nonsevere`$rowname,
                                            dge_vsT1_PerTP_perCT_perCat_down$`CD16 Monocytes_T3_nonsevere`$rowname,
                                            dge_vsT1_PerTP_perCT_perCat_down$`CD16 Monocytes_T4_nonsevere`$rowname,
                                            dge_vsT1_PerTP_perCT_perCat_down$`CD16 Monocytes_T5_nonsevere`$rowname,
                                            dge_vsHC_PerTP_perCT_perCat_up$`CD16 Monocytes_T2_nonsevere`$rowname,
                                            dge_vsHC_PerTP_perCT_perCat_up$`CD16 Monocytes_T3_nonsevere`$rowname,
                                            dge_vsHC_PerTP_perCT_perCat_up$`CD16 Monocytes_T4_nonsevere`$rowname))



cd16Mono_nonsevereGenes_LC = cd16Mono_nonsevereGenes_LC[!grepl("IGL",cd16Mono_nonsevereGenes_LC)]
cd16Mono_nonsevereGenes_LC = cd16Mono_nonsevereGenes_LC[!grepl("IGH",cd16Mono_nonsevereGenes_LC)]
cd16Mono_nonsevereGenes_LC = cd16Mono_nonsevereGenes_LC[!grepl("^AC[0-9].+",cd16Mono_nonsevereGenes_LC)]
cd16Mono_nonsevereGenes_LC = cd16Mono_nonsevereGenes_LC[!grepl("^AL[0-9].+",cd16Mono_nonsevereGenes_LC)]
cd16Mono_nonsevereGenes_LC = cd16Mono_nonsevereGenes_LC[!grepl("^AP[0-9].+",cd16Mono_nonsevereGenes_LC)]
cd16Mono_nonsevereGenes_LC = cd16Mono_nonsevereGenes_LC[!grepl("^RPL",cd16Mono_nonsevereGenes_LC)]
cd16Mono_nonsevereGenes_LC = cd16Mono_nonsevereGenes_LC[!grepl("^RPS",cd16Mono_nonsevereGenes_LC)]
cd16Mono_nonsevereGenes_LC = cd16Mono_nonsevereGenes_LC[!grepl("^IGK",cd16Mono_nonsevereGenes_LC)]
cd16Mono_nonsevereGenes_LC = cd16Mono_nonsevereGenes_LC[!grepl("-AS1",cd16Mono_nonsevereGenes_LC)]

cd16Mono = subset(bigObj,subset = CellType == "CD16 Monocytes")
cd16Mono$AcuteCovidStatus_TimePoint2 = paste0(cd16Mono$AcuteCovidStatus,"_",cd16Mono$TimePointAll)
cd16Mono$AcuteCovidStatus_TimePoint2 = ifelse(cd16Mono$AcuteCovidStatus_TimePoint2 %in% c("R_T3","R_T4","R_T5"), "R_Recov", as.character(cd16Mono$AcuteCovidStatus_TimePoint2))
cd16MonogeneMat = AverageExpression(cd16Mono,assay = "RNA",features = cd16Mono_nonsevereGenes_LC,group.by = "AcuteCovidStatus_TimePoint2")
cd16MonogeneMat = as.matrix(cd16MonogeneMat$RNA)

cd16MonogeneMat = cd16MonogeneMat[,c(1,12,7,2,8,9,10,11,3,4,5,6)]

cd16MonogeneMat_scaled = scale_rows(cd16MonogeneMat)

##########
nkcells_nonsevereGenes_LC = Reduce(union,list(dge_vsR_PerTP_perCT_perCat_up$`NK cells_T2_nonsevere`$rowname,
                                           dge_vsR_PerTP_perCT_perCat_up$`NK cells_T3_nonsevere`$rowname,
                                           dge_vsR_PerTP_perCT_perCat_up$`NK cells_T4_nonsevere`$rowname,
                                           dge_vsT1_PerTP_perCT_perCat_up$`NK cells_T2_nonsevere`$rowname,
                                           dge_vsT1_PerTP_perCT_perCat_up$`NK cells_T3_nonsevere`$rowname,
                                           dge_vsT1_PerTP_perCT_perCat_up$`NK cells_T4_nonsevere`$rowname,
                                           dge_vsT1_PerTP_perCT_perCat_down$`NK cells_T2_nonsevere`$rowname,
                                           dge_vsT1_PerTP_perCT_perCat_down$`NK cells_T3_nonsevere`$rowname,
                                           dge_vsT1_PerTP_perCT_perCat_down$`NK cells_T4_nonsevere`$rowname,
                                           dge_vsT1_PerTP_perCT_perCat_down$`NK cells_T5_nonsevere`$rowname,
                                           dge_vsHC_PerTP_perCT_perCat_up$`NK cells_T2_nonsevere`$rowname,
                                           dge_vsHC_PerTP_perCT_perCat_up$`NK cells_T3_nonsevere`$rowname,
                                           dge_vsHC_PerTP_perCT_perCat_up$`NK cells_T4_nonsevere`$rowname))


nkcells_nonsevereGenes_LC = nkcells_nonsevereGenes_LC[!grepl("IGL",nkcells_nonsevereGenes_LC)]
nkcells_nonsevereGenes_LC = nkcells_nonsevereGenes_LC[!grepl("IGH",nkcells_nonsevereGenes_LC)]
nkcells_nonsevereGenes_LC = nkcells_nonsevereGenes_LC[!grepl("^AC[0-9].+",nkcells_nonsevereGenes_LC)]
nkcells_nonsevereGenes_LC = nkcells_nonsevereGenes_LC[!grepl("^AL[0-9].+",nkcells_nonsevereGenes_LC)]
nkcells_nonsevereGenes_LC = nkcells_nonsevereGenes_LC[!grepl("^AP[0-9].+",nkcells_nonsevereGenes_LC)]
nkcells_nonsevereGenes_LC = nkcells_nonsevereGenes_LC[!grepl("^RPL",nkcells_nonsevereGenes_LC)]
nkcells_nonsevereGenes_LC = nkcells_nonsevereGenes_LC[!grepl("^RPS",nkcells_nonsevereGenes_LC)]
nkcells_nonsevereGenes_LC = nkcells_nonsevereGenes_LC[!grepl("^IGK",nkcells_nonsevereGenes_LC)]
nkcells_nonsevereGenes_LC = nkcells_nonsevereGenes_LC[!grepl("-AS1",nkcells_nonsevereGenes_LC)]

nkcells = subset(bigObj,subset = CellType == "NK cells")
nkcells$AcuteCovidStatus_TimePoint2 = paste0(nkcells$AcuteCovidStatus,"_",nkcells$TimePointAll)
nkcells$AcuteCovidStatus_TimePoint2 = ifelse(nkcells$AcuteCovidStatus_TimePoint2 %in% c("R_T3","R_T4","R_T5"), "R_Recov", as.character(nkcells$AcuteCovidStatus_TimePoint2))
nkcellsgeneMat = AverageExpression(nkcells,assay = "RNA",features = nkcells_nonsevereGenes_LC,group.by = "AcuteCovidStatus_TimePoint2")
nkcellsgeneMat = as.matrix(nkcellsgeneMat$RNA)

nkcellsgeneMat = nkcellsgeneMat[,c(1,12,7,2,8,9,10,11,3,4,5,6)]

nkcellsgeneMat_scaled = scale_rows(nkcellsgeneMat)
nkMarkersForPlotting = c("BCOR","TMEME67","CRTAM","FEZ1","CFAP20","KDM6B","IFFO2","GCH1","NFKBIA",
                         "PPP1R15A","CD83","MLF1","TSIX","HES4","LDLR","U2AF1","VIM","HDAC5","PAK1","SREBF1","ZC3H12A","ITGA5",
                         "MAN2B1","RFX1","ECE1","EVI5","DENND5A","GAS7","MICAL3","TAGLN2","TNIP1","HLA-DQB1","MYADM")


#######
bcells_nonsevereGenes_LC = Reduce(union,list(dge_vsR_PerTP_perCT_perCat_up$`B cells_T2_nonsevere`$rowname,
                                          dge_vsR_PerTP_perCT_perCat_up$`B cells_T3_nonsevere`$rowname,
                                          dge_vsR_PerTP_perCT_perCat_up$`B cells_T4_nonsevere`$rowname,
                                          dge_vsT1_PerTP_perCT_perCat_up$`B cells_T2_nonsevere`$rowname,
                                          dge_vsT1_PerTP_perCT_perCat_up$`B cells_T3_nonsevere`$rowname,
                                          dge_vsT1_PerTP_perCT_perCat_up$`B cells_T4_nonsevere`$rowname,
                                          dge_vsT1_PerTP_perCT_perCat_down$`B cells_T2_nonsevere`$rowname,
                                          dge_vsT1_PerTP_perCT_perCat_down$`B cells_T3_nonsevere`$rowname,
                                          dge_vsT1_PerTP_perCT_perCat_down$`B cells_T4_nonsevere`$rowname,
                                          dge_vsT1_PerTP_perCT_perCat_down$`B cells_T5_nonsevere`$rowname,
                                          dge_vsHC_PerTP_perCT_perCat_up$`B cells_T2_nonsevere`$rowname,
                                          dge_vsHC_PerTP_perCT_perCat_up$`B cells_T3_nonsevere`$rowname,
                                          dge_vsHC_PerTP_perCT_perCat_up$`B cells_T4_nonsevere`$rowname))


bcells_nonsevereGenes_LC = bcells_nonsevereGenes_LC[!grepl("^AC[0-9].+",bcells_nonsevereGenes_LC)]
bcells_nonsevereGenes_LC = bcells_nonsevereGenes_LC[!grepl("^AL[0-9].+",bcells_nonsevereGenes_LC)]
bcells_nonsevereGenes_LC = bcells_nonsevereGenes_LC[!grepl("^RPL",bcells_nonsevereGenes_LC)]
bcells_nonsevereGenes_LC = bcells_nonsevereGenes_LC[!grepl("^RPS",bcells_nonsevereGenes_LC)]
bcells_nonsevereGenes_LC = bcells_nonsevereGenes_LC[!grepl("^AP[0-9].+",bcells_nonsevereGenes_LC)]
bcells_nonsevereGenes_LC = bcells_nonsevereGenes_LC[!grepl("-AS1",bcells_nonsevereGenes_LC)]


bcells = subset(bigObj,subset = CellType == "B cells")
bcells$AcuteCovidStatus_TimePoint2 = paste0(bcells$AcuteCovidStatus,"_",bcells$TimePointAll)
bcells$AcuteCovidStatus_TimePoint2 = ifelse(bcells$AcuteCovidStatus_TimePoint2 %in% c("R_T3","R_T4","R_T5"), "R_Recov", as.character(bcells$AcuteCovidStatus_TimePoint2))
bcellsgeneMat = AverageExpression(bcells,assay = "RNA",features = bcells_nonsevereGenes_LC,group.by = "AcuteCovidStatus_TimePoint2")
bcellsgeneMat = as.matrix(bcellsgeneMat$RNA)

bcellsgeneMat = bcellsgeneMat[,c(1,12,7,2,8,9,10,11,3,4,5,6)]

bcellsgeneMat_scaled = scale_rows(bcellsgeneMat)

#########
cd8T_nonsevereGenes_LC = Reduce(union,list(dge_vsR_PerTP_perCT_perCat_up$`CD8 T cells_T2_nonsevere`$rowname,
                                        dge_vsR_PerTP_perCT_perCat_up$`CD8 T cells_T3_nonsevere`$rowname,
                                        dge_vsR_PerTP_perCT_perCat_up$`CD8 T cells_T4_nonsevere`$rowname,
                                        dge_vsT1_PerTP_perCT_perCat_up$`CD8 T cells_T2_nonsevere`$rowname,
                                        dge_vsT1_PerTP_perCT_perCat_up$`CD8 T cells_T3_nonsevere`$rowname,
                                        dge_vsT1_PerTP_perCT_perCat_up$`CD8 T cells_T4_nonsevere`$rowname,
                                        dge_vsT1_PerTP_perCT_perCat_down$`CD8 T cells_T2_nonsevere`$rowname,
                                        dge_vsT1_PerTP_perCT_perCat_down$`CD8 T cells_T3_nonsevere`$rowname,
                                        dge_vsT1_PerTP_perCT_perCat_down$`CD8 T cells_T4_nonsevere`$rowname,
                                        dge_vsT1_PerTP_perCT_perCat_down$`CD8 T cells_T5_nonsevere`$rowname,
                                        dge_vsHC_PerTP_perCT_perCat_up$`CD8 T cells_T2_nonsevere`$rowname,
                                        dge_vsHC_PerTP_perCT_perCat_up$`CD8 T cells_T3_nonsevere`$rowname,
                                        dge_vsHC_PerTP_perCT_perCat_up$`CD8 T cells_T4_nonsevere`$rowname))


cd8T_nonsevereGenes_LC = cd8T_nonsevereGenes_LC[!grepl("IGL",cd8T_nonsevereGenes_LC)]
cd8T_nonsevereGenes_LC = cd8T_nonsevereGenes_LC[!grepl("IGH",cd8T_nonsevereGenes_LC)]
cd8T_nonsevereGenes_LC = cd8T_nonsevereGenes_LC[!grepl("^AC[0-9].+",cd8T_nonsevereGenes_LC)]
cd8T_nonsevereGenes_LC = cd8T_nonsevereGenes_LC[!grepl("^AL[0-9].+",cd8T_nonsevereGenes_LC)]
cd8T_nonsevereGenes_LC = cd8T_nonsevereGenes_LC[!grepl("^RPL",cd8T_nonsevereGenes_LC)]
cd8T_nonsevereGenes_LC = cd8T_nonsevereGenes_LC[!grepl("^RPS",cd8T_nonsevereGenes_LC)]
cd8T_nonsevereGenes_LC = cd8T_nonsevereGenes_LC[!grepl("^IGK",cd8T_nonsevereGenes_LC)]
cd8T_nonsevereGenes_LC = cd8T_nonsevereGenes_LC[!grepl("^AP[0-9].+",cd8T_nonsevereGenes_LC)]
cd8T_nonsevereGenes_LC = cd8T_nonsevereGenes_LC[!grepl("-AS1",cd8T_nonsevereGenes_LC)]


cd8T = subset(bigObj,subset = CellType == "CD8 T cells")
cd8T$AcuteCovidStatus_TimePoint2 = paste0(cd8T$AcuteCovidStatus,"_",cd8T$TimePointAll)
cd8T$AcuteCovidStatus_TimePoint2 = ifelse(cd8T$AcuteCovidStatus_TimePoint2 %in% c("R_T3","R_T4","R_T5"), "R_Recov", as.character(cd8T$AcuteCovidStatus_TimePoint2))
cd8TgeneMat = AverageExpression(cd8T,assay = "RNA",features = cd8T_nonsevereGenes_LC,group.by = "AcuteCovidStatus_TimePoint2")
cd8TgeneMat = as.matrix(cd8TgeneMat$RNA)

cd8TgeneMat = cd8TgeneMat[,c(1,12,7,2,8,9,10,11,3,4,5,6)]

cd8TgeneMat_scaled = scale_rows(cd8TgeneMat)

#######

cd4T_nonsevereGenes_LC = Reduce(union,list(dge_vsR_PerTP_perCT_perCat_up$`CD4 T cells_T2_nonsevere`$rowname,
                                        dge_vsR_PerTP_perCT_perCat_up$`CD4 T cells_T3_nonsevere`$rowname,
                                        dge_vsR_PerTP_perCT_perCat_up$`CD4 T cells_T4_nonsevere`$rowname,
                                        dge_vsT1_PerTP_perCT_perCat_up$`CD4 T cells_T2_nonsevere`$rowname,
                                        dge_vsT1_PerTP_perCT_perCat_up$`CD4 T cells_T3_nonsevere`$rowname,
                                        dge_vsT1_PerTP_perCT_perCat_up$`CD4 T cells_T4_nonsevere`$rowname,
                                        dge_vsT1_PerTP_perCT_perCat_down$`CD4 T cells_T2_nonsevere`$rowname,
                                        dge_vsT1_PerTP_perCT_perCat_down$`CD4 T cells_T3_nonsevere`$rowname,
                                        dge_vsT1_PerTP_perCT_perCat_down$`CD4 T cells_T4_nonsevere`$rowname,
                                        dge_vsT1_PerTP_perCT_perCat_down$`CD4 T cells_T5_nonsevere`$rowname,
                                        dge_vsHC_PerTP_perCT_perCat_up$`CD4 T cells_T2_nonsevere`$rowname,
                                        dge_vsHC_PerTP_perCT_perCat_up$`CD4 T cells_T3_nonsevere`$rowname,
                                        dge_vsHC_PerTP_perCT_perCat_up$`CD4 T cells_T4_nonsevere`$rowname))


cd4T_nonsevereGenes_LC = cd4T_nonsevereGenes_LC[!grepl("IGL",cd4T_nonsevereGenes_LC)]
cd4T_nonsevereGenes_LC = cd4T_nonsevereGenes_LC[!grepl("IGH",cd4T_nonsevereGenes_LC)]
cd4T_nonsevereGenes_LC = cd4T_nonsevereGenes_LC[!grepl("^AC[0-9].+",cd4T_nonsevereGenes_LC)]
cd4T_nonsevereGenes_LC = cd4T_nonsevereGenes_LC[!grepl("^AL[0-9].+",cd4T_nonsevereGenes_LC)]
cd4T_nonsevereGenes_LC = cd4T_nonsevereGenes_LC[!grepl("^RPL",cd4T_nonsevereGenes_LC)]
cd4T_nonsevereGenes_LC = cd4T_nonsevereGenes_LC[!grepl("^RPS",cd4T_nonsevereGenes_LC)]
cd4T_nonsevereGenes_LC = cd4T_nonsevereGenes_LC[!grepl("^IGK",cd4T_nonsevereGenes_LC)]
cd4T_nonsevereGenes_LC = cd4T_nonsevereGenes_LC[!grepl("^AP[0-9].+",cd4T_nonsevereGenes_LC)]
cd4T_nonsevereGenes_LC = cd4T_nonsevereGenes_LC[!grepl("-AS1",cd4T_nonsevereGenes_LC)]

cd4T = subset(bigObj,subset = CellType == "CD4 T cells")
cd4T$AcuteCovidStatus_TimePoint2 = paste0(cd4T$AcuteCovidStatus,"_",cd4T$TimePointAll)
cd4T$AcuteCovidStatus_TimePoint2 = ifelse(cd4T$AcuteCovidStatus_TimePoint2 %in% c("R_T3","R_T4","R_T5"), "R_Recov", as.character(cd4T$AcuteCovidStatus_TimePoint2))
cd4TgeneMat = AverageExpression(cd4T,assay = "RNA",features = cd4T_nonsevereGenes_LC,group.by = "AcuteCovidStatus_TimePoint2")
cd4TgeneMat = as.matrix(cd4TgeneMat$RNA)

cd4TgeneMat = cd4TgeneMat[,c(1,12,7,2,8,9,10,11,3,4,5,6)]

cd4TgeneMat_scaled = scale_rows(cd4TgeneMat)


