library(Seurat)
library(tidyverse)

RomanovDev = readRDS("RomanovDevNeuronsHuman.rds")
KimDev = readRDS("KimDevNeuronsHuman.rds")
EdKaZhou = readRDS("EdKaZhouHypoNeurons_mt10_integrated.rds")
EdKaZhou@meta.data$Dataset = ifelse(EdKaZhou@meta.data$Age %in% c(29, 42, 50), "Siletti", "Zhou")
EdKaZhou@meta.data$Dataset = ifelse(EdKaZhou@meta.data$sample %in% c("CS22_2_hypo", "CS22_hypo", "GW16_hypo", "GW18_hypo", "GW20_34_hypo", "GW22T_hypo1", "GW25_3V_hypo", "GW19_hypo"), "Herb", EdKaZhou@meta.data$Dataset)
EdKaZhouNeurons = SplitObject(EdKaZhou, split.by = "Dataset")
rm(EdKaZhou)
hypoMap = readRDS("hypoMapHuman.rds")
#hypoMap = readRDS("~/Dropbox/Columbia/Sci Adv/hypoMapHuman.rds")
Idents(hypoMap) = "Dataset"
hypoMap = subset(hypoMap, idents = c("Dowsett10xnuc", "RomanovDev10x", "KimDev10x"), invert=T)
hypoMapData = SplitObject(hypoMap, split.by = "Dataset")
rm(hypoMap)

Hypo.list = c("Romanov" = RomanovDev, "Kim" = KimDev, EdKaZhouNeurons, hypoMapData)

for(n in names(Hypo.list)){
  DefaultAssay(Hypo.list[[n]]) = "RNA"
  Hypo.list[[n]]  = DietSeurat(Hypo.list[[n]])
  MT.genes = grep(pattern="^MT-", x = row.names(Hypo.list[[n]]), value=T)
  percent.mt = Matrix::colSums(Hypo.list[[n]]@assays[["RNA"]][MT.genes, ])/Matrix::colSums(Hypo.list[[n]]@assays[["RNA"]])
  Hypo.list[[n]] = AddMetaData(Hypo.list[[n]], percent.mt, "percent.mt")
  Hypo.list[[n]] <- NormalizeData(Hypo.list[[n]])
  Hypo.list[[n]] <- FindVariableFeatures(Hypo.list[[n]], selection.method = "vst", nfeatures = 2000)
}
feat <- SelectIntegrationFeatures(object.list = Hypo.list)

for(n in names(Hypo.list)){
  Hypo.list[[n]] = ScaleData(Hypo.list[[n]], features = feat, verbose = FALSE)
  Hypo.list[[n]] = RunPCA(Hypo.list[[n]], features = feat, verbose = FALSE)
}

Hypo.anchors <- FindIntegrationAnchors(object.list = Hypo.list, anchor.features = feat, reduction = "rpca")
Hypo.combined <- IntegrateData(anchorset = Hypo.anchors)

#saveRDS(Hypo.combined, "~/HypoMapHUMAN.rds")

Hypo.combined = ScaleData(Hypo.combined, vars.to.regress = "percent.mt")
Hypo.combined = RunPCA(Hypo.combined)
Hypo.combined = RunUMAP(Hypo.combined, dims = 1:30)
saveRDS(Hypo.combined, "~/HypoMapHUMANDev.rds")