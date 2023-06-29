library(Seurat)
library(tidyverse)


#HUMAN DATA
EdKaZhou = readRDS("EdKaZhouHypoNeurons_mt10_integrated.rds")
EdKaZhou@meta.data$Dataset = ifelse(EdKaZhou@meta.data$Age %in% c(29, 42, 50), "Siletti", "Zhou")
EdKaZhou@meta.data$Dataset = ifelse(EdKaZhou@meta.data$sample %in% c("CS22_2_hypo", "CS22_hypo", "GW16_hypo", "GW18_hypo", "GW20_34_hypo", "GW22T_hypo1", "GW25_3V_hypo", "GW19_hypo"), "Herb", EdKaZhou@meta.data$Dataset)
DefaultAssay(EdKaZhou) = "RNA"
EdKaZhouNeurons = SplitObject(EdKaZhou, split.by = "Dataset")


#HYPOMAP DATA
hypoMap = readRDS("~/hypoMap.rds")
hypoMap = ConvertGeneNames(hypoMap, reference.names = row.names(Neuro.combined), homolog.table = "https://seurat.nygenome.org/azimuth/references/homologs.rds")
#saveRDS(hypoMap, "~/hypoMapHuman.rds")
#hypoMap = readRDS("hypoMapHuman.rds")
Idents(hypoMap) = "Dataset"
hypoMap = subset(hypoMap, idents = c("Dowsett10xnuc", "RomanovDev10x", "KimDev10x"), invert=T)
Idents(hypoMap) = "C2_named"
hypoMap = subset(hypoMap, idents = c("C2-1: Neurons"))
hypoMapData = SplitObject(hypoMap, split.by = "Dataset")

#KIM DATA
Kim_Meta = read.csv("~/Downloads/GSE132355_E10-P45_umap_data.csv")
Kim_MetaNeurons = subset(Kim_Meta, Kim_Meta$Cluster %in% c("GE/POA", "VMH", "POA/SCN", "SMN", "PVH & SON", "MMN", "PMN", "LH", "ARC", "POA", "SCN", "DMH/ID", "Hypothalamic neurons"))
Kim_Data = readRDS("~/Downloads/GSE132355_E10-P45_log2cpm.rds")
Kim_Data = CreateSeuratObject(Kim_Data, project= "KimDev")
Kim_Data = subset(Kim_Data, cells = Kim_MetaNeurons$X.1)
DefaultAssay(Kim_Data) = "RNA"
Kim_Data = ConvertGeneNames(Kim_Data, reference.names = row.names(EdKaZhou), homolog.table = "https://seurat.nygenome.org/azimuth/references/homologs.rds")
Kim_Data@meta.data$Dataset = "KimDev"
#saveRDS(Kim_Data, "~/KimDevNeuronsHuman.rds")
#KimDev = readRDS("KimDevNeuronsHuman.rds")


#ROMANOV DATA
Romanov_Data = readRDS("~/Downloads/GSE132730_TractNAE_integrated.rds")
Idents(Romanov_Data) = "Class"
Romanov_Data = subset(Romanov_Data, idents = "Neurons")
DefaultAssay(Romanov_Data) = "RNA"
Romanov_Data = ConvertGeneNames(Romanov_Data, reference.names = row.names(EdKaZhou), homolog.table = "https://seurat.nygenome.org/azimuth/references/homologs.rds")
Romanov_Data@meta.data$Dataset = "RomanovDev"
#saveRDS(Romanov_Data, "~/RomanovDevNeuronsHuman.rds")
#RomanovDev = readRDS("RomanovDevNeuronsHuman.rds")


#TRANSFER HYPOMAP DATA TO HUMAN
Hypo.anchors <- FindTransferAnchors(reference = hypoMap, query = EdKaZhou,  dims = 1:30)
predictions <- TransferData(anchorset = Hypo.anchors, refdata = list("Dataset" = hypoMap$Dataset, "C2_named" = hypoMap$C2_named, "C7_named" = hypoMap$C7_named, "C25_named" = hypoMap$C25_named, "C66_named" = hypoMap$C66_named, "C185_named" = hypoMap$C185_named, "C286_named" = hypoMap$C286_named, "C465_named" = hypoMap$C465_named, "Region_predicted" = hypoMap$Region_predicted, "Region_summarized" = hypoMap$Region_summarized),   dims = 1:30)
EdKaZhou <- AddMetaData(EdKaZhou, metadata = predictions)
saveRDS(predictions, "predictions_refined.rds")


Hypo.list = c("Romanov" = RomanovDev, "Kim" = KimDev, EdKaZhouNeurons, hypoMapData)
rm(hypoMap)
rm(EdKaZhou)

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

#saveRDS(Hypo.combined, "~/HypoMapHUMANDev.rds")

Hypo.combined = ScaleData(Hypo.combined, vars.to.regress = "percent.mt")
Hypo.combined = RunPCA(Hypo.combined)
Hypo.combined = RunUMAP(Hypo.combined, dims = 1:30)
saveRDS(Hypo.combined, "~/HypoMapHUMANDev.rds")