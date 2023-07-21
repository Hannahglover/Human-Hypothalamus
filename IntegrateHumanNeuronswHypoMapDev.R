library(Seurat)
library(tidyverse)


#HUMAN DATA
EdKaZhou = readRDS("~/Dropbox/LabMac/EdKaZhouHypoNeurons_mt10_integrated.rds")
EdKaZhou@meta.data$Dataset = ifelse(EdKaZhou@meta.data$Age %in% c(29, 42, 50), "Siletti", "Zhou")
EdKaZhou@meta.data$Dataset = ifelse(EdKaZhou@meta.data$sample %in% c("CS22_2_hypo", "CS22_hypo", "GW16_hypo", "GW18_hypo", "GW20_34_hypo", "GW22T_hypo1", "GW25_3V_hypo", "GW19_hypo"), "Herb", EdKaZhou@meta.data$Dataset)
DefaultAssay(EdKaZhou) = "RNA"
EdKaZhouNeurons = SplitObject(EdKaZhou, split.by = "Dataset")

#HYPOMAP DATA
hypoMap = readRDS("hypoMap.rds")
hypoMap = ConvertGeneNames(hypoMap, reference.names = row.names(Neuro.combined), homolog.table = "https://seurat.nygenome.org/azimuth/references/homologs.rds")
#saveRDS(hypoMap, "~/hypoMapHuman.rds")
#hypoMap = readRDS("~/Dropbox/Columbia/Sci Adv/hypoMapHuman.rds")
#hypoMap = readRDS("~/Dropbox/Columbia/Sci Adv/hypoMapHuman.rds")
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
#RomanovDev = readRDS("~/RomanovDevNeuronsHuman.rds")

Hypo.list = c("Romanov" = RomanovDev, "Kim" = KimDev, EdKaZhouNeurons, hypoMapData)

CommonGenes = intersect(row.names(RomanovDev@assays$RNA@data), row.names(KimDev@assays$RNA@data))
CommonGenes = intersect(CommonGenes, row.names(EdKaZhouNeurons@assays$RNA@data))
CommonGenes = intersect(CommonGenes, row.names(hypoMapData@assays$RNA@data))
write.csv(as.data.frame(CommonGenes), "CommonGenes_Fig3.csv")
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
Hypo.combined = RunHarmony(Hypo.combined, group.by.vars = "DatasetClean")
Hypo.combined = RunUMAP(Hypo.combined, reduction = "harmony", dims = 1:30)
saveRDS(Hypo.combined, "~/HypoMapHUMANDev.rds")



### TRANSFER HYPOMAP ANNOTATIONS
Hypo.combined@meta.data$Source = ifelse(Hypo.combined@meta.data$Dataset %in% c("Herb", "Siletti", "Zhou", "KimDev", "RomanovDev"), "NonHypoMap", "HypoMap")
DefaultAssay(Hypo.combined) = "integrated"
Hypo.Split = SplitObject(Hypo.combined, split.by = "Source")

Hypo.anchors <- FindTransferAnchors(reference = Hypo.Split[["HypoMap"]], query = Hypo.Split[["NonHypoMap"]],  dims = 1:30, features = CommonGenes$CommonGenes)
predictions <- TransferData(anchorset = Hypo.anchors, refdata = list("Dataset" = Hypo.Split[["HypoMap"]]$Dataset, "C2_named" = Hypo.Split[["HypoMap"]]$C2_named, "C7_named" = Hypo.Split[["HypoMap"]]$C7_named, "C25_named" = Hypo.Split[["HypoMap"]]$C25_named, "C66_named" = Hypo.Split[["HypoMap"]]$C66_named, "C185_named" = Hypo.Split[["HypoMap"]]$C185_named, "C286_named" = Hypo.Split[["HypoMap"]]$C286_named, "C465_named" = Hypo.Split[["HypoMap"]]$C465_named, "Region_predicted" = Hypo.Split[["HypoMap"]]$Region_predicted, "Region_summarized" = Hypo.Split[["HypoMap"]]$Region_summarized),   dims = 1:30)
saveRDS(predictions, "predictions_SplitObjRNA.rds")

for(x in names(predictions)){
  PullAnnots = predictions[[x]]
  PullAnnotations = PullAnnots %>% dplyr::select(predicted.id)
  Hypo.combined = AddMetaData(Hypo.combined, PullAnnotations, paste("HypoMapPredicted_", x, sep=""))
  PullAnnotations = PullAnnots %>% dplyr::select(prediction.score.max)
  Hypo.combined = AddMetaData(Hypo.combined, PullAnnotations, paste("HypoMapPredictionScore_", x, sep=""))
}


#Clean up dataset annotation
Hypo.combined@meta.data$DatasetClean = gsub("Affinati10x", "Affinati [Mouse Adult VMH]", gsub("Anderson10x", "Liu  [Mouse Adult VMH]", gsub("CampbellDropseq", "Campbell [Mouse Adult ARC]", gsub("ChenDropseq", "Chen [Mouse Adult]", gsub("Dowsett10xnuc", "Dowsett [Mouse Adult Hindbrain]", gsub("Flynn10x", "Mickelsen [Mouse Adult VPH]", gsub("Herb", "Herb [Human Fetal]", gsub("Kim10x", "Kim [Mouse Adult VMH]", gsub("KimDev", "Kim [Mouse Development]", gsub("LeeDropseq", "Lee [Mouse Adult]", gsub("Mickelsen10x", "Mickelsen [Mouse Adult LH]", gsub("Moffit10x", "Moffit [Mouse Adult PO]", gsub("Morris10x", "Morris [Mouse Adult SCN]", gsub("Mousebrainorg10x", "Zeisel [Mouse Adult]", gsub("RomanovDev", "Romanov [Mouse Development]", gsub("RossiDropseq", "Rossi [Mouse Adult LH", gsub("Rupp10x", "Rupp [Mouse Adult Lepr+]", gsub("Siletti", "Siletti [Human Adult]", gsub("Wen10x", "Wen [Mouse Adult SCN]", gsub("wenDropseq", "Wen [Mouse Adult SCN]", gsub("Zhou", "Zhou [Human Fetal]", Hypo.combined@meta.data$Dataset)))))))))))))))))))))


saveRDS(Hypo.combined, "~/HypoMapHUMANDev.rds")

### 


library(Seurat)
library(tidyverse)

Hypo.combined = readRDS("Hypothalamus_HumanMouseInt.rds")

HumanFetalAnnots = read.csv("HUMAN_FETAL_ASSIGNMENTS_4JUL23.csv", row.names = 1)
HumanK560 = HumanFetalAnnots %>% dplyr::select(K560)
Hypo.combined = AddMetaData(Hypo.combined, HumanK560, "HumanFetalAnnot")

Hypo.combined@meta.data$Species = ifelse(Hypo.combined@meta.data$Dataset %in% c("Herb", "Siletti", "Zhou"), "Human", "Mouse")
DefaultAssay(Hypo.combined) = "integrated"
Hypo.Split = SplitObject(Hypo.combined, split.by = "Species")

Hypo.anchors <- FindTransferAnchors(reference = Hypo.Split[["Human"]], query = Hypo.Split[["Mouse"]],  dims = 1:30)
predictions <- TransferData(anchorset = Hypo.anchors, refdata = list("HumanFetalAnnot" = Hypo.Split[["Human"]]$HumanFetalAnnot), dims = 1:30)
saveRDS(predictions, "predictions_388toHypoMap.rds")
sys.time()



MRTreeRes = readRDS("tree80_L15_treeTrim.rds")
Resolutions = as.data.frame(MRTreeRes)
Resolutions$K560 = gsub(558, 11, gsub(554, 544, gsub(551, 348,  Resolutions$K560)))
Resolutions$K494 = gsub(508, 29,  Resolutions$K494)
Resolutions$Barcodes = row.names(Resolutions)

for(x in colnames(Resolutions)){
  PullMeta = Resolutions %>% dplyr::select(x)  
  Hypo.combined = AddMetaData(Hypo.combined, PullMeta, x)  
}

Hypo.combined@meta.data$Human = ifelse(Hypo.combined@meta.data$Dataset %in% c("Herb", "Zhou"), "HumanDevelopment", "None")
Hypo.combined@meta.data$Human = ifelse(Hypo.combined@meta.data$Dataset %in% c("Siletti"), "HumanAdult", Hypo.combined@meta.data$Human)

DefaultAssay(Hypo.combined) = "integrated"
Hypo.Split = SplitObject(Hypo.combined, split.by = "Human")
Hypo.Split[["None"]] = NULL
Hypo.anchors <- FindTransferAnchors(reference = Hypo.Split[["HumanAdult"]], query = Hypo.Split[["HumanDevelopment"]],  dims = 1:30)
predictions <- TransferData(anchorset = Hypo.anchors, refdata = list("K6" = Hypo.Split[["HumanAdult"]]$K6, 
                                                                     "K16" = Hypo.Split[["HumanAdult"]]$K16, 
                                                                     "K28" = Hypo.Split[["HumanAdult"]]$K28, 
                                                                     "K40" = Hypo.Split[["HumanAdult"]]$K40, 
                                                                     "K55" = Hypo.Split[["HumanAdult"]]$K55, 
                                                                     "K116" = Hypo.Split[["HumanAdult"]]$K116, 
                                                                     "K169" = Hypo.Split[["HumanAdult"]]$K169, 
                                                                     "K199" = Hypo.Split[["HumanAdult"]]$K199, 
                                                                     "K229" = Hypo.Split[["HumanAdult"]]$K229,
                                                                     "K285" = Hypo.Split[["HumanAdult"]]$K285, 
                                                                     "K341" = Hypo.Split[["HumanAdult"]]$K341, 
                                                                     "K391" = Hypo.Split[["HumanAdult"]]$K391, 
                                                                     "K444" = Hypo.Split[["HumanAdult"]]$K444, 
                                                                     "K494" = Hypo.Split[["HumanAdult"]]$K494, 
                                                                     "K560" = Hypo.Split[["HumanAdult"]]$K560),   dims = 1:30)
saveRDS(predictions, "predictions_370toFetal.rds")