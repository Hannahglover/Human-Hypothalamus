library(Seurat)
library(tidyverse)

#HUMAN DATA
NeoFetalNeurons = readRDS("~/NeoFetalNeurons.rds") #All Fetal + Neo Neurons
NeoFetalNeurons$Study[is.na(NeoFetalNeurons$Study)] = "Huang"
NeoFetalNeurons$StudyTime = paste(NeoFetalNeurons$Study, NeoFetalNeurons$Timepoint, sep="_")
DefaultAssay(NeoFetalNeurons) = "RNA"
HumNeoFet <- SplitObject(NeoFetalNeurons, split.by = "StudyTime")
HumNeoFet[["Huang_79D"]] = NULL
rm(NeoFetalNeurons)


EdKaZhouHypoNeurons = readRDS("~/EdKaZhouHypoNeurons_mt10_integrated_Diet.rds") #All Neurons (inc adult and fetal)
Idents(EdKaZhouHypoNeurons) = "Donor"
EdKaZhouHypoNeurons2 = subset(EdKaZhouHypoNeurons, idents = c("H18.30.002", "H19.30.001", "H19.30.002"))
DefaultAssay(EdKaZhouHypoNeurons2) = "RNA"
HumAd <- SplitObject(EdKaZhouHypoNeurons2, split.by = "Donor")

#HYPOMAP DATA
#hypoMap = readRDS("hypoMap.rds")
#hypoMap = ConvertGeneNames(hypoMap, reference.names = row.names(Neuro.combined), homolog.table = "https://seurat.nygenome.org/azimuth/references/homologs.rds")
#saveRDS(hypoMap, "~/hypoMapHuman.rds")
hypoMap = readRDS("~/hypoMapHuman.rds")
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
saveRDS(Hypo.combined, "~/HypoMapHUMANDev_30AUG23.rds")

Hypo.anchors <- FindIntegrationAnchors(object.list = Hypo.list, anchor.features = feat, reduction = "rpca")
Hypo.combined <- IntegrateData(anchorset = Hypo.anchors)

saveRDS(Hypo.combined, "~/HypoMapHUMANDev_30AUG23.rds")

Hypo.combined = read.csv("~/HypoMapHUMANDev_30AUG23.rds")
Hypo.combined = ScaleData(Hypo.combined, vars.to.regress = c("percent.mt", "nCountRNA"))
Hypo.combined = RunPCA(Hypo.combined)
Hypo.combined$DatasetClean =ifelse(Hypo.combined$Study %in% c("Huang", "Kriegstein", "Zhou"), Hypo.combined$Study, Hypo.combined$DatasetClean) 
#Hypo.combined = RunHarmony(Hypo.combined, group.by.vars = "DatasetClean")
Hypo.combined = RunUMAP(Hypo.combined, reduction = "pca", dims = 1:30)
saveRDS(Hypo.combined, "~/HypoMapHUMANDev_30AUG23_nCountScaled.rds")

#Clean up dataset annotation
Hypo.combined@meta.data$DatasetClean = gsub("Affinati10x", "Affinati [Mouse Adult VMH]", gsub("Anderson10x", "Liu  [Mouse Adult VMH]", gsub("CampbellDropseq", "Campbell [Mouse Adult ARC]", gsub("ChenDropseq", "Chen [Mouse Adult]", gsub("Dowsett10xnuc", "Dowsett [Mouse Adult Hindbrain]", gsub("Flynn10x", "Mickelsen [Mouse Adult VPH]", gsub("Herb", "Herb [Human Fetal]", gsub("Kim10x", "Kim [Mouse Adult VMH]", gsub("KimDev", "Kim [Mouse Development]", gsub("LeeDropseq", "Lee [Mouse Adult]", gsub("Mickelsen10x", "Mickelsen [Mouse Adult LH]", gsub("Moffit10x", "Moffit [Mouse Adult PO]", gsub("Morris10x", "Morris [Mouse Adult SCN]", gsub("Mousebrainorg10x", "Zeisel [Mouse Adult]", gsub("RomanovDev", "Romanov [Mouse Development]", gsub("RossiDropseq", "Rossi [Mouse Adult LH", gsub("Rupp10x", "Rupp [Mouse Adult Lepr+]", gsub("Siletti", "Siletti [Human Adult]", gsub("Wen10x", "Wen [Mouse Adult SCN]", gsub("wenDropseq", "Wen [Mouse Adult SCN]", gsub("Zhou", "Zhou [Human Fetal]", Hypo.combined@meta.data$Dataset)))))))))))))))))))))


saveRDS(Hypo.combined, "~/HypoMapHUMANDev_30AUG23.rds")
