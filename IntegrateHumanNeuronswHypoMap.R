library(Seurat)
library(Azimuth)

hypoMap = readRDS("~/hypoMap.rds")
hypoMap = ConvertGeneNames(hypoMap, reference.names = row.names(Neuro.combined), homolog.table = "https://seurat.nygenome.org/azimuth/references/homologs.rds")
saveRDS(hypoMap, "~/hypoMapHuman.rds")

#hypoMap = readRDS("~/hypoMapHuman.rds")
Idents(hypoMap) = "C2_named"
hypoMapNeuro = subset(hypoMap, idents = c("C2-1: Neurons"))
rm(hypoMap)
EdKaZhouHypo = readRDS("EdKaZhouHypoNeurons_mt10_integrated.rds")

Hypo.anchors <- FindTransferAnchors(reference = hypoMapNeuro, query = EdKaZhouHypo,  dims = 1:30)
predictions <- TransferData(anchorset = Hypo.anchors, refdata = list("C2_named" = hypoMapNeuro$C2_named, "C7_named" = hypoMapNeuro$C7_named, "C25_named" = hypoMap$C25_named, "C66_named" = hypoMap$C66_named, "C185_named" = hypoMap$C185_named, "C286_named" = hypoMap$C286_named, "C465_named" = hypoMap$C465_named, "Region_predicted" = hypoMap$Region_predicted, "Region_summarized" = hypoMap$Region_summarized),   dims = 1:30)
EdKaZhouHypo <- AddMetaData(EdKaZhouHypo, metadata = predictions)
saveRDS(predictions, "predictions.rds")

Hypo.list = c("Mouse" = EdKaZhouHypo, "Human" = hypoMapNeuro)
rm(EdKaZhouHypo); rm(hypoMapNeuro)

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
saveRDS(Hypo.combined, "~/HypoMapHUMAN.rds")