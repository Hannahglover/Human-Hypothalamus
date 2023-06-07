##HypoCombined- WORKING on X1-16XL

library(Seurat)

KaZhouAll = readRDS("~/KaZhouAll_mt10_integrated.rds") #All Fetal
KaZhouAll@meta.data$sample2 = gsub("_.*", "", gsub("22T", "22", KaZhouAll@meta.data$sample))
KaZhouAll@meta.data$StudySample = paste(KaZhouAll@meta.data$Study, KaZhouAll@meta.data$sample2, sep="_") #17 clusters
DefaultAssay(KaZhouAll) = "RNA"
Hypo.list <- SplitObject(KaZhouAll, split.by = "StudySample")
rm(KaZhouAll)

for(n in names(Hypo.list)){
  Hypo.list[[n]] = subset(Hypo.list[[n]])  
  Hypo.list[[n]]  = DietSeurat(Hypo.list[[n]])
  MT.genes = grep(pattern="^MT-", x = row.names(Hypo.list[[n]]), value=T)
  percent.mt = Matrix::colSums(Hypo.list[[n]]@assays[["RNA"]][MT.genes, ])/Matrix::colSums(Hypo.list[[n]]@assays[["RNA"]])
  Hypo.list[[n]] = AddMetaData(Hypo.list[[n]], percent.mt, "percent.mt")
  Hypo.list[[n]] <- NormalizeData(Hypo.list[[n]])
  Hypo.list[[n]] <- FindVariableFeatures(Hypo.list[[n]], selection.method = "vst", nfeatures = 2000)
}
#saveRDS(Hypo.list, "~/Dropbox/Columbia/Hypo.listNeuro_23MAY23.rds")


EdKaZhouHypoNeurons = readRDS("~/EdKaZhouHypoNeurons_mt10_integrated.rds") #All Neurons (inc adult and fetal)
Idents(EdKaZhouHypoNeurons) = "Donor"
EdKaZhouHypoNeurons2 = subset(EdKaZhouHypoNeurons, idents = c("H18.30.002", "H19.30.001", "H19.30.002"))
DefaultAssay(EdKaZhouHypoNeurons2) = "RNA"
Hypo.list2 <- SplitObject(EdKaZhouHypoNeurons2, split.by = "Donor")

rm(EdKaZhouHypoNeurons2); rm(EdKaZhouHypoNeurons)
for(n in names(Hypo.list2)){
  Hypo.list[[n]]  = DietSeurat(Hypo.list2[[n]])
  MT.genes = grep(pattern="^MT-", x = row.names(Hypo.list[[n]]), value=T)
  percent.mt = Matrix::colSums(Hypo.list[[n]]@assays[["RNA"]][MT.genes, ])/Matrix::colSums(Hypo.list[[n]]@assays[["RNA"]])
  Hypo.list[[n]] = AddMetaData(Hypo.list[[n]], percent.mt, "percent.mt")
  Hypo.list[[n]] <- NormalizeData(Hypo.list[[n]])
  Hypo.list[[n]] <- FindVariableFeatures(Hypo.list[[n]], selection.method = "vst", nfeatures = 2000)
}
saveRDS(Hypo.list, "~/Hypo.list_1JUN23.rds")


feat <- SelectIntegrationFeatures(object.list = Hypo.list)

for(n in names(Hypo.list)){
  Hypo.list[[n]] = ScaleData(Hypo.list[[n]], features = feat, verbose = FALSE)
  Hypo.list[[n]] = RunPCA(Hypo.list[[n]], features = feat, verbose = FALSE)
}

Hypo.anchors <- FindIntegrationAnchors(object.list = Hypo.list, anchor.features = feat, reduction = "rpca")
Hypo.combined <- IntegrateData(anchorset = Hypo.anchors)

saveRDS(Hypo.combined, "~/Hypo.combined_1JUNE23.rds")