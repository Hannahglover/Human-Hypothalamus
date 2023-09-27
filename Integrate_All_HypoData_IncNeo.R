##HypoCombined- WORKING on X1-16XL

library(Seurat)
NeoFetalNeurons = readRDS("~/NeoFetalNeurons.rds") #All Fetal + Neo Neurons
NeoFetalNeurons$Study[is.na(NeoFetalNeurons$Study)] = "Huang"
NeoFetalNeurons$StudyTime = paste(NeoFetalNeurons$Study, NeoFetalNeurons$Timepoint, sep="_")
DefaultAssay(NeoFetalNeurons) = "RNA"
Hypo.list <- SplitObject(NeoFetalNeurons, split.by = "StudyTime")
rm(NeoFetalNeurons)

for(n in names(Hypo.list)){
  Hypo.list[[n]] = subset(Hypo.list[[n]])  
  Hypo.list[[n]]  = DietSeurat(Hypo.list[[n]])
  MT.genes = grep(pattern="^MT-", x = row.names(Hypo.list[[n]]), value=T)
  percent.mt = Matrix::colSums(Hypo.list[[n]]@assays[["RNA"]][MT.genes, ])/Matrix::colSums(Hypo.list[[n]]@assays[["RNA"]])
  Hypo.list[[n]] = AddMetaData(Hypo.list[[n]], percent.mt, "percent.mt")
  Hypo.list[[n]] <- NormalizeData(Hypo.list[[n]])
  Hypo.list[[n]] <- FindVariableFeatures(Hypo.list[[n]], selection.method = "vst", nfeatures = 2000)
}


EdKaZhouHypoNeurons = readRDS("~/EdKaZhouHypoNeurons_mt10_integrated_Diet.rds") #All Neurons (inc adult and fetal)
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
saveRDS(Hypo.list, "~/HypoNeoNeur.list_29AUG23.rds")
#Hypo.list = readRDS("~/HypoNeoNeur.list_29AUG23.rds")

Hypo.list[["Huang_79D"]] = NULL
feat <- SelectIntegrationFeatures(object.list = Hypo.list)

for(n in names(Hypo.list)){
  Hypo.list[[n]] = ScaleData(Hypo.list[[n]], features = feat, verbose = FALSE)
  Hypo.list[[n]] = RunPCA(Hypo.list[[n]], features = feat, verbose = FALSE)
}

Hypo.anchors <- FindIntegrationAnchors(object.list = Hypo.list, anchor.features = feat, reduction = "rpca")
Hypo.combined <- IntegrateData(anchorset = Hypo.anchors, k.weight = 50)

saveRDS(Hypo.combined, "~/HypoNeoNeur_29AUG23.rds")

Hypo.combined = ScaleData(Hypo.combined,  vars.to.regress = c("percent.mt", "nCount_RNA"))
Hypo.combined <- RunPCA(Hypo.combined, verbose = FALSE)
Hypo.combined <- RunUMAP(Hypo.combined, reduction = "pca", dims = 1:30)

saveRDS(Hypo.combined, "~/HypoNeoNeur_29AUG23.rds")
