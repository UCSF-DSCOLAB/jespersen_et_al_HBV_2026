library(Seurat)
library(tidyverse)
library(harmony)
library(future)


options(future.globals.maxSize = 100 * 1024^3) # Setting maxSize to 100Gb. https://github.com/satijalab/seurat/issues/1845
plan("multiprocess", workers = 10)

out_dir = "/krummellab/data1/DSCoLab/XVIR1/10x/SCG2__52_py_RKP/rna_integration_wnnUMAPClust/"
xygenes_file = "/krummellab/data1/rpatel5/data/10X_ref_v5.0.0/chrXY_genes.tsv"
print_message <- function(message) {
  cat("[", format(Sys.time()), "]", message, "\n")
}

sobjs_all_adt_corr = readRDS("/krummellab/data1/DSCoLab/XVIR1/10x/SCG2__52_py_RKP/adt_integration/sobj_adt_corr_v1.1.rds")

## Preparing cell-cycle genes. cc.genes is loaded with Seurat.
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Perform batch correction in RNA using harmony
print_message("RNA normalization, finding var. feat., scaling, PCA and harmony-based batch-correction")
DefaultAssay(sobjs_all_adt_corr) <- "RNA"
sobjs_all_both_corr <- NormalizeData(sobjs_all_adt_corr) %>%
                FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

xygenes = read.csv(xygenes_file, header=F, sep="\t")
var_feat = VariableFeatures(sobjs_all_both_corr)
print("Variable features: chrX genes")
var_feat[ var_feat %in% xygenes$V2[ xygenes$V1 == "chrX"] ]
print("Variable features: chrY genes")
var_feat[ var_feat %in% xygenes$V2[ xygenes$V1 == "chrY"] ]
print("Variable features: mitochondrial genes")
grep("^MT-", var_feat, value=T)
print("Variable features: ribosomal genes")
grep("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", var_feat, value=T)
print("Variable features: starting with HB")
grep("^HB", var_feat, value=T)
print("Variable features: starting with IG")
grep("^IG", var_feat, value=T)

# Decided to not remove the chrX/Y genes since they didn't seem to affect the sex-based bias in clustering much.
#print_message("Removing chrX/chrY genes from the list of variable features")
#print(paste0("Number of variable features before: ", length(var_feat)))
#var_feat = var_feat[ ! var_feat %in% xygenes$V2  ]
#VariableFeatures(sobjs_all_both_corr) = var_feat
#print(paste0("Number of variable features after: ", length(var_feat)))

sobjs_all_both_corr <- sobjs_all_both_corr %>%
                CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE) %>%
                ScaleData(vars.to.regress = c("percent.mt","percent.ribo","nCount_RNA","nFeature_RNA","S.Score", "G2M.Score") ) %>%
                RunPCA(features = VariableFeatures(object = .), npcs=50) %>% 
                RunHarmony("orig.ident")
rm(sobjs_all_adt_corr)
gc()
print_message("Saving object after batch-correction")
saveRDS(sobjs_all_both_corr, file = file.path(out_dir,"sobjs_all_both_corr_v1.rds"))


# FindMultiModalNeighbors
# Using prune.SNN = 1/20 which prevents too small clusters. https://github.com/satijalab/seurat/issues/4793
print_message("Finding multi-modal neighbors")
sobjs_all_both_corr <- FindMultiModalNeighbors(
  sobjs_all_both_corr, reduction.list = list("harmony", "int.apca"), 
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight", weighted.nn.name = "weighted.nn", prune.SNN = 1/20)
print_message("Saving object after finding neighbors")
saveRDS(sobjs_all_both_corr, file = file.path(out_dir,"sobjs_all_both_corr_v1.rds"))


# Clustering
print_message("Clustering")
ress = c(0.5,0.7,1,1.2,1.5,1.8,2,2.2,2.4,2.6,3,3.5,4)
for(res in ress) {
  print_message(paste0("resolution: ", res))
  sobjs_all_both_corr = FindClusters(sobjs_all_both_corr, graph.name = "wsnn", algorithm = 3, resolution = res, verbose = TRUE, n.start = 10)
  gc()
}
print_message("Saving object after clustering")
saveRDS(sobjs_all_both_corr, file = file.path(out_dir,"sobjs_all_both_corr_clust_v1.rds"))

print_message("UMAP calculations")
sobjs_all_both_corr <- RunUMAP(sobjs_all_both_corr, reduction = "harmony", reduction.name = "harmony.umap", reduction.key = "harmony.UMAP_", dims = 1:30)
sobjs_all_both_corr <- RunUMAP(sobjs_all_both_corr, reduction = "int.apca", reduction.name = "int.aumap", reduction.key = "int.aUMAP_", dims = 1:18)
sobjs_all_both_corr <- RunUMAP(sobjs_all_both_corr, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
saveRDS(sobjs_all_both_corr, file = file.path(out_dir,"sobjs_all_both_corr_clust_v1.rds"))


#DimPlot(sobjs_all_adt_corr, reduction = "wnn.umap", group.by = "orig.ident", pt.size = 0.01, split.by = "orig.ident", ncol = 6)
#DimPlot(sobjs_all_adt_corr, reduction = "harmony.umap", group.by = "orig.ident", pt.size = 0.01, split.by = "orig.ident", ncol = 6)
#DimPlot(sobjs_all_adt_corr, reduction = "int.aumap", group.by = "orig.ident", pt.size = 0.01, split.by = "orig.ident", ncol = 6)


DefaultAssay(sobjs_all_both_corr) <- "raw.ADT"
sobjs_all_both_corr <- NormalizeData(sobjs_all_both_corr, normalization.method = 'CLR', margin = 2)
saveRDS(sobjs_all_both_corr, file = file.path(out_dir,"sobjs_all_both_corr_clust_rawadtNormed_v1.rds"))


