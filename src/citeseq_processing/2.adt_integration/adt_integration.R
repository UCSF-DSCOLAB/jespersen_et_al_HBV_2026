library(tidyverse)
library(Seurat)
library(future)

out_dir = "/krummellab/data1/DSCoLab/XVIR1/10x/SCG2__52_py_RKP/adt_integration/"

print_message <- function(message) {
  cat("[", format(Sys.time()), "]", message, "\n")
}

all_dsb_objs = paste0("/krummellab/data1/DSCoLab/XVIR1/10x/SCG2__52_py_RKP/dsb_normalization/objs_XVIR1-POOL-SCG",c(2:26,28:51),".rdata")

sobjs.list = list()

print_message("Reading dsb-normed objects")
for(f in all_dsb_objs) {
  print(f)
  load(f)
  # Delete the SCT assay from doubletFinder to reduce the object sizes.
  #DefaultAssay(sobjs[[1]]) = 'ADT'
  #sobjs[[1]][['SCT']] = NULL

  sobjs.list = append( sobjs.list, sobjs)
}

#GC
sobjs = sobjs.list
rm(sobjs.list)
gc()


## Merge data.
print_message("Merging dsb-normed objects")
sobjs_all = merge( sobjs[[1]], y = sobjs[-1], project = "XVIR1", add.cell.ids = names(sobjs) )

# GC
rm(sobjs)
gc()

sobjs_all

# Remove unnecessary columns from metadata.
sobjs_all@meta.data[, grep("DF.classi|_ADT|_SCT", colnames(sobjs_all@meta.data), value=T)] = NULL

## I see that the pool-specific effect comes only from the ADT data (pools 4 and 5 (preped by Vrinda) are very different from other 5 pools). Hence, integrating the ADT dat. Then correct the pool-effect in RNA using harmony, and use the PCA from integrated ADT data and harmony from RNA data for neighbors and UMAP/clustering.

##### Integrate ADT data across libraries (orid.ident)
DefaultAssay(sobjs_all) <- "ADT"
sobjs.list <- SplitObject(sobjs_all, split.by = "orig.ident")

# Select features that are repeatedly variable across datasets for integration
#features <- SelectIntegrationFeatures(object.list = sobjs.list)
# Above command fails since the normalized data is imported from DSB and not calculate inside Seurat. Hence just using all markers used in DSB for the following analyses.
features <- rownames(sobjs_all[["ADT"]])

# GC
rm(sobjs_all)
gc()


# For ImmunoMicrobiome, I used first library from each pool as reference to speed up the process. Also using "rpca" speeds up the process. https://github.com/satijalab/seurat/discussions/3999
# However, since for some pools there was significant heterogeneity between libraries in terms of the QC metrics, I am going to try run ADT integration without 'reference'.

## For "rpca", Data needs to be scaled and PCA needs to be calculated for each individual object.
print_message("Scaling and PCAing each library")
sobjs.list <- lapply(X = sobjs.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = TRUE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})


#anchors <- FindIntegrationAnchors(object.list = sobjs.list, anchor.features = features, reference = grep("SCG1", names(sobjs.list)), reduction = "rpca" )
print_message("Calculating anchors")
options(future.globals.maxSize = 100 * 1024^3) # Setting maxSize to 100Gb. https://github.com/satijalab/seurat/issues/1845
plan("multiprocess", workers = 100)
anchors <- FindIntegrationAnchors(object.list = sobjs.list, anchor.features = features, reduction = "rpca" )
print_message("Saving anchors")
saveRDS(anchors, file = file.path(out_dir,paste0("anchors_v1.rds")))

# GC
rm(sobjs.list)
gc()

# this command creates an 'integrated' data assay
print_message("Performing integration")
sobj_all_adt_corr <- IntegrateData(anchorset = anchors, new.assay.name = "integrated.ADT")
print_message("Saving ADT-integrated data")
saveRDS(sobj_all_adt_corr, file = file.path(out_dir,paste0("sobj_adt_corr_v1.rds")))

# GC
rm(anchors)
gc()


# Calculate S and G2M scores.
# Preparing cell-cycle genes. cc.genes is loaded with Seurat.
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
DefaultAssay(sobj_all_adt_corr) = 'RNA'
sobj_all_adt_corr = CellCycleScoring(sobj_all_adt_corr, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)


# Calculate PCA for the integrated.ADT
DefaultAssay(sobj_all_adt_corr) <- "integrated.ADT"
sobj_all_adt_corr <- ScaleData(sobj_all_adt_corr, vars.to.regress = c("nCount_ADT","nFeature_ADT","S.Score", "G2M.Score") ) %>%
                RunPCA(npcs = 30, verbose = FALSE, reduction.name = "int.apca")
saveRDS(sobj_all_adt_corr, file = file.path(out_dir,paste0("sobj_adt_corr_v1.1.rds")))

