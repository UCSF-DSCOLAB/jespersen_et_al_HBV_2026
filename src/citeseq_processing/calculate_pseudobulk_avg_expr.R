library(Seurat)
library(tidyverse)
library(limma)
library(edgeR)


xvir = readRDS("sobjs_all_both_corr_clust_rawadtNormed_v1.rds")
meta = readRDS("GSE291286_CITEseq_cell_metadata.rds") # Available on GEO (GSE291286)
xvir = AddMetaData(xvir, meta)


# Keep cells other than CD4/CD8/B/Monocytes/cDCs.
xvir = subset(xvir, cells = Cells(xvir)[ !is.na(xvir$Coarse_label) ] )


celltype_levels = unique(as.character(xvir$Fine_label))
xvir$key = paste0(xvir$CoLabs_patient, "#", xvir$draw_shifted, "#", xvir$Fine_label)
avg = AverageExpression(xvir, assays='RNA', slot = "data", group.by="key")
avg_list = list()
for(ctype in celltype_levels) {
        avg_list[[ctype]] = avg$RNA %>% as.data.frame() %>% dplyr::select(ends_with(paste0("#",ctype)))
}
print("Saving the average expression")
avg_list = avg_list[c("Non-Naive CD4","Non-Naive CD8")]
saveRDS(avg_list, "../../data/citeseq_processing/avg_exp_all_genes_fine_label_subset-NNCD4-NNCD8.rds")

