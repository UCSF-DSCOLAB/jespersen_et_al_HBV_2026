library(tidyverse)
library(Seurat)
library(dsb)
args = commandArgs(trailingOnly=TRUE)
s_idx = as.numeric(args[1])
proc_dir = "/krummellab/data1/immunox/XVIR1/data/single_cell_GEX/processed/"
out_dir = "/krummellab/data1/DSCoLab/XVIR1/10x/SCG2__52_py_RKP/dsb_normalization/"
pools = paste0("XVIR1-POOL-SCG", c(2:26,28:51))
s = pools[s_idx]
sobjs = list()
stat = list()
#meta_clr = read.csv("/krummellab/data1/DSCoLab/XVIR1/10x/SCG2__52_py_CL/data_files/cafs_RNA_only_UMAP_v2_metadata.csv.gz")
#sample_meta = meta_clr[,c("SCG_ID_scg","SCG_ID_scc","BEST.GUESS","CoLabs_patient","CoLabs_sample","batch_scg","batch_scc","transition_base","transition_prioritize_iVR","transition_prioritize_FirstALT","transition_prioritize_chronic","study_group_broad","study_group_fine","LIBRARY")] %>% unique()
sample_meta = read.csv("/krummellab/data1/DSCoLab/XVIR1/10x/SCG2__52_py_CL/data_files/cafs_RNA_only_UMAP_v2_sample_metadata.csv")
sample_meta$key = paste0(sample_meta$SCG_ID_scg, "#", sample_meta$BEST.GUESS)

final_junk_markers = c("CD202b--Tie2--Tek","CD338--ABCG2","mast-cell-tryptase","CD9.1","CD140b","CD140a","CD326--EpCAM","CD138--Syndecan-1","XCR1.1","TCR-G--TCR-D","CD79a--IGA","C5L2","CD336--NKp44","CD1a","CD61","TCR-V-A-24-J-A-18--iNKT cell","CD135","CD102","CD49f","CD18","CD29","CD45","CD112--Nectin-2","CD52.1","CD48.1","CD47.1","CD90--Thy1","CD269--BCMA")

pdf(paste0(out_dir,"/plots_", s, ".pdf"), width=10, height=10)
print(paste0("Processing data for ", s))
stat_local = list()

sobj_local = readRDS( file.path(proc_dir, s, paste0("/qc_and_addModalities/",s,"_seurat_object_qcAndAddModalities.rds")) )

# Using the second round of analysis, we selected 9 clusters, potential platelets and potential RBCs to remove. Removing those cell barcodes here.
br = read.csv(file.path(out_dir,"barcode_to_remove.csv"))
sobj_local = subset(sobj_local, cells=Cells(sobj_local)[ ! paste0(s, "_", Cells(sobj_local)) %in% br$x ] )

# Remove unwanted metadata variables.
sobj_local@meta.data[,c("nCount_SCT","nFeature_SCT","SCT_snn_res.0.8","seurat_clusters","pANN_0.25_0.01_6710","DF.classifications_0.25_0.01_6710","doubletFinderCalls")] = NULL

# Delete the SCT assay from doubletFinder to reduce the object sizes.
DefaultAssay(sobj_local) = 'ADT'
sobj_local[['SCT']] = NULL


meta_local = sobj_local@meta.data %>% mutate(key = paste0(orig.ident, "#", BEST.GUESS))
sample_meta_idx = match(meta_local$key, sample_meta$key)
sobj_local = AddMetaData( sobj_local,
  sample_meta[sample_meta_idx,] %>% `rownames<-`(Cells(sobj_local))
  )

# Perform dsb normalization for ADT data.
# Select 'background' droplets
raw_data <- Read10X_h5( file.path(proc_dir, s, "/cellranger/raw_feature_bc_matrix.h5") )
empty_drops <- colnames(raw_data$`Gene Expression`)[colSums(raw_data$`Gene Expression`) < 100 & colSums(raw_data$`Antibody Capture`) > 10 ]
empty_drops <- empty_drops[! empty_drops %in% Cells(sobj_local)]
ADT_background <- raw_data$`Antibody Capture`[, empty_drops]
rownames(ADT_background) = gsub("_", "-", rownames(ADT_background))
ADT_counts <- sobj_local[['ADT']]@counts

# Remove the selected markers from DSB normalization
ADT_background = ADT_background[ ! rownames(ADT_background) %in% final_junk_markers,]
ADT_counts <- ADT_counts[ ! rownames(ADT_counts) %in% final_junk_markers,]

# The feature_reference.csv used by Arjun for cellranger included extra markers which were not used for staining, but they have small counts. Here keeping only those markers that were used for staining.
used_markers = read.csv("/krummellab/data1/DSCoLab/XVIR1/10x/SCG2__52_py_RKP/dsb_normalization/feature_reference_nameMatchedToSObj_final.csv")$name
ADT_background = ADT_background[ rownames(ADT_background) %in% used_markers,]
ADT_counts <- ADT_counts[ rownames(ADT_counts) %in% used_markers,]


#print(hist(sort(apply(ADT_counts == 0, 1, sum)/ncol(ADT_counts)), breaks=100, main=s))
#ADT_counts = ADT_counts[ apply(ADT_counts, 1, mean) != 0, ]
#ADT_counts = ADT_counts[ apply(ADT_counts == 0, 1, sum)/ncol(ADT_counts) < 0.999, ]
#ADT_counts = ADT_counts[ apply(ADT_counts, 1, max) > 4, ]
#ADT_background = ADT_background[ rownames(ADT_background) %in% rownames(ADT_counts), ]

#Preserving the raw ADT data
sobj_local[['raw.ADT']] = sobj_local[['ADT']]

stat_local[[ 'empty_drops_count' ]] <- length(empty_drops)

adt_norm = DSBNormalizeProtein(
  cell_protein_matrix = ADT_counts,
  empty_drop_matrix = ADT_background,
  denoise.counts = TRUE,
  use.isotype.control = FALSE
)

adt_norm = apply(adt_norm,
                 2,
                 function(x){ ifelse(test = x < -10, yes = 0, no = x)})

sobj_local[['ADT']] = CreateAssayObject(data = adt_norm)

# Remove the cells for the samples that need to be removed. Using the sample metadata from Chris's RNA-only object for this.
sobj_local = subset(sobj_local, cells = Cells(sobj_local)[!is.na(sobj_local$key)] )

p = VlnPlot(sobj_local, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 2, pt.size = 0, group.by="CoLabs_sample")
print(p)

print(dim(sobj_local))
stat_local[[ "dim_raw" ]] <- dim(sobj_local)
stat_local[[ "droplet.types" ]] = table(sobj_local$DROPLET.TYPE.FINAL)

stat[[ s ]] = stat_local
sobjs[[ s ]] = sobj_local
dev.off()
save(sobjs,stat, file= file.path(out_dir,paste0("objs_",s,".rdata")))

