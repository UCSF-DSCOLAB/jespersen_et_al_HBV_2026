library(Seurat)
library(tidyverse)
library(limma)
library(edgeR)


xvir = readRDS("sobjs_all_both_corr_clust_rawadtNormed_v1.rds")
meta = readRDS("GSE291286_CITEseq_cell_metadata.rds") # Available on GEO (GSE291286)
xvir = AddMetaData(xvir, meta)


out_dir = "../../data/dge_result/"



# 4 samples were loaded into two pools. Here are the samples to be removed that Lia selected based on the sample quality
# https://colabs-workspace.slack.com/archives/D05NB06JBSS/p1709241019593849
remove_samples = data.frame(
  CoLabs_sample = c("XVIR1-HS41-SNEGB1","XVIR1-HS41-SNEGB6","XVIR1-HS42-SNEGB3","XVIR1-HS7-SNEGB4"),
  batch_scg = c(11,7,4,1)
  )

# Run paired analysis using limma-voom strategy for each celltype in each annotation level and for each timepoint compared to baseline.
dge_res = list()
dge_obj = list()
dge_y = list()
for(annot_clm in c("Fine_label")) {
  dge_res[[annot_clm]] = list()
  dge_obj[[annot_clm]] = list()
  dge_y[[annot_clm]] = list()
  for(tp_clm in c("transition_ivr","transition_peakALT")) {
    dge_res[[annot_clm]][[tp_clm]] = list()
    dge_obj[[annot_clm]][[tp_clm]] = list()
    dge_y[[annot_clm]][[tp_clm]] = list()
    for(ctype in na.omit(unique(xvir@meta.data[,annot_clm]))) {
      dge_res[[annot_clm]][[tp_clm]][[ctype]] = list()
      dge_res[[annot_clm]][[tp_clm]][[ctype]][["wInteraction"]] = list()
      dge_res[[annot_clm]][[tp_clm]][[ctype]][["woInteraction"]] = list()
      dge_obj[[annot_clm]][[tp_clm]][[ctype]] = list()
      print(paste(annot_clm, tp_clm, ctype))
      # Extract counts and pseudobulk for each cell type and timepoint comparison
      sub = subset(xvir, cells = Cells(xvir)[ !is.na(xvir@meta.data[,tp_clm]) & xvir@meta.data[,annot_clm] == ctype] )

      # Remove one of the duplicated samples that Lia selected based on the sample quality
      for(i in 1:nrow(remove_samples)) {
        sub = subset(sub, cells = Cells(sub)[!(sub$CoLabs_sample == remove_samples[i, "CoLabs_sample"] & sub$batch_scg == remove_samples[i, "batch_scg"])] )
      }

      counts = as.matrix(sub[['RNA']]@counts)
      pb_counts = t(rowsum(t(counts), sub$CoLabs_sample))
      sub$patient = sub$CoLabs_patient
      sub$timepoint = sub@meta.data[,tp_clm]
      sub$outcome = sub$study_group_broad

      # Keep only those patients that have both timepoints
      sample_meta = unique(sub@meta.data[,c("CoLabs_sample","patient","timepoint","outcome","batch_scg")])
      valid_patients = names(which(table(sample_meta$patient) == 2))
      sample_meta = sample_meta[ sample_meta$patient %in% valid_patients , ]
      pb_counts = pb_counts[ , colnames(pb_counts) %in% sample_meta$CoLabs_sample ]
      sample_meta = sample_meta[ match(colnames(pb_counts), sample_meta$CoLabs_sample) , ]

      # Prepare DGE object and remove lowly expressed genes, i.e. genes with CPM > 1 in > 10% of samples.
      y <- DGEList(
          counts = pb_counts,
          samples = sample_meta
      )
      y = y[ rowSums(cpm(y) > 1) > (ncol(y)*0.1) ,]
      colnames(y) = paste(y$samples$patient, y$samples$timepoint, y$samples$outcome, sep="_")
      y$samples$lib.size<-colSums(y$counts)
      y = calcNormFactors(y)
      dge_y[[annot_clm]][[tp_clm]][[ctype]] = y

      # With interaction term
      design <- model.matrix(~ batch_scg + timepoint * outcome, data=y$samples)
      v <- voom(y, design)

      dupcor <- duplicateCorrelation(v, design, block=y$samples$patient)
      print(dupcor$consensus)
      if(dupcor$consensus <= 0) {
        warning("The dupcor$consensus is not positive. This needs an attention. See Mark Segal's email.\n")
      }
      fit <- lmFit(v, design, block=y$samples$patient, correlation=dupcor$consensus)
      fit <- eBayes(fit)
      dge_res[[annot_clm]][[tp_clm]][[ctype]][["wInteraction"]][["intercept"]] <- topTable(fit, coef=1, sort.by = "P", n = Inf)
      dge_res[[annot_clm]][[tp_clm]][[ctype]][["wInteraction"]][["timepoint"]] <- topTable(fit, coef=colnames(design)[ncol(design)-2], sort.by = "P", n = Inf)
      dge_res[[annot_clm]][[tp_clm]][[ctype]][["wInteraction"]][["outcome"]] <- topTable(fit, coef=colnames(design)[ncol(design)-1], sort.by = "P", n = Inf)
      dge_res[[annot_clm]][[tp_clm]][[ctype]][["wInteraction"]][["interaction"]] <- topTable(fit, coef=colnames(design)[ncol(design)], sort.by = "P", n = Inf)
      dge_obj[[annot_clm]][[tp_clm]][[ctype]][["wInteraction"]] = list(design=design, fit=fit, consensus=dupcor$consensus)
      
      # Without interaction term
      design <- model.matrix(~ batch_scg + timepoint + outcome, data=y$samples)
      v <- voom(y, design)

      dupcor <- duplicateCorrelation(v, design, block=y$samples$patient)
      print(dupcor$consensus)
      if(dupcor$consensus <= 0) {
        warning("The dupcor$consensus is not positive. This needs an attention. See Mark Segal's email.\n")
      }
      fit <- lmFit(v, design, block=y$samples$patient, correlation=dupcor$consensus)
      fit <- eBayes(fit)
      dge_res[[annot_clm]][[tp_clm]][[ctype]][["woInteraction"]][["intercept"]] <- topTable(fit, coef=1, sort.by = "P", n = Inf)
      dge_res[[annot_clm]][[tp_clm]][[ctype]][["woInteraction"]][["timepoint"]] <- topTable(fit, coef=colnames(design)[ncol(design)-1], sort.by = "P", n = Inf)
      dge_res[[annot_clm]][[tp_clm]][[ctype]][["woInteraction"]][["outcome"]] <- topTable(fit, coef=colnames(design)[ncol(design)], sort.by = "P", n = Inf)
      dge_obj[[annot_clm]][[tp_clm]][[ctype]][["woInteraction"]] = list(design=design, fit=fit, consensus=dupcor$consensus)
    }
    saveRDS(dge_res, paste0(out_dir, "dge_res_v1.rds"))
    saveRDS(dge_obj, paste0(out_dir, "dge_obj_v1.rds"))
  }
  saveRDS(dge_y, paste0(out_dir, "dge_y_v1.rds"))
}

