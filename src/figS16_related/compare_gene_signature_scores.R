library(tidyverse)
t = read.csv("../../data/figS16_related/pathways_fig6hi.csv")
avg_exp = readRDS("../../data/citeseq_processing/avg_exp_all_genes_fine_label_subset-NNCD4-NNCD8.rds")
meta = readRDS("GSE291286_CITEseq_cell_metadata.rds") # Available on GEO (GSE291286)
tp_meta = data.frame(pat=meta$CoLabs_patient, draw_shifted=meta$draw_shifted, ivr=meta$transition_ivr, peakalt=meta$transition_peakALT, outcome=meta$study_group_broad) %>% unique()
tp_meta$transition = unlist(apply(tp_meta, 1, function(x){ val=unique(na.omit(x[3:4])); ifelse(length(val) > 0, val, NA) }))

grp_levels = c("clearance_baseline", "persistence_baseline", "clearance_iVR", "persistence_iVR", "clearance_peakALT", "persistence_peakALT")

plotlist = list()
for(i in 1:nrow(t)) {
  # Get leading edge genes
  tp=t[i, 1]
  fine=t[i, 2]
  comp=t[i, 3]
  path=t[i, 4]
  file=t[i, 5]
  gsea = read.csv(file)
  le = unlist(str_split(gsea[gsea$pathway == path,]$leadingEdge, ", "))
  # Calculate scores for Baseline
  sample = as.data.frame(tp_meta) %>% filter(transition == "Baseline") %>% mutate(id=paste0(pat,"#",draw_shifted,"#",fine)) %>% select(id, outcome)
  sample = sample[ sample$id %in% colnames(avg_exp[[fine]]), ]
  sample$grp = paste0(ifelse(sample$outcome == "SNEGB", "clearance", "persistence"), "_baseline" )
  #sample$score = rowMeans(apply(log1p(avg_exp[[fine]][le, sample$id]), 1, scale))
  df = sample
  # Calculate scores for tp
  sample = as.data.frame(tp_meta) %>% filter(transition == tp) %>% mutate(id=paste0(pat,"#",draw_shifted,"#",fine)) %>% select(id, outcome)
  sample = sample[ sample$id %in% colnames(avg_exp[[fine]]), ]
  sample$grp = paste0(ifelse(sample$outcome == "SNEGB", "clearance", "persistence"), "_", tp )
  #sample$score = rowMeans(apply(log1p(avg_exp[[fine]][le, sample$id]), 1, scale))
  df = rbind(df, sample)
  df$score = rowMeans(apply(log1p(avg_exp[[fine]][le, df$id]), 1, scale))
  grp_lvl_tmp = grp_levels[ grp_levels %in% df$grp ]
  df$grp = factor(df$grp, levels = grp_lvl_tmp)
  comps = list()
  for(a in 1:length(grp_lvl_tmp[-1])){ for(b in (a+1):length(grp_lvl_tmp)) { comps = append(comps, list(c(grp_lvl_tmp[a], grp_lvl_tmp[b]))) } }
  comps2 = list(c("clearance_baseline","persistence_baseline"), c(paste0("clearance_",tp), paste0("persistence_",tp)))
  plotlist[[path]] = ggplot(df, aes(grp, score)) +
                      geom_boxplot(outliers=F) +
                      geom_jitter() + theme_classic() +
                      ggsignif::geom_signif(comparisons=comps2, step_increase=0.1, test="wilcox.test") +
                      theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
                      labs(y=path, title=paste0(fine))
}
pdf("../../data/figS16_related/mean_zscore_boxplot_v3_scaleAcrossBL-TP.pdf", width=8, height=11)
ggpubr::ggarrange(plotlist=plotlist, ncol=4, nrow=2)
dev.off()

