library(tidyverse)
library(grid)
library(ggpubr)
xvir = readRDS("sobjs_all_both_corr_clust_rawadtNormed_v1.rds")
xvir = xvir@meta.data

calc_freq <- function(xvir = xvir@meta.data, annotation_column = "Coarse_label", return_tidy_df = FALSE) {
  SNEGBs <- xvir$CoLabs_patient[xvir$study_group_broad == "SNEGB"] %>% unique()
  NCTLBs <- xvir$CoLabs_patient[xvir$study_group_broad == "NCTLB"] %>% unique()
  freq_BL = as.data.frame(t(prop.table(table( xvir[ xvir$transition_base == "Baseline", annotation_column, drop=T], xvir$CoLabs_patient[ xvir$transition_base == "Baseline"], useNA = "ifany"),2) * 100))
  colnames(freq_BL) <- c("Subject", "Cluster", "Freq_BL")
  freq_iVR = as.data.frame(t(prop.table(table( xvir[ xvir$transition_ivr == "iVR", annotation_column, drop=T], xvir$CoLabs_patient[ xvir$transition_ivr == "iVR"], useNA = "ifany"),2) * 100))
  colnames(freq_iVR) <- c("Subject", "Cluster", "Freq_iVR")
  #freq_firstALT = as.data.frame(t(prop.table(table( xvir[ xvir$transition_firstALT == "firstALT", annotation_column, drop=T], xvir$CoLabs_patient[ xvir$transition_firstALT == "firstALT"], useNA = "ifany"),2) * 100))
  #colnames(freq_firstALT) <- c("Subject", "Cluster", "Freq_firstALT")
  freq_peakALT = as.data.frame(t(prop.table(table( xvir[ xvir$transition_peakALT == "peakALT", annotation_column, drop=T], xvir$CoLabs_patient[ xvir$transition_peakALT == "peakALT"], useNA = "ifany"),2) * 100))
  colnames(freq_peakALT) <- c("Subject", "Cluster", "Freq_peakALT")
  freq <- merge(freq_BL, freq_iVR, by = c("Subject", "Cluster"), suffixes = c("",""), all = TRUE)
  #freq <- merge(freq, freq_firstALT, by = c("Subject", "Cluster"), suffixes = c("",""), all = TRUE)
  freq <- merge(freq, freq_peakALT, by = c("Subject", "Cluster"), suffixes = c("",""), all = TRUE)
  freq$Outcome <- NA
  freq$Outcome[freq$Subject %in% SNEGBs] <- "SNEGB"
  freq$Outcome[freq$Subject %in% NCTLBs] <- "NCTLB"
  freq <- freq[,c("Subject", "Outcome", "Cluster", "Freq_BL", "Freq_iVR", "Freq_firstALT", "Freq_peakALT")]
  freq <- freq[!is.na(freq$Cluster),]
  freq <- freq[!is.na(freq$Subject),]
  freq <- freq[rowSums(is.na(freq)) == 0,]

  freq_tidy <- freq %>% select(Subject, Outcome, Cluster, Freq_BL)
  colnames(freq_tidy) <- c("Subject", "Outcome", "Cluster", "Freq")
  freq_tidy$Timepoint <- "BL"
  temp_df <- freq %>% select(Subject, Outcome, Cluster, Freq_iVR)
  colnames(temp_df) <- c("Subject", "Outcome", "Cluster", "Freq")
  temp_df$Timepoint <- "iVR"
  freq_tidy <- rbind(freq_tidy, temp_df)
  temp_df <- freq %>% select(Subject, Outcome, Cluster, Freq_firstALT)
  colnames(temp_df) <- c("Subject", "Outcome", "Cluster", "Freq")
  temp_df$Timepoint <- "firstALT"
  freq_tidy <- rbind(freq_tidy, temp_df)
  temp_df <- freq %>% select(Subject, Outcome, Cluster, Freq_peakALT)
  colnames(temp_df) <- c("Subject", "Outcome", "Cluster", "Freq")
  temp_df$Timepoint <- "peakALT"
  freq_tidy <- rbind(freq_tidy, temp_df)
  freq_tidy$Timepoint <- factor(freq_tidy$Timepoint, levels = c("BL", "iVR", "firstALT", "peakALT"))
  freq_tidy$Outcome <- factor(freq_tidy$Outcome, levels = c("SNEGB", "NCTLB"))

  if(return_tidy_df) {
    return(freq_tidy)
  }
  return(freq)
}

xvir$wsnn_res.2_wsub_mod = ifelse(is.na(xvir$Fine_label), NA, xvir$wsnn_res.2_wsub)
freq_cluster_tidy = calc_freq(xvir, "wsnn_res.2_wsub_mod", T)
freq_fine_tidy = calc_freq(xvir, "Fine_label", T)
freq_pop_tidy = calc_freq(xvir, "Population_label", T)


sig_box_plots = list()

for(g in c("fine","pop","cluster")) {
  # Calculate wilcoxon pvalues
  freq_tidy = get(paste0("freq_", g,"_tidy"))
  df = data.frame()
  for(tp in unique(freq_tidy$Timepoint)) {
    freq_tidy_tp = freq_tidy[ freq_tidy$Timepoint == tp, ]
    for(ctype in unique(freq_tidy_tp$Cluster)) {
      freq_tidy_tp_ct = freq_tidy_tp[ freq_tidy_tp$Cluster == ctype, ]
      clearers = freq_tidy_tp_ct$Freq[ freq_tidy_tp_ct$Out=="SNEGB" ]
      persisters = freq_tidy_tp_ct$Freq[ freq_tidy_tp_ct$Out=="NCTLB" ]
      df = rbind(
        df,
        data.frame(tp=tp, ctype=ctype, log2FC = log2(mean(clearers)/mean(persisters)), pval=wilcox.test(clearers, persisters)$p.value)
      )
    }
  }

  # Calculate adjusted pvalues for each timepoint separately
  df$padj = NA
  for(tp in unique(freq_tidy$Timepoint)) {
    df$padj[ df$tp == tp ] = p.adjust(df$pval[ df$tp == tp ], "BH")
  }


  signif = df %>% mutate(signif=
    case_when(
      padj < 0.0001 ~ "****",
      padj < 0.001 ~ "***",
      padj < 0.01 ~ "**",
      padj < 0.1 ~ "*",
      TRUE ~ ""
      )
    ) %>% dplyr::select(tp,ctype,signif) %>% pivot_wider(names_from="tp", values_from="signif") %>% column_to_rownames("ctype")

  df_sig = df[df$padj < 0.1,]
  for(i in 1:nrow(df_sig)) {
    tp = df_sig[i,"tp"]
    ctype = df_sig[i,"ctype"]
    padj = df_sig[i,"padj"]
    freq_tidy_tmp = dplyr::filter(freq_tidy, Timepoint == tp, Cluster == ctype)
    sig_box_plots[[paste(g,tp,ctype)]] = freq_tidy_tmp %>% 
    ggplot(aes(Outcome, Freq, color=Outcome)) + geom_boxplot(aes(fill=Outcome), outlier.shape = NA) + geom_jitter() + 
    theme_classic() + labs(y="% of total cells", title=paste0(g, "-level\n", tp, "|", ctype)) + 
    ylim(0, max(freq_tidy_tmp$Freq)*1.1) + theme(plot.title = element_text(size=10), axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
    stat_pvalue_manual(data.frame(group1="SNEGB",group2="NCTLB",padj=format.pval(padj,2), y.position=max(freq_tidy_tmp$Freq)*1.05), label = "padj") +
    scale_color_manual(values=c("NCTLB"="#0E72B3","SNEGB"="#000000")) +
    scale_fill_manual(values=c("NCTLB"="#55B4EA","SNEGB"="#6C6C6C"))
  }

  pdf(paste0("../../data/diff_freq_analysis/diff_freq_wilcox_heatmap_",g,".pdf"), width=4, height=10)
  p = df %>% dplyr::select(tp,ctype,log2FC) %>% pivot_wider(names_from="tp", values_from="log2FC") %>% column_to_rownames("ctype") %>% ComplexHeatmap::Heatmap(cluster_columns = F,
      cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%s", signif[i, j]), x, y, gp = gpar(fontsize = 10))
          }
  )
  print(p)
  dev.off()
}

pdf("../../data/diff_freq_analysis/diff_freq_wilcox_significant_boxplots.pdf", width=7, height=10)
ggarrange(plotlist=sig_box_plots)
dev.off()

