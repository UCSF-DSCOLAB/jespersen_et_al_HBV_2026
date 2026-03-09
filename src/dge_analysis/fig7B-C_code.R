library(tidyverse)

padj_cutoff = 0.1

ctype = "Non-Naive CD4"
term = "outcome"
data = data.frame()
degs = list()
for(tp in c("ivr","firstALT","peakALT")) {
    f = paste0("../../data/dge_result/v1_res/Fine_label/transition_",tp,"/",ctype,"/woInteraction/",term,".csv")
    res = read.csv(f)
    degs[[paste0(tp,"-up")]] = res$X[ res$logFC > 0 & res$adj.P.Val < padj_cutoff ]
    degs[[paste0(tp,"-down")]] = res$X[ res$logFC < 0 & res$adj.P.Val < padj_cutoff ]
    l = res$logFC[res$adj.P.Val < padj_cutoff]
    l_grp_tbl = case_when(
        l > 1 ~ "1",
        l > 0.75 ~ "0.75",
        l > 0.5 ~ "0.5",
        l > 0.25 ~ "0.25",
        l >= 0 ~ ">0",
        l < - 1 ~ "-1",
        l < -0.75 ~ "-0.75",
        l < -0.5 ~ "-0.5",
        l < -0.25 ~ "-0.25",
        l < 0 ~ "<0"
    ) %>% table()
    l_grp_tbl = data.frame(grp = names(l_grp_tbl), n = ifelse(grepl("-|<", names(l_grp_tbl)), -1*l_grp_tbl, l_grp_tbl), tp=tp)
    data = rbind(
        data,
        l_grp_tbl
    )
}
grp_levels = c("1","0.75","0.5","0.25",">0","<0","-0.25","-0.5","-0.75","-1")
gene_counts = group_by(data, tp) %>% summarise(c = sum(abs(n))) %>% column_to_rownames("tp")

# Add zeros for the groups that are missing in specific comparisons
uniq_grps = unique(data$grp)
for(tp in c("ivr","firstALT","peakALT")) {
  uniq_grp_tp = unique(data$grp[data$tp==tp])
  if(length(uniq_grps) != length(uniq_grp_tp)) {
    data = rbind(
      data,
      data.frame(grp = setdiff(uniq_grps, uniq_grp_tp), n = 0, tp = tp)
    )
  }
}

pdf("../../data/dge_result/deg_counts_magnitude_nonnaive_cd4.pdf", width=6, height=4)
data %>% mutate(grp = factor(grp, levels = grp_levels ), tp = factor(tp, levels = c("ivr","firstALT","peakALT"))) %>%
  ggplot(aes(tp, n, fill=grp)) + geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c(rev(pals::brewer.reds(5)), pals::brewer.blues(5)) %>% setNames(grp_levels)) +
  theme_classic() +
  scale_x_discrete(
    labels=c(
      paste0("ivr (", gene_counts["ivr","c"], ")" ),
      paste0("firstALT (", gene_counts["firstALT","c"], ")"),
      paste0("peakALT (", gene_counts["peakALT","c"], ")")))

m = make_comb_mat(degs[grep("ivr|peakALT", names(degs))])
UpSet(m,
      top_annotation = upset_top_annotation(m, add_numbers = TRUE),
    right_annotation = upset_right_annotation(m, add_numbers = TRUE))
dev.off()


