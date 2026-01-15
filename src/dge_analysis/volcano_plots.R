library(tidyverse)
library(ggrepel)
t = readRDS("../../data/dge_result/dge_res_v1.rds")
res = res = t$Fine_label$transition_ivr$`Non-Naive CD4`$woInteraction$outcome
# Add significance status
padj_cutoff = 0.1
lfc_cutoff = 0.25
res$test_status = ifelse(
  abs(res$logFC) > lfc_cutoff & res$adj.P.Val < padj_cutoff,
  "sig_padj",
  "NS"
  )
res = arrange(res, test_status)
ginfo = read.csv("10X_ref_v5.0.0/gene_info.tsv", sep="\t", header=F)

pdf("../final_figs/volcano_nonnaiveCD4_ivr_outcome.pdf", width=12, height=10)
for(mark_top_n in c(10,20,30,40,50)) {
plot_title = paste0(
  "iVR|Non-Naive CD4|woInteraction|outcome",
  "\nn_padj:", sum(res$test_status == "sig_padj"),
  "\nlabelled_top:", mark_top_n
                )

res_for_lbl =
arrange(res, adj.P.Val) %>% rownames_to_column("X") %>%
  dplyr::filter(test_status != "NS") %>% head(mark_top_n) %>%
  mutate( gtype =
    ifelse(X %in% ginfo$V3[ginfo$V1 %in% c("chrX","chrY")],
      "XY",
      ifelse(X %in% ginfo$V3[ginfo$V4 == "lncRNA"],
        "lnc",
        "pc"
        )
      )
    )

 p = ggplot(res, aes(logFC, -log10(P.Value))) +
  geom_point(aes(color = test_status), shape=16, size=3.4) +
  scale_color_manual(values = c("NS"="grey", "sig_padj"="red", "sig_pval"="orange", "XY"=pals::brewer.reds(11)[4], "lnc"="orange", "pc"=pals::brewer.reds(11)[9])) +
  geom_text_repel(data = res_for_lbl,
                  aes(label = X, color=gtype),
                  max.overlaps=20) +
  #scale_color_manual(values = c("XY"=pals::brewer.reds(11)[4], "lnc"="orange", "pc"=pals::brewer.reds(11)[9])) +
  xlab("log2FC") +
  ylab("-log10(P-value)") +
  ggtitle(plot_title) +
  theme_classic()
  print(p)
}
dev.off()







library(tidyverse)
library(ggrepel)
t = readRDS("../../data/dge_result/dge_res_v1.rds")
res = res = t$Fine_label$transition_peakALT$`Non-Naive CD8`$woInteraction$outcome
# Add significance status
padj_cutoff = 0.1
lfc_cutoff = 0.25
res$test_status = ifelse(
  abs(res$logFC) > lfc_cutoff & res$adj.P.Val < padj_cutoff,
  "sig_padj",
  "NS"
  )
res = arrange(res, test_status)
ginfo = read.csv("10X_ref_v5.0.0/gene_info.tsv", sep="\t", header=F)

pdf("../final_figs/volcano_nonnaiveCD8_peakALT_outcome.pdf", width=12, height=10)
for(mark_top_n in c(10,20,30,40,50)) {
plot_title = paste0(
  "peakALT|Non-Naive CD8|woInteraction|outcome",
  "\nn_padj:", sum(res$test_status == "sig_padj"),
  "\nlabelled_top:", mark_top_n
                )

res_for_lbl =
arrange(res, adj.P.Val) %>% rownames_to_column("X") %>%
  dplyr::filter(test_status != "NS") %>% head(mark_top_n) %>%
  mutate( gtype =
    ifelse(X %in% ginfo$V3[ginfo$V1 %in% c("chrX","chrY")],
      "XY",
      ifelse(X %in% ginfo$V3[ginfo$V4 == "lncRNA"],
        "lnc",
        "pc"
        )
      )
    )

 p = ggplot(res, aes(logFC, -log10(P.Value))) +
  geom_point(aes(color = test_status), shape=16, size=3.4) +
  scale_color_manual(values = c("NS"="grey", "sig_padj"="red", "sig_pval"="orange", "XY"=pals::brewer.reds(11)[4], "lnc"="orange", "pc"=pals::brewer.reds(11)[9])) +
  geom_text_repel(data = res_for_lbl,
                  aes(label = X, color=gtype),
                  max.overlaps=20) +
  #scale_color_manual(values = c("XY"=pals::brewer.reds(11)[4], "lnc"="orange", "pc"=pals::brewer.reds(11)[9])) +
  xlab("log2FC") +
  ylab("-log10(P-value)") +
  ggtitle(plot_title) +
  theme_classic()
  print(p)
}
dev.off()

