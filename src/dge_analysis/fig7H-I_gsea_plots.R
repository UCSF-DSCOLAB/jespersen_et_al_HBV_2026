library(tidyverse)
library(genekitr)

my_geneset = gmtFile_to_Term2Gene("msigdb_v2023.2.Hs_GMTs/c7.all.v2023.2.Hs.symbols.gmt")
#my_geneset = readRDS("my_geneset_c7.rds")


use_db = "c7"
fgres = read.csv("../../data/dge_result/v1_res/Fine_label/transition_ivr/Non-Naive CD4/woInteraction/outcome_gsea.csv") %>% dplyr::filter(db == use_db)
fgres_new = fgres[c(2,2,8,6,7,3,4,4)] %>% setNames(c("ID","Description","setSize","enrichmentScore","NES","pvalue","p.adjust","qvalue"))
fgres_new$rank = 1
fgres_new$leading_edge = 1
fgres_new$geneID = unname(sapply(fgres$leadingEdge, function(x) gsub(", ","/",x)))
fgres_new$Count = 1

res = read.csv("../../data/dge_result/v1_res/Fine_label/transition_ivr/Non-Naive CD4/woInteraction/outcome.csv")
my_genelist = data.frame(ID=res$X , logfc=res$logFC) %>% arrange(-logfc)

my_gse = list(gsea_df = fgres_new, geneset = my_geneset, genelist = my_genelist, exponent = 1, org = "hsapiens")

paths = list(
  "cd4_diff" = c("GSE26928_NAIVE_VS_EFF_MEMORY_CD4_TCELL_DN","GSE3982_EFF_MEMORY_VS_CENT_MEMORY_CD4_TCELL_UP","GSE3982_EFF_MEMORY_VS_CENT_MEMORY_CD4_TCELL_DN","GSE26928_NAIVE_VS_EFF_MEMORY_CD4_TCELL_UP"),
  "tfh" = c("GSE21379_TFH_VS_NON_TFH_CD4_TCELL_UP")
)
plotlist_cd4 = list()
for(p in names(paths)) {
  # Get genes to plot: select leading-edge genes that are also significantly DE.
  deg = res$X[ abs(res$logFC) > 0.25 & res$adj.P.Val < 0.1 ]
  le = unique(unlist(str_split(my_gse$gsea_df$geneID[ my_gse$gsea_df$ID %in% paths[[p]] ], "\\/")))
  #plotlist_cd4[[p]] = plotGSEA(my_gse, plot_type = "classic", show_pathway = paths[[p]], wrap_length = 20, show_gene = intersect(deg, le) )
  plotlist_cd4[[p]] = plotGSEA(my_gse, plot_type = "classic", show_pathway = paths[[p]], wrap_length = 20 )
}
pdf("../../data/dge_result/gsea_plots_nncd4_ivr.pdf", width=10, height=15); print(plotlist_cd4); dev.off()



use_db = "c7"
fgres = read.csv("../../data/dge_result/v1_res/Fine_label/transition_peakALT/Non-Naive CD8/woInteraction/outcome_gsea.csv") %>% dplyr::filter(db == use_db)
fgres_new = fgres[c(2,2,8,6,7,3,4,4)] %>% setNames(c("ID","Description","setSize","enrichmentScore","NES","pvalue","p.adjust","qvalue"))
fgres_new$rank = 1
fgres_new$leading_edge = 1
fgres_new$geneID = unname(sapply(fgres$leadingEdge, function(x) gsub(", ","/",x)))
fgres_new$Count = 1

res = read.csv("../../data/dge_result/v1_res/Fine_label/transition_peakALT/Non-Naive CD8/woInteraction/outcome.csv")
my_genelist = data.frame(ID=res$X , logfc=res$logFC) %>% arrange(-logfc)

my_gse = list(gsea_df = fgres_new, geneset = my_geneset, genelist = my_genelist, exponent = 1, org = "hsapiens")

paths = list(
  "cd8_diff_1" = c("GSE26495_NAIVE_VS_PD1LOW_CD8_TCELL_DN","GSE8678_IL7R_LOW_VS_HIGH_EFF_CD8_TCELL_DN","KAECH_NAIVE_VS_DAY8_EFF_CD8_TCELL_UP"),
  "cd8_diff_2" = c("GSE26495_PD1HIGH_VS_PD1LOW_CD8_TCELL_DN","GSE26495_PD1HIGH_VS_PD1LOW_CD8_TCELL_UP"),
  "cd8_diff_3" = c("GSE8678_IL7R_LOW_VS_HIGH_EFF_CD8_TCELL_UP"),
  "cd8_diff_4" = c("KAECH_NAIVE_VS_DAY8_EFF_CD8_TCELL_DN")
)

plotlist_cd8 = list()
for(p in names(paths)) {
  # Get genes to plot: select leading-edge genes that are also significantly DE.
  deg = res$X[ abs(res$logFC) > 0.25 & res$adj.P.Val < 0.1 ]
  le = unique(unlist(str_split(my_gse$gsea_df$geneID[ my_gse$gsea_df$ID %in% paths[[p]] ], "\\/")))
  #plotlist_cd8[[p]] = plotGSEA(my_gse, plot_type = "classic", show_pathway = paths[[p]], wrap_length = 20, show_gene = intersect(deg, le) )
  plotlist_cd8[[p]] = plotGSEA(my_gse, plot_type = "classic", show_pathway = paths[[p]], wrap_length = 20)
}
pdf("../../data/dge_result/gsea_plots_nncd8_peakALT.pdf", width=10, height=15); print(plotlist_cd8); dev.off()

