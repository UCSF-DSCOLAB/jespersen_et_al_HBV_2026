library(tidyverse)

make_dotplot_of_select_genes <- function(files = files, ctype = "Non-Naive CD8", grp_gene_df, preserve_input_gene_order = FALSE, custom_color_range = NULL) {
  ## Collect all tp data
  all_res = data.frame()
  for(tp in c("ivr","peakALT")) {
    f = grep(tp,
             grep( paste0("\\/",ctype,"\\/"),
                  files,
                  value = T
                  ),
             value=T
             )
    print(f)
    if(length(f) > 1) {
      warning("Found multiple files matching the criteria. Exiting the loop. Fix this and rerun.")
      break
    }
    res = read.csv(f)
    res$tp = tp
    all_res = rbind(all_res, res[,c(1,2,6,8)])
  }

  ## Check that all genes in grp_gene_df are in all_res
  missing_gene_idx = which(! grp_gene_df$g %in% all_res$X)
  if(length(missing_gene_idx) > 0) {
    stop("Following gene names do not exist in the data. If there are typos fix them and otherwise remove them from grp_gene_df and run again.\n" , paste(grp_gene_df$g[missing_gene_idx], collapse = ","))
  }

  ## Get the DGE data for input genes
  data = data.frame()
  for(tp in c("ivr","peakALT")) {
    all_res_tmp = all_res[ all_res$tp == tp, ]
    data = rbind(
      data,
      cbind(grp_gene_df[,"grp",drop=F], all_res_tmp[ match(grp_gene_df$g, all_res_tmp$X), ])
    )
  }

  ## Arrange the genes by hierarchical clustering
  # Extract lfc values
  lfc = mutate(data, key = paste(tp)) %>%
      dplyr::select(c("X", "key", "logFC")) %>%
      unique() %>% # To remove duplicated genes from calculating the gene order
      pivot_wider(names_from = key, values_from = logFC) %>%
      column_to_rownames("X")

  # Order genes based on lfc
  ch = ComplexHeatmap::Heatmap(lfc)
  gene_order = (rownames(lfc)[ ComplexHeatmap::row_order(ch) ])

  # Arrange the gene levels based on the order from the clustered heatmap
  data = data %>%
      mutate(X = factor(X, levels = gene_order),
             tp = factor(tp, levels = rev(c("ivr","peakALT")))
      )

  if(preserve_input_gene_order) {
      data = data %>%
      mutate(X = factor(X, levels = grp_gene_df$g),
             grp = factor(grp, levels = unique(grp_gene_df$grp) )
      )
  }

  if(! is.null(custom_color_range)) {
    low_val = custom_color_range[1]
    high_val = custom_color_range[2]
  } else {
    low_val = min(data$logFC)
    high_val = max(data$logFC)
  }

  ## Make plot
  p1 = data %>%
      ggplot(aes(tp, X, size=abs(logFC), color=logFC)) +
      geom_point(shape=1) + theme_classic() + theme(axis.text.x = element_text(angle=45,hjust=1), strip.text.x = element_text(size = 7)) +
      geom_point(shape=16, data = data %>% dplyr::filter(adj.P.Val < padj_cutoff)) +
      scale_colour_gradient2(low=pals::brewer.rdbu(11)[10], mid="white", high=pals::brewer.rdbu(11)[2], midpoint = 0, limit=c(low_val, high_val)) +
      facet_grid(. ~ grp, scales = "free", space = "free") + coord_flip()

  return(p1)
}

padj_cutoff = 0.1
term = "outcome"
model_type = "woInteraction"
files = list.files("../../data/dge_result/v1_res/", full.names = T, pattern = paste0(term, ".csv"), recursive = T)
files = grep(model_type, files, value = T)



## Prepare a data.frame of grp and gene mapping for Non-Naive CD8s main figure
nncd8_gene_main = data.frame()
nncd8_gene_list_file = read.csv("../../data/dge_result/CD8_gene_list_mainFig_v3.csv", check.names = F)

# prepare grp_gene_df
for(clm in colnames(nncd8_gene_list_file)) {
    nncd8_gene_main = rbind(
      nncd8_gene_main,
      data.frame("grp"=clm, "g" = setdiff(nncd8_gene_list_file[,clm],""))
    )
}

p1=make_dotplot_of_select_genes(files = files, ctype = "Non-Naive CD8", grp_gene_df = nncd8_gene_main, preserve_input_gene_order = TRUE, custom_color_range =
c(-1.6,2.8))
pdf("../../data/dge_result/Fig7G_select_gene_lfc_plots_nonnaive_cd8_main_fig_v3.pdf", width=8.5, height=2.5)
p1
dev.off()



## Prepare a data.frame of grp and gene mapping for Non-Naive CD4s main figure
nncd4_gene_main = data.frame()
nncd4_gene_list_file = read.csv("../../data/dge_result/CD4_gene_list_mainFig_v3.csv", check.names = F)

# prepare grp_gene_df
for(clm in colnames(nncd4_gene_list_file)) {
    nncd4_gene_main = rbind(
      nncd4_gene_main,
      data.frame("grp"=clm, "g" = setdiff(nncd4_gene_list_file[,clm],""))
    )
}

p1=make_dotplot_of_select_genes(files = files, ctype = "Non-Naive CD4", grp_gene_df = nncd4_gene_main, preserve_input_gene_order = TRUE, custom_color_range = c(-1.6,2.8))
pdf("../../data/dge_result/Fig7F_select_gene_lfc_plots_nonnaive_cd4_main_fig_v3.pdf", width=8.5, height=2.5)
p1
dev.off()

