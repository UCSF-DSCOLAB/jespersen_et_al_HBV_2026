library(tidyverse)
library(ggrepel)
library(edgeR)
t = readRDS("../../data/dge_result/dge_res_v1.rds")
ys = readRDS("../../data/dge_result/dge_y_v1.rds")
df = data.frame()
pval_cutoff = 0.05
padj_cutoff = 0.1
lfc_cutoff = 0.25
mark_top_n = 100    # Top n genes to label in volcano plots
top_n_genes = 100   # Top n genes to plot in heatmap
outdir = "../../data/dge_result/v1_res/"
dir.create(outdir)
for(i in names(t)) {
  dir.create(paste(outdir,i,sep="/"), showWarnings=F)
  for(j in names(t[[i]])) {
    dir.create(paste(outdir,i,j,sep="/"), showWarnings=F)
    for(k in names(t[[i]][[j]])) {
      dir.create(paste(outdir,i,j,k,sep="/"), showWarnings=F)
      for(l in names(t[[i]][[j]][[k]])) {
        dir.create(paste(outdir,i,j,k,l,sep="/"), showWarnings=F)
        volcano_plots = list()
        heatmap_plots = list()
        for(m in names(t[[i]][[j]][[k]][[l]])) {
          #dir.create(paste(outdir,i,j,k,l,m,sep="/"), showWarnings=F)
          res = t[[i]][[j]][[k]][[l]][[m]]
          if(m == "intercept")
            next
          write.csv(res, paste0(paste(outdir,i,j,k,l,m,sep="/"),".csv") )

          # Add significance status
          res$test_status = ifelse(
            abs(res$logFC) > lfc_cutoff & res$adj.P.Val < padj_cutoff,
            "sig_padj",
            ifelse(
              abs(res$logFC) > lfc_cutoff & res$P.Value < pval_cutoff,
              "sig_pval",
              "NS"
              )
            )


          plot_title = paste0(
                j,"|",k,"|",l,"|",m,
                "\nn_pval:", sum(res$test_status == "sig_pval"),
                "\nn_padj:", sum(res$test_status == "sig_padj")
              )
          # Prepare volcano plot
          volcano_plots[[m]] <- ggplot(res, aes(logFC, -log10(P.Value))) +
            geom_point(aes(col = test_status)) +
            scale_color_manual(values = c("NS"="grey", "sig_padj"="red", "sig_pval"="orange")) +
            geom_text_repel(data = arrange(res, P.Value) %>% rownames_to_column("X") %>% dplyr::filter(test_status != "NS") %>% head(mark_top_n),
                            aes(label = X),
                            max.overlaps=20) +
            xlab("log2FC") +
            ylab("-log10(P-value)") +
            ggtitle(plot_title) +
            theme_classic()

          # Prepare heatmap
          # Prepare matrix of top_n_genes
          y = ys[[i]][[j]][[k]]
          mat = cpm(y)[ arrange(res[res$test_status != "NS",], "P.value") %>% head(top_n_genes) %>% rownames() ,]
          top_gene_lfc = arrange(res[res$test_status != "NS",], "P.value") %>% head(top_n_genes) %>% pull(logFC)
          # Calculate the log2FC of selected timepoint and Baseline
          colname_splits = str_split(colnames(mat), "_") %>% as.data.frame() %>% t()
          pats = unique(colname_splits[,c(1,3)])
          tps = unique(colname_splits[,2])
          tps_numerator = setdiff(tps, "Baseline")
          lfc_mat = data.frame()
          for(idx in 1:nrow(pats)) {
            numerator = paste0(pats[idx, 1], "_", tps_numerator, "_", pats[idx, 2])
            denominator = paste0(pats[idx, 1], "_Baseline_", pats[idx, 2])
            lfcs = log2((mat[, numerator]+1) / (mat[, denominator]+1))
            lfc_mat = rbind(
              lfc_mat,
              data.frame(lfcs) %>% t() %>% `rownames<-`(paste0(pats[idx, 1], "_", pats[idx, 2]))
              )
          }
          lfc_mat = t(lfc_mat)
          # Prepare ComplexHeatmap
          heatmap_plots[[m]] = ComplexHeatmap::Heatmap(lfc_mat,
            heatmap_legend_param = list( title = paste0("log2(",tps_numerator,"/Baseline)")),
            top_annotation = ComplexHeatmap::HeatmapAnnotation(
              outcome = str_extract(colnames(lfc_mat), "[^_]+$"),
              col = list(outcome = c("SNEGB"="black","NCTLB"="lightgrey"))
              ),
            right_annotation = ComplexHeatmap::rowAnnotation(
              limma_log2FC = top_gene_lfc,
              col = list(limma_log2FC = circlize::colorRamp2(c(-1*max(abs(top_gene_lfc)), 0, max(abs(top_gene_lfc))), c("blue", "white", "red")))
              ),
            clustering_method_columns = "ward.D2",
            clustering_method_rows = "ward.D2",
            column_title = plot_title
          )

          # Gather the DEG counts
          df<-rbind(
            df,
            data.frame(
              t(c(i,j,k,l,m,"pval",sum(res$P.Value < pval_cutoff & abs(res$logFC) > lfc_cutoff)))
            ),
            data.frame(
              t(c(i,j,k,l,m,"padj",sum(res$adj.P.Val < padj_cutoff & abs(res$logFC) > lfc_cutoff)))
            )
          )
        }
        pdf(paste0(paste(outdir,i,j,k,l,sep="/"),"/volcanos.pdf"), width=10, height=10)
        print(volcano_plots)
        dev.off()
        pdf(paste0(paste(outdir,i,j,k,l,sep="/"),"/heatmaps.pdf"), width=10, height=20)
        print(heatmap_plots)
        dev.off()
      }
    }
  }
}


df[df$X5 != "intercept",] %>% dplyr::filter(X7 > 0) %>% mutate(x=paste(X1, X3), y=paste(X2,X5,X4)) %>% ggplot(aes(x, y, size=as.numeric(X7), color=X6)) + geom_point(shape=1) + theme_classic() + theme(axis.text.x = element_text(angle=45, hjust=1))

