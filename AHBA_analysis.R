
library(ComplexHeatmap)
library(corrplot)

curr.dir <- "./data/"

ahba <- load(paste0(curr.dir, "AHBA"))

brainExprNorm <- sapply(brainExpr, function(x) t(scale(t(x))), simplify = FALSE) # Z-score expression across samples

regions <- c("subthalamus", "claustrum", "amygdala", "hippocampal formation", "epithalamus", "cingulate gyrus",
             "thalamus", "insula", "parahippocampal gyrus", "temporal lobe", "frontal lobe", "globus pallidus",
             "parietal lobe", "mesencephalon", "basal forebrain", "striatum", "myelencephalon", "occipital lobe",
             "white matter", "hypothalamus", "pons", "cerebellum")
region_id <- sapply(regions, sample.ids, simplify = FALSE)

# samples per donor
region_id <- lapply(donorNames, function(d){
  lapply(region_id, function(r){
    e <- brainExprNorm[[d]]
    intersect(r, colnames(e))
  })
})
roi_size <- sapply(region_id, function(x)sapply(x, length))
# write.table(roi_size, file = "roi_size.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

# median expression of mg genes in selected regions
expr <- lapply(donorNames, function(d){# expression per donor
  e <- sapply(region_id[[d]], function(r){
    e <- brainExprNorm[[d]][mg_genes, ]
    if (length(r)>1){
      apply(e[, r], 1, median)  
    }
    else if (length(r)==1){
      e[, r]
    } else {
      rep(NA, nrow(e))
    }
  })
})
expr <- apply(simplify2array(expr), c(1,2), function(x) median(x, na.rm = TRUE)) # median expression across donors
rownames(expr) <- entrezId2Name(rownames(expr))
expr <- t(expr)

tiff(paste0(curr.dir, "ahba_heatmap_am_v2.tiff"), width = 1800, height = 2500, units = "px", res = 300)
Heatmap(expr, cluster_rows = TRUE, cluster_columns = TRUE,
        row_names_gp = gpar(fontsize = 8), row_names_side = "left",
        heatmap_legend_param = list(title = "Z-score(expression)"), column_names_side = "top",
        row_dend_side = "right", row_dend_width = unit(15, "mm"),
        clustering_method_rows = "average", clustering_method_columns = "average", column_names_gp = gpar(fontsize = 8),
        rect_gp = gpar(col = "white", lwd = 1))
dev.off()

#save eps
setEPS()
postscript(paste0(curr.dir, "ahba_heatmap_am_v2.eps")) 
Heatmap(expr, cluster_rows = TRUE, cluster_columns = TRUE,
        row_names_gp = gpar(fontsize = 8), row_names_side = "left",
        heatmap_legend_param = list(title = "Z-score(expression)"), column_names_side = "top",
        row_dend_side = "right", row_dend_width = unit(15, "mm"),
        clustering_method_rows = "average", clustering_method_columns = "average", column_names_gp = gpar(fontsize = 8),
        rect_gp = gpar(col = "white", lwd = 1))
dev.off()

# generate heat map with AGRN 208 isoform (probe A_24_P358462, probe id = 1025471)
row.names(probeExprNorm$donor9861)

# median expression of AGRN 208 in selected regions
iso_probe_id <- c("1025472", "1025471", "1025470", "1025469")

expr_AGRNiso <- lapply(donorNames, function(d){# expression per donor
  e <- sapply(region_id[[d]], function(r){
    e <- probeExprNorm[[d]][iso_probe_id, ]
    if (length(r)>1){
      apply(e[, r], 1, median)  
    }
    else if (length(r)==1){
      e[, r]
    } else {
      rep(NA, nrow(e))
    }
  })
})
expr_AGRNiso <- apply(simplify2array(expr_AGRNiso), c(1,2), function(x) median(x, na.rm = TRUE)) # median expression across donors
rownames(expr_AGRNiso) <- entrezId2Name(rownames(expr_AGRNiso))
expr_AGRNiso <- t(expr_AGRNiso)

expr_AGRN208 <- expr_AGRNiso[,2]
expr <- cbind(expr,expr_AGRN208)
colnames(expr)[9] <- "AGRN-208"
expr <- expr[,-2]

tiff(paste0(curr.dir, "ahba_heatmap_AGRN208_v2.tiff"), width = 1800, height = 2500, units = "px", res = 300)
Heatmap(expr, cluster_rows = TRUE, cluster_columns = TRUE,
        row_names_gp = gpar(fontsize = 8), row_names_side = "left",
        heatmap_legend_param = list(title = "Z-score(expression)"), column_names_side = "top",
        row_dend_side = "right", row_dend_width = unit(15, "mm"),
        clustering_method_rows = "average", clustering_method_columns = "average", column_names_gp = gpar(fontsize = 8),
        rect_gp = gpar(col = "white", lwd = 1))
dev.off()

#save eps
setEPS()
postscript(paste0(curr.dir, "ahba_heatmap_AGRN208_v2.eps")) 
Heatmap(expr, cluster_rows = TRUE, cluster_columns = TRUE,
        row_names_gp = gpar(fontsize = 8), row_names_side = "left",
        heatmap_legend_param = list(title = "Z-score(expression)"), column_names_side = "top",
        row_dend_side = "right", row_dend_width = unit(15, "mm"),
        clustering_method_rows = "average", clustering_method_columns = "average", column_names_gp = gpar(fontsize = 8),
        rect_gp = gpar(col = "white", lwd = 1))
dev.off()

#save pdf
pdf(paste0(curr.dir, "ahba_heatmap_AGRN208_v2.pdf")) 
Heatmap(expr, cluster_rows = TRUE, cluster_columns = TRUE,
        row_names_gp = gpar(fontsize = 8), row_names_side = "left",
        heatmap_legend_param = list(title = "Z-score(expression)"), column_names_side = "top",
        row_dend_side = "right", row_dend_width = unit(15, "mm"),
        clustering_method_rows = "average", clustering_method_columns = "average", column_names_gp = gpar(fontsize = 8),
        rect_gp = gpar(col = "white", lwd = 1))
dev.off()

#plot the expression of AGRN isoforms
colnames(expr_AGRNiso) <- probes$probe_name[probes$probe_id %in% iso_probe_id]
pdf(paste0(curr.dir, "heat_map_AGRN_probes_v2.pdf"))
Heatmap(expr_AGRNiso, cluster_rows = TRUE, cluster_columns = TRUE,
        row_names_gp = gpar(fontsize = 8), row_names_side = "left",
        heatmap_legend_param = list(title = "Z-score(expression)"), column_names_side = "top",
        row_dend_side = "right", row_dend_width = unit(15, "mm"),
        clustering_method_rows = "average", clustering_method_columns = "average", column_names_gp = gpar(fontsize = 8),
        rect_gp = gpar(col = "white", lwd = 1))
dev.off()