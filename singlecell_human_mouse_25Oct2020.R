### MG genes in the Allen Mouse cortex and Hipp data 10x

library(ComplexHeatmap)
library("RColorBrewer")
library(wesanderson)


curr.dir <- "./data/"


### Human-----
mg.genes <- c("ACHE","AGRN","CHRNA1","COLQ","DOK7","LRP4","MUSK","RAPSN")

human <- read.csv(paste0(curr.dir, "trimmedmean_human_multiple_cortical_areas_SS2.csv"))

human.metadata <- read.csv(paste0(curr.dir, "metadata_human.csv"))

human.subclass <- unique(human.metadata$subclass_label)
human.subclass <- human.subclass[-which(human.subclass == "")]
human.subclass <- sort(human.subclass)
human.cluster <- unique(human.metadata$cluster_label)

human.subclass_order <- vector()
human.subclass_class <- vector()
for (i in 1:length(human.subclass)){
  human.subclass_order[i] <- unique(human.metadata$subclass_order[human.metadata$subclass_label==human.subclass[i]])
  human.subclass_class <- c(human.subclass_class, as.character(unique(human.metadata$class_label[human.metadata$subclass_label==human.subclass[i]])))
}
names(human.subclass_order) <- human.subclass
names(human.subclass_class) <- human.subclass

human.subclass <- human.subclass[order(human.subclass_order)]
human.subclass_class <- human.subclass_class[order(human.subclass_order)]

human.cluster.info <- data.frame()
for (i in 1 : length(human.cluster)){
  tmp <- human.metadata[human.metadata$cluster_label == human.cluster[i],]
  human.cluster.info[i,1:4] <- tmp[1,c("cluster_label","class_label","subclass_label","region_label")]
}

# select MG genes
data_human <- human[,mg.genes]
row.names(data_human) <- human$cluster_label

#summarize data to subclasses
human.avg_subclass <- matrix(data=NA,length(human.subclass),dim(data_human)[2])
for (i in 1:length(human.subclass)){
  human.avg_subclass[i,] <- colMeans(data_human[row.names(data_human) %in% human.cluster.info$cluster_label[human.cluster.info$subclass_label==human.subclass[i]],])
}
row.names(human.avg_subclass) <- human.subclass
colnames(human.avg_subclass) <- colnames(data_human)

row_ha = rowAnnotation(Class = human.subclass_class, col = list(Class = c("GABAergic" = "#1b9e77", "Glutamatergic" = "#d95f02", "Non-neuronal" = "#7570b3")))

tiff(paste0(curr.dir, "sc_human_subclass_v4.tiff"), width = 1800, height = 2500, units = "px", res = 300)
Heatmap(log2(human.avg_subclass+1), cluster_rows = FALSE, cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 8), row_names_side = "left",
        heatmap_legend_param = list(title = "Log2(TPM+1)"), column_names_side = "top",
        row_dend_side = "right", row_dend_width = unit(15, "mm"),
        clustering_method_rows = "average", clustering_method_columns = "average", column_names_gp = gpar(fontsize = 8),
        right_annotation = row_ha, rect_gp = gpar(col = "white", lwd = 1), width = ncol(human.avg_subclass)*unit(5, "mm"), 
        height = nrow(human.avg_subclass)*unit(5, "mm"))
dev.off()

#save eps
setEPS()
postscript(paste0(curr.dir, "sc_human_subclass_v4.eps")) 
Heatmap(log2(human.avg_subclass+1), cluster_rows = FALSE, cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 8), row_names_side = "left",
        heatmap_legend_param = list(title = "Log2(TPM+1)"), column_names_side = "top",
        row_dend_side = "right", row_dend_width = unit(15, "mm"),
        clustering_method_rows = "average", clustering_method_columns = "average", column_names_gp = gpar(fontsize = 8),
        right_annotation = row_ha, rect_gp = gpar(col = "white", lwd = 1), width = ncol(human.avg_subclass)*unit(5, "mm"), 
        height = nrow(human.avg_subclass)*unit(5, "mm"))
dev.off()

#save pdf
pdf(paste0(curr.dir, "sc_human_subclass_v4.pdf"))
Heatmap(log2(human.avg_subclass+1), cluster_rows = FALSE, cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 8), row_names_side = "left",
        heatmap_legend_param = list(title = "Log2(TPM+1)"), column_names_side = "top",
        row_dend_side = "right", row_dend_width = unit(15, "mm"),
        clustering_method_rows = "average", clustering_method_columns = "average", column_names_gp = gpar(fontsize = 8),
        right_annotation = row_ha, rect_gp = gpar(col = "white", lwd = 1), width = ncol(human.avg_subclass)*unit(5, "mm"), 
        height = nrow(human.avg_subclass)*unit(5, "mm"))

dev.off()

### Mouse-----
mg.genes_mouse <- c("Ache","Agrn","Chrna1","Colq","Dok7","Lrp4","Musk","Rapsn")

mouse <- read.csv(paste0(curr.dir, "trimmedmean_mouse_wholecortex_hippocampus_SS2.csv"))

mouse.metadata <- read.csv(paste0(curr.dir, "metadata_mouse.csv"))

#remove HIP cells
mouse.metadata <- mouse.metadata[mouse.metadata$region_label != "HIP",]

mouse.subclass <- unique(mouse.metadata$subclass_label)
mouse.subclass <- mouse.subclass[-which(mouse.subclass == "")]
mouse.subclass <- sort(mouse.subclass)
mouse.cluster <- unique(mouse.metadata$cluster_label)

mouse.subclass_order <- vector()
mouse.subclass_class <- vector()
for (i in 1:length(mouse.subclass)){
  mouse.subclass_order[i] <- unique(mouse.metadata$subclass_order[mouse.metadata$subclass_label==mouse.subclass[i]])
  mouse.subclass_class <- c(mouse.subclass_class, as.character(unique(mouse.metadata$class_label[mouse.metadata$subclass_label==mouse.subclass[i]])))
}
names(mouse.subclass_order) <- mouse.subclass
names(mouse.subclass_class) <- mouse.subclass

mouse.subclass <- mouse.subclass[order(mouse.subclass_order)]
mouse.subclass_class <- mouse.subclass_class[order(mouse.subclass_order)]

mouse.cluster.info <- data.frame()
for (i in 1 : length(mouse.cluster)){
  tmp <- mouse.metadata[mouse.metadata$cluster_label == mouse.cluster[i],]
  mouse.cluster.info[i,1:4] <- tmp[1,c("cluster_label","class_label","subclass_label","region_label")]
}

# select MG genes
data_mouse <- mouse[,mg.genes_mouse]
row.names(data_mouse) <- mouse$cluster_label

#summarize data to subclasses
mouse.avg_subclass <- matrix(data=NA,length(mouse.subclass),dim(data_mouse)[2])
for (i in 1:length(mouse.subclass)){
  mouse.avg_subclass[i,] <- colMeans(data_mouse[row.names(data_mouse) %in% mouse.cluster.info$cluster_label[mouse.cluster.info$subclass_label==mouse.subclass[i]],])
}
row.names(mouse.avg_subclass) <- mouse.subclass
colnames(mouse.avg_subclass) <- colnames(data_mouse)

# remove subclasses not present in human data
subclass_toremove <- c("CR", "Meis2", "L2/3 IT PPP", "L2 IT ENTI", "L2 IT RHP", "L2/3 ENTI", "L5 IT TPE-ENTI",
                       "L6 IT ENTI", "L5 PPP", "SUB-ProS", "CA1-ProS", "NP SUB", "CT SUB", "L6b/CT ENT", "NP PPP",
                       "V3d", "L3 RSP-ACA")

mouse.avg_subclass <- mouse.avg_subclass[-which(row.names(mouse.avg_subclass) %in% subclass_toremove),]
mouse.subclass_class <- mouse.subclass_class[-which(names(mouse.subclass_class) %in% subclass_toremove)]

row_ha = rowAnnotation(Class = mouse.subclass_class, col = list(Class = c("GABAergic" = "#1b9e77", "Glutamatergic" = "#d95f02", "Non-Neuronal" = "#7570b3")))

tiff(paste0(curr.dir, "sc_mouse_subclass_v4.tiff"), width = 1800, height = 2500, units = "px", res = 300)
Heatmap(log2(mouse.avg_subclass+1), cluster_rows = FALSE, cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 8), row_names_side = "left",
        heatmap_legend_param = list(title = "Log2(TPM+1)"), column_names_side = "top",
        row_dend_side = "right", row_dend_width = unit(15, "mm"),
        clustering_method_rows = "average", clustering_method_columns = "average", column_names_gp = gpar(fontsize = 8),
        right_annotation = row_ha, rect_gp = gpar(col = "white", lwd = 1), width = ncol(mouse.avg_subclass)*unit(5, "mm"), 
        height = nrow(mouse.avg_subclass)*unit(5, "mm"))
dev.off()

#save eps
setEPS()
postscript(paste0(curr.dir, "sc_mouse_subclass_v4.eps")) 
Heatmap(log2(mouse.avg_subclass+1), cluster_rows = FALSE, cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 8), row_names_side = "left",
        heatmap_legend_param = list(title = "Log2(TPM+1)"), column_names_side = "top",
        row_dend_side = "right", row_dend_width = unit(15, "mm"),
        clustering_method_rows = "average", clustering_method_columns = "average", column_names_gp = gpar(fontsize = 8),
        right_annotation = row_ha, rect_gp = gpar(col = "white", lwd = 1), width = ncol(mouse.avg_subclass)*unit(5, "mm"), 
        height = nrow(mouse.avg_subclass)*unit(5, "mm"))
dev.off()

#save pdf
pdf(paste0(curr.dir, "sc_mouse_subclass_v4.pdf")) 
Heatmap(log2(mouse.avg_subclass+1), cluster_rows = FALSE, cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 8), row_names_side = "left",
        heatmap_legend_param = list(title = "Log2(TPM+1)"), column_names_side = "top",
        row_dend_side = "right", row_dend_width = unit(15, "mm"),
        clustering_method_rows = "average", clustering_method_columns = "average", column_names_gp = gpar(fontsize = 8),
        right_annotation = row_ha, rect_gp = gpar(col = "white", lwd = 1), width = ncol(mouse.avg_subclass)*unit(5, "mm"), 
        height = nrow(mouse.avg_subclass)*unit(5, "mm"))
dev.off()