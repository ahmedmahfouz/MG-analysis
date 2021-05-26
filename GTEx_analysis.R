
library(ComplexHeatmap)
library(corrplot)

curr.dir <- "./data/"

gtex <- read.table(paste0(curr.dir, "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"), sep = "\t", comment.char = "#", header = TRUE, skip = 2)
head(gtex)
rownames(gtex) <- gtex$Name

gtex_subset <- gtex[gtex$Description %in% c("ACHE","AGRN","CHRNA1","COLQ","DOK7","LRP4","MUSK","RAPSN"),]
rownames(gtex_subset) <- gtex_subset$Description
gtex_subset <- gtex_subset[c("ACHE","AGRN","CHRNA1","COLQ","DOK7","LRP4","MUSK","RAPSN"),-c(1,2)]
head(gtex_subset)

#remove samples
colnames(gtex_subset)
gtex_subset <- gtex_subset[,-c(9,11,14,22,23)]
dim(gtex_subset)

tissue_annot <- c("Adipose","Adipose","Adrenal.Gland","Artery","Artery","Artery","Bladder",
                  "Brain","Brain","Brain","Brain","Brain","Brain","Brain",
                  "Brain","Brain","Brain","Brain","Brain","Brain","Breast",
                  "Cells","Cells","Cervix","Cervix","Colon","Colon",
                  "Esophagus","Esophagus","Esophagus","Fallopian.Tube","Heart","Heart",
                  "Kidney","Kidney","Liver","Lung","Minor.Salivary.Gland","Muscle","Nerve",
                  "Ovary","Pancreas","Pituitary","Prostate","Skin","Skin","Small.Intestine",
                  "Spleen","Stomach","Testis","Thyroid","Uterus","Vagina","Whole.Blood")
tissue_annot <- tissue_annot[-c(9,11,14,22,23)]

tiff(paste0(curr.dir, "gtex_heatmap_v5_full_dendrogram.tiff"), width = 1800, height = 2500, units = "px", res = 300)
Heatmap(t(log10(gtex_subset+1)), cluster_rows = TRUE, cluster_columns = TRUE,
        row_names_gp = gpar(fontsize = 8), row_names_side = "left",
        heatmap_legend_param = list(title = "Log10(TPM)"), column_names_side = "top",
        row_dend_side = "right", row_dend_width = unit(15, "mm"),
        clustering_method_rows = "average", clustering_method_columns = "average", column_names_gp = gpar(fontsize = 8),
        rect_gp = gpar(col = "white", lwd = 1))
dev.off()

#save eps
setEPS()
postscript(paste0(curr.dir, "gtex_heatmap_v5_full_dendrogram.eps")) 
Heatmap(t(log10(gtex_subset+1)), cluster_rows = TRUE, cluster_columns = TRUE,
        row_names_gp = gpar(fontsize = 8), row_names_side = "left",
        heatmap_legend_param = list(title = "Log10(TPM)"), column_names_side = "top",
        row_dend_side = "right", row_dend_width = unit(15, "mm"),
        clustering_method_rows = "average", clustering_method_columns = "average", column_names_gp = gpar(fontsize = 8),
        rect_gp = gpar(col = "white", lwd = 1))
dev.off()

# write.csv(gtex_subset, file = paste0(curr.dir, "mg_genes_gtex.csv"))

### process transcript data (AGRN)
# awk '$1=="ENST00000620552.4"' GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct >output_AGRN_208.txt
# head -3 GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct >output_colnames.txt

### process transcript data (ACHE)
# awk '$1=="ENST00000428317.5"' GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct >output_ACHE_207.txt

AGRN_208_expr <- read.table(paste0(curr.dir, "output_AGRN_208.txt"))
AGRN_208_expr <- as.matrix(AGRN_208_expr[,-c(1,2)])

transcript_header <- read.table(paste0(curr.dir, "output_colnames.txt"), skip = 2)
transcript_header <- as.matrix(transcript_header[,-c(1,2)])

ACHE_207_expr <- read.table(paste0(curr.dir, "output_ACHE_207.txt"))
ACHE_207_expr <- as.matrix(ACHE_207_expr[,-c(1,2)])

sample_attr <- read.delim(paste0(curr.dir, "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"), as.is = TRUE)

tissues <- unique(sample_attr$SMTSD)

tissues <- tissues[-c(2, 35, 37, 38, 46, 50)]

A <- sort(tissues)
B <- sort(colnames(gtex_subset))
C <- data.frame(A,B)

med.expr.AGRN208 <- matrix()
med.expr.ACHE207 <- matrix()
for (i in 1:length(tissues)){
  med.expr.AGRN208[i] <- median(AGRN_208_expr[,which(transcript_header %in% sample_attr$SAMPID[sample_attr$SMTSD == tissues[i] ])])
  med.expr.ACHE207[i] <- median(ACHE_207_expr[,which(transcript_header %in% sample_attr$SAMPID[sample_attr$SMTSD == tissues[i] ])])
}
names(med.expr.AGRN208) <- tissues
med.expr.AGRN208 <- as.data.frame(t(med.expr.AGRN208))

names(med.expr.ACHE207) <- tissues
med.expr.ACHE207 <- as.data.frame(t(med.expr.ACHE207))

gtex_subset <- gtex_subset[,B]
colnames(gtex_subset)

med.expr.AGRN208 <- med.expr.AGRN208[,A]
colnames(med.expr.AGRN208)
colnames(med.expr.AGRN208) <- colnames(gtex_subset)

med.expr.ACHE207 <- med.expr.ACHE207[,A]
colnames(med.expr.ACHE207)
colnames(med.expr.ACHE207) <- colnames(gtex_subset)

new.gtex.w.isoforms <- rbind(gtex_subset, med.expr.AGRN208, med.expr.ACHE207)
row.names(new.gtex.w.isoforms)[c(9,10)] <- c("AGRN-208","ACHE-207")

# plot heatmap with isoform only
new.gtex.w.isoforms <- new.gtex.w.isoforms[-c(1,2),] # remove ACHE and AGRN genes

tiff(paste0(curr.dir, "gtex_heatmap_isoforms_full_dendrogram_v4.tiff"), width = 1800, height = 2500, units = "px", res = 300)
Heatmap(t(log10(new.gtex.w.isoforms+1)), cluster_rows = TRUE, cluster_columns = TRUE,
        row_names_gp = gpar(fontsize = 8), row_names_side = "left",
        heatmap_legend_param = list(title = "Log10(TPM)"), column_names_side = "top",
        row_dend_side = "right", row_dend_width = unit(15, "mm"),
        clustering_method_rows = "average", clustering_method_columns = "average", 
        column_names_gp = gpar(fontsize = 8), rect_gp = gpar(col = "white", lwd = 1))
dev.off()

#save eps
setEPS()
postscript(paste0(curr.dir, "gtex_heatmap_isoforms_full_dendrogram_v4.eps")) 
Heatmap(t(log10(new.gtex.w.isoforms+1)), cluster_rows = TRUE, cluster_columns = TRUE,
        row_names_gp = gpar(fontsize = 8), row_names_side = "left",
        heatmap_legend_param = list(title = "Log10(TPM)"), column_names_side = "top",
        row_dend_side = "right", row_dend_width = unit(15, "mm"),
        clustering_method_rows = "average", clustering_method_columns = "average", 
        column_names_gp = gpar(fontsize = 8), rect_gp = gpar(col = "white", lwd = 1))
dev.off()

#save pdf
pdf(paste0(curr.dir, "gtex_heatmap_isoforms_full_dendrogram_v4.pdf")) 
Heatmap(t(log10(new.gtex.w.isoforms+1)), cluster_rows = TRUE, cluster_columns = TRUE,
        row_names_gp = gpar(fontsize = 8), row_names_side = "left",
        heatmap_legend_param = list(title = "Log10(TPM)"), column_names_side = "top",
        row_dend_side = "right", row_dend_width = unit(15, "mm"),
        clustering_method_rows = "average", clustering_method_columns = "average", 
        column_names_gp = gpar(fontsize = 8), rect_gp = gpar(col = "white", lwd = 1))
dev.off()

library(ggplot2)
library(RColorBrewer)
fte_theme <- function() {
  
  # Generate the colors for the chart procedurally with RColorBrewer
  
  palette <- brewer.pal("Greys", n=9)
  # color.background = palette[2]
  # color.grid.major = palette[3]
  color.background = palette[1]
  color.grid.major = palette[1]
  color.axis.text = palette[6]
  color.axis.title = palette[7]
  color.title = palette[9]
  
  # Begin construction of chart
  
  theme_bw(base_size=9) +
    
    # Set the entire chart region to a light gray color
    
    theme(panel.background=element_rect(fill=color.background, color=color.background)) +
    theme(plot.background=element_rect(fill=color.background, color=color.background)) +
    theme(panel.border=element_rect(color=color.background)) +
    
    # Format the grid
    
    theme(panel.grid.major=element_line(color=color.grid.major,size=.25)) +
    theme(panel.grid.minor=element_blank()) +
    theme(axis.ticks=element_blank()) +
    
    # Format the legend, but hide by default
    
    theme(legend.position="none") +
    theme(legend.background = element_rect(fill=color.background)) +
    theme(legend.text = element_text(size=7,color=color.axis.title)) +
    
    # Set title and axis labels, and format these and tick marks
    
    theme(plot.title=element_text(color=color.title, size=15, vjust=1.25)) +
    # theme(axis.text.x=element_text(size=7,color=color.axis.text,angle = 90, hjust = 1)) +
    theme(axis.text.x=element_text(size=10,color=color.axis.text)) +
    theme(axis.text.y=element_text(size=10,color=color.axis.text)) +
    theme(axis.title.x=element_text(size=12,color=color.axis.title, vjust=0)) +
    theme(axis.title.y=element_text(size=12,color=color.axis.title, vjust=1.25)) +
    
    # Plot margins
    
    theme(plot.margin = unit(c(0.35, 0.2, 0.3, 0.35), "cm"))
}
fte_theme_v <- function() {
  
  # Generate the colors for the chart procedurally with RColorBrewer
  
  palette <- brewer.pal("Greys", n=9)
  # color.background = palette[2]
  # color.grid.major = palette[3]
  color.background = palette[1]
  color.grid.major = palette[1]
  color.axis.text = palette[6]
  color.axis.title = palette[7]
  color.title = palette[9]
  
  # Begin construction of chart
  
  theme_bw(base_size=9) +
    
    # Set the entire chart region to a light gray color
    
    theme(panel.background=element_rect(fill=color.background, color=color.background)) +
    theme(plot.background=element_rect(fill=color.background, color=color.background)) +
    theme(panel.border=element_rect(color=color.background)) +
    
    # Format the grid
    
    theme(panel.grid.major=element_line(color=color.grid.major,size=.25)) +
    theme(panel.grid.minor=element_blank()) +
    theme(axis.ticks=element_blank()) +
    
    # Format the legend, but hide by default
    
    theme(legend.position="none") +
    theme(legend.background = element_rect(fill=color.background)) +
    theme(legend.text = element_text(size=7,color=color.axis.title)) +
    
    # Set title and axis labels, and format these and tick marks
    
    theme(plot.title=element_text(color=color.title, size=15, vjust=1.25)) +
    theme(axis.text.x=element_text(size=10,color=color.axis.text,angle = 90, hjust = 1)) +
    # theme(axis.text.x=element_text(size=10,color=color.axis.text)) +
    theme(axis.text.y=element_text(size=10,color=color.axis.text)) +
    theme(axis.title.x=element_text(size=12,color=color.axis.title, vjust=0)) +
    theme(axis.title.y=element_text(size=12,color=color.axis.title, vjust=1.25)) +
    
    # Plot margins
    
    theme(plot.margin = unit(c(0.35, 0.2, 0.3, 0.35), "cm"))
}

AGRN_208_expr <- read.table(paste0(curr.dir, "output_AGRN_208.txt"))
AGRN_208_expr <- AGRN_208_expr[,-c(1,2)]

sample_attr <- read.delim(paste0(curr.dir, "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"), as.is = TRUE)

data <- data.frame(t(AGRN_208_expr), sample_attr[which(sample_attr$SAMPID %in% transcript_header),])
head(data)

tiff(paste0(curr.dir, "gtex_boxplot_AGRN208.tiff"), width = 1800, height = 2500, units = "px", res = 300)
ggplot(data, aes(x=SMTSD, y=t.AGRN_208_expr.)) + 
  geom_boxplot(fill="#c0392b", alpha=0.75) +
  coord_flip() +
  fte_theme_v() +
  labs(x = "", y = "TPM")
dev.off()



