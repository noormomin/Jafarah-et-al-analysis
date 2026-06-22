# AHA Plots 
library(ggplot2)
library(dplyr)
library(Seurat)
DimPlot(integrated, reduction = "umap", split.by ='batchid', ncol = 4, group.by  = "integrated_snn_res.0.5") # umap of each dataset in a panel
#  heat map and dot plot to show clusters----
cluster_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
pdf("heatmap_merge_integrated_0.5.pdf", width = 12, height = 20)
DoHeatmap(integrated, features = top10$gene) + NoLegend()+
  theme(axis.text.y = element_text(size = 8)) 
dev.off()

#cols = "RdYlBu") +

top3_markers <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

integrated$cluster_ordered <- factor(
  integrated$integrated_snn_res.0.5,
  levels = c("0", "2", "3", "10", "15", "16", "9", "11", "1", "13", "4", "7", "14", "5", "6", "8", "12")
)

Idents(integrated) <- "cluster_ordered"

dot_plot1 <- DotPlot(
  integrated,
  features = unique(top3_markers$gene),
  dot.scale = 6,
  cols = c("RdYlBu"),  # custom color gradient
) +
 
  scale_size(range = c(1, 6), name = "% Expressed") +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 9),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    panel.grid = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(
    title = "Marker Expression Across Annotated Cell Types",
    x = NULL,
    y = NULL
  )

# Print the plot
print(dot_plot1)
levels = 
new.cluster.ids <- c("Macrophages", "Macrophages", "Macrophages", "Macrophages", "Macrophages", "Macrophages", "Mono/Mac", "Mono/Mac", "Monocytes", "Monocytes", "Neutrophils", "Neutrophils","DC",	 "T cells/NK", "B cells", "Fibroblsts",  "other")
names(x = new.cluster.ids) <- levels(x = integrated)
integrated <- RenameIdents(object = integrated, new.cluster.ids)

dot_plot2 <- DotPlot(
  integrated,
  features = unique(top3_markers$gene),
  dot.scale = 6,
  cols = c("RdYlBu")  # custom color gradient
) +
  
  scale_size(range = c(1, 6), name = "% Expressed") +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 9),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    panel.grid = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(
    title = "Marker Expression Across Annotated Cell Types",
    x = NULL,
    y = NULL
  )

# Print the plot
print(dot_plot2)


plot_grid(dot_plot1, dot_plot2, labels = c('A', 'B'), label_size = 12)

##cluster composition in each sample----
library(ggplot2)
library(dplyr)
library(tidyr)

cluster_composition <- table(Idents(integrated), integrated$timepoint) %>%
  as.data.frame()
colnames(cluster_composition) <- c("Cluster", "Sample", "Count")

# Step 3: Convert to percent composition
cluster_composition <- cluster_composition %>%
  group_by(Sample) %>%
  mutate(Fraction = Count / sum(Count) * 100)
cluster_composition$Sample <- factor(cluster_composition$Sample, levels = c("D0", "D1", "D3/4","D5","D7","D14"))

# Step 4: Plot stacked bar
barplot_d <- ggplot(cluster_composition, aes(x = Sample, y = Fraction, fill = Cluster)) +
  geom_bar(stat = "identity", width = 0.8) +
  theme_minimal() +
  ylab("Fraction (%)") +
  xlab("") +
  ggtitle("Cluster Composition per Sample") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(barplot_d)
##DEG per cluster----
# Step 1: Filter by adjusted p-value
degs_filtered <- cluster_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  summarise(num_DEGs = n()) %>%
  ungroup()

# Step 2: Plot
barplot_e <- ggplot(degs_filtered, aes(x = reorder(cluster, num_DEGs), y = num_DEGs)) +
  geom_bar(stat = "identity", width = 0.6) +
  coord_flip() +
  theme_minimal() +
  xlab("Cluster") +
  ylab("Number of DEGs (FDR < 0.05)") +
  ggtitle("DEG Count per Cluster") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(barplot_e)


###finding human orthologs ----
filtr_upreg_neu_D0_1<-up_deg_neutrophils_D0_D1[up_deg_neutrophils_D0_D1$p_val_adj < 0.05 & abs(up_deg_neutrophils_D0_D1$avg_log2FC) > 1.5, ]#using cluster DEG
mouse_degs <- rownames(filtr_upreg_neu_D0_1) #using cluster DEG

filtr_upreg_dds <-res_1v0_clean[res_1v0_clean$log2FoldChange>0 & res_1v0_clean$padj < 0.05,]#using ddseq2
dds_up_mouse_degsrownames<-filtr_upreg_dds$gene##using ddseq2


library(biomaRt)
#depends on if ensembl si working this sday :/
#mouse <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", mirror = "uswest")
#human <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "uswest")

# Retry getLDS
#orthologs <- getLDS( attributes = c("mgi_symbol"), filters = "mgi_symbol", values = mouse_degs, mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows = TRUE)


library(orthogene)

ortholog_df <- convert_orthologs(
  gene_df = mouse_degs,
  input_species = "mouse",
  output_species = "human"
)

dds_ortholog_df<- convert_orthologs(
  gene_df = dds_up_mouse_degsrownames,
  input_species = "mouse",
  output_species = "human"
)
write.csv(dds_ortholog_df, file="dds_ortholog_results.csv")
head(ortholog_df)
# Check mapped genes only (non-NA human symbols)
ortholog_df_clean <- ortholog_df[!is.na(ortholog_df$ortholog_gene), ]

# View some results
head(ortholog_df_clean)

library(ggplot2)

# Count how many DEGs have vs. don’t have orthologs
n_total <- length(mouse_degs)
n_with_ortholog <- length(unique(ortholog_df$input_gene))
n_without_ortholog <- n_total - n_with_ortholog

# Create a data frame for plotting
df <- data.frame(
  Category = c("With Human Ortholog", "Without Human Ortholog"),
  Count = c(n_with_ortholog, n_without_ortholog)
)

# Compute fraction
df$fraction = df$Count / sum(df$Count)

# Add labels
df$label <- paste0(df$Category, "\n", round(df$fraction * 100), "%")

# Donut chart
ggplot(df, aes(x = 2, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  xlim(0.5, 2.5) +  # Creates the "hole" effect
  theme_void() +
  theme(legend.position = "none") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c("With Human Ortholog" = "tomato2", 
                               "Without Human Ortholog" = "#C0C0C0")) +
  ggtitle("Mouse DEGs with Human Orthologs")


####cell type per time point---- 

neu_cluster_composition <- table(neutrophils_integrated$timepoint, Idents(neutrophils_integrated)) %>%
  as.data.frame()
colnames(neu_cluster_composition) <- c("Timepoint", "Cluster", "Count")

# Step 3: Convert to percent composition
neu_cluster_composition <- neu_cluster_composition %>%
  group_by(Timepoint) %>%
  mutate(Fraction = Count / sum(Count) * 100)
neu_cluster_composition$Timepoint <- factor(neu_cluster_composition$Timepoint, levels = c("D0", "D1", "D3/4","D5","D7","D14"))

# Step 4: Plot stacked bar
barplot_a <- ggplot(neu_cluster_composition, aes(x = Timepoint, y = Fraction, fill = Cluster)) +
  geom_bar(stat = "identity", width = 0.8) +
  theme_minimal() +
  ylab("Fraction (%)") +
  xlab("") +
  ggtitle("Cluster Composition Over Time from MI") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(barplot_a)

VlnPlot(neutrophils_integrated, features = c("Retnlg", "S100a8"))


### volcano plot -----
# Extract the data you need for the volcano plot
library(EnhancedVolcano)
df2 <- data.frame(
  log2FC = as.numeric(as.character(res_1v0_clean$log2FoldChange)),
  pvalue = as.numeric(as.character(res_1v0_clean$pvalue)),
  padj = as.numeric(as.character(res_1v0_clean$padj)),
  gene_name = res_1v0_clean$gene
)
keyvals <- ifelse(
  df2$log2FC < -1, 'grey',
  ifelse(df2$log2FC > 1, '#561A37',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == '#561A37'] <- 'high'
names(keyvals)[keyvals == 'black'] <- 'mid'
names(keyvals)[keyvals == 'grey'] <- 'low'
lab_italics <- paste0("italic('", df2$gene_name, "')")
selectLab_italics = paste0(
  "italic('",
  c('Cxcl3','Retnlg','Ngp', 'Wfdc17','Wfdc21','Chil3','Lcn2','S100a8','Spp1','Thbs1'),
  "')")

labeled_volcano<- EnhancedVolcano(df2 ,
                lab=lab_italics,
                x = "log2FC",
                y = "padj",
                selectLab = selectLab_italics,
                title = "Neutrophils: D1 vs D0",
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                pCutoff = 0.5,
                FCcutoff = 1,
                pointSize = 4,
                labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = FALSE,
                parseLabels = TRUE,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                lengthConnectors = unit(1, 'cm'),
                arrowheads = FALSE,
                colCustom = keyvals
                )
unlabeled_volcano <-  EnhancedVolcano(df2 ,
                                      lab=NA,
                                      x = "log2FC",
                                      y = "padj",
                                      gridlines.major = FALSE,
                                      gridlines.minor = FALSE,
                                      pCutoff = 0.5,
                                      FCcutoff = 1,
                                      pointSize = 4,
                                      labSize = 6.0,
                                      labCol = 'black',
                                      labFace = 'bold',
                                      title = NULL,
                                      subtitle = NULL,
                                      caption = NULL,
                                      legendPosition = "none",
                                      colCustom = keyvals
) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
ggsave(
  "unlabeled_volcano.png",
  plot = unlabeled_volcano,
  width = 6,
  height = 6,
  dpi = 600
)

ggsave(
  "labeled_volcano.png",
  plot = labeled_volcano,
  width = 6,
  height = 6,
  dpi = 600
)