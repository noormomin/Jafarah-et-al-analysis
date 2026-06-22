
#Introduction------------------------------------
#We are analyzing datasets from three different papers (PMID: 35023642, 32811295, and 33525909).
#These datasets are harvested from Mus musculus hearts, sorted for CD45+ cells, then underwent sngle cell RNAseq.
#All the data has been dowloaded in the preprocessed format, meaning that it underwent cellranger and some QC analysis before the following analysis.
#HJ_Momin_Lab 
#20250316





#Load required libraries----
library(Seurat)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(ggrepel)
library(patchwork)
library(scCustomize)
library(glmGamPoi)
library(SingleR)
library(celldex) 



#load preprocessed data ----

#GSE135310
GSE135310_SS<- Read10X_GEO("~/Desktop/Hilda/transcriptomics/mouse/preprocessed_data/GSE135310_RAW/GSM4005123_F8") #this is steady state (sham)
GSE135310_D7<- Read10X_GEO("~/Desktop/Hilda/transcriptomics/mouse/preprocessed_data/GSE135310_RAW/GSM4005124_F9") #this is day 7
GSE135310_D1<- Read10X_GEO("~/Desktop/Hilda/transcriptomics/mouse/preprocessed_data/GSE135310_RAW/GSM4005125_F11") #this is day 01 
GSE135310_D3<- Read10X_GEO("~/Desktop/Hilda/transcriptomics/mouse/preprocessed_data/GSE135310_RAW/GSM4005126_F12")#this is day 03
GSE135310_D5<- Read10X_GEO("~/Desktop/Hilda/transcriptomics/mouse/preprocessed_data/GSE135310_RAW/GSM4005127_G1")#this is day 05

#GSE157244
GSE157244_D4<- Read10X_GEO("~/Desktop/Hilda/transcriptomics/mouse/preprocessed_data/GSE157244_RAW/GSM4762807_HRT_WT_D4")
GSE157244_D0<- Read10X_GEO("~/Desktop/Hilda/transcriptomics/mouse/preprocessed_data/GSE157244_RAW/GSM4762818_HRT_WT_D0")

#GSE163465
GSE163465_sham<- Read10X_GEO("~/Desktop/Hilda/transcriptomics/mouse/preprocessed_data/GSE163465_RAW/GSM4985022_Sham")
GSE163465_D3<- Read10X_GEO("~/Desktop/Hilda/transcriptomics/mouse/preprocessed_data/GSE163465_RAW/GSM4985023_Day_3")
GSE163465_D7<- Read10X_GEO("~/Desktop/Hilda/transcriptomics/mouse/preprocessed_data/GSE163465_RAW/GSM4985024_Day_7")
GSE163465_D14<- Read10X_GEO("~/Desktop/Hilda/transcriptomics/mouse/preprocessed_data/GSE163465_RAW/GSM4985025_Day_14")

#Now we can create the seurat objects, as we create the seurat objects each original matrix gets removed
GSE135310_SS <- CreateSeuratObject(counts=GSE135310_SS, project="GSE135310", min.cells = 3, min.features = 200)
GSE135310_SS$timepoint <-"D0"
GSE135310_SS$batchid <-"GSE135310_SS"

GSE135310_D7<-CreateSeuratObject(counts=GSE135310_D7, project="GSE135310", min.cells = 3, min.features = 200)
GSE135310_D7$timepoint <-"D7"
GSE135310_D7$batchid <-"GSE135310_D7"
GSE135310_D1<-CreateSeuratObject(counts=GSE135310_D1, project="GSE135310", min.cells = 3, min.features = 200)
GSE135310_D1$timepoint <-"D1"
GSE135310_D1$batchid <-"GSE135310_D1"
GSE135310_D3<-CreateSeuratObject(counts=GSE135310_D3, project="GSE135310", min.cells = 3, min.features = 200)
GSE135310_D3$timepoint <-"D3"
GSE135310_D3$batchid <-"GSE135310_D3"
GSE135310_D5<-CreateSeuratObject(counts=GSE135310_D5, project="GSE135310", min.cells = 3, min.features = 200)
GSE135310_D5$timepoint <-"D5"
GSE135310_D5$batchid <-"GSE135310_D5"

GSE157244_D4<-CreateSeuratObject(counts=GSE157244_D4, project="GSE157244", min.cells = 3, min.features = 200)
GSE157244_D4$timepoint <-"D4"
GSE157244_D4$batchid <-"GSE157244_D4"
GSE157244_D0<-CreateSeuratObject(counts=GSE157244_D0, project="GSE157244", min.cells = 3, min.features = 200)
GSE157244_D0$timepoint <-"D0"
GSE157244_D0$batchid <-"GSE157244_D0"

GSE163465_sham<-CreateSeuratObject(counts=GSE163465_sham, project="GSE163465", min.cells = 3, min.features = 200)
GSE163465_sham$timepoint <-"D0"
GSE163465_sham$batchid <-"GSE163465_D0"
GSE163465_D3<-CreateSeuratObject(counts=GSE163465_D3, project="GSE163465", min.cells = 3, min.features = 200)
GSE163465_D3$timepoint <-"D3"
GSE163465_D3$batchid <-"GSE163465_D3"
GSE163465_D7<-CreateSeuratObject(counts=GSE163465_D7, project="GSE163465", min.cells = 3, min.features = 200)
GSE163465_D7$timepoint <-"D7"
GSE163465_D7$batchid <-"GSE163465_D7"
GSE163465_D14<-CreateSeuratObject(counts=GSE163465_D14, project="GSE163465", min.cells = 3, min.features = 200)
GSE163465_D14$timepoint <-"D14"
GSE163465_D14$batchid <-"GSE163465_D14"

# to preform psuedobulk rawcounts must be normalized and QC
seurat_list <- list(
  D0_A= GSE135310_SS,
  D1=GSE135310_D1,
  D3_A= GSE135310_D3,
  D5=GSE135310_D5,
  D7_A= GSE135310_D7,
  D4= GSE157244_D4, 
  D0_B=GSE157244_D0,
  D0_C=GSE163465_sham,
  D3_B= GSE163465_D3,
  D7_B=GSE163465_D7,
  D14= GSE163465_D14
)

# Calculate Mitochondrial Gene Percentage----

# For each seurat, compute the percentage of mitochondrial transcripts (using pattern "^mt-")
 
for (i in 1:length(seurat_list)) {

  seurat_list[[i]][["percent.mito"]] <- PercentageFeatureSet(seurat_list[[i]], pattern = "^mt-")
}

# Using lapply to generate FeatureScatter plots for each object in seurat_list
feature_scatter_plots <- lapply(names(seurat_list), function(name) {
  obj <- seurat_list[[name]]
  # Create FeatureScatter plot and add the object name as the title
  p2<-FeatureScatter(obj, "nCount_RNA", "percent.mito") + ggtitle(name)
  p <- FeatureScatter(obj, "nCount_RNA", "nFeature_RNA") + ggtitle(name)
  return(p+p2)
})

# Combine all plots into one using patchwork

combined_plot <- wrap_plots(feature_scatter_plots)
print(combined_plot)

#filter using QC terms
for (i in 1:length(seurat_list)) {
  seurat_list[[i]] <- subset(seurat_list[[i]], subset = nFeature_RNA > 200  & percent.mito < 5)
}


# Process each Seurat object (Normalization & Variable Feature Detection)---------

for (i in 1:length(seurat_list)) {
  seurat_list[[i]] <- NormalizeData(seurat_list[[i]])
  seurat_list[[i]] <- FindVariableFeatures(seurat_list[[i]], selection.method = "vst", nfeatures = 2000)
}

feature_scatter_plots_postQC <- lapply(names(seurat_list), function(name) {
  obj <- seurat_list[[name]]
  # Create FeatureScatter plot and add the object name as the title
  p2<-FeatureScatter(obj, "nCount_RNA", "percent.mito") + ggtitle(name)
  p <- FeatureScatter(obj, "nCount_RNA", "nFeature_RNA") + ggtitle(name)
  return(p+p2)
})

combined_plot_postQC <- wrap_plots(feature_scatter_plots_postQC)
print(combined_plot_postQC)

    experiment.aggregate <- merge(x = seurat_list[[1]], y = seurat_list[-1])

#Find the common genes across all objects in list ----

    # Split the merged object into a list by sample batch
    list <- SplitObject(experiment.aggregate, split.by = "batchid")
    #( I added this bc it kept giving me that the vectors are not the same length-HJ)
    common_genes <- Reduce(intersect, lapply(list, function(x) rownames(x)))
    #common_genes_1 <- Reduce(intersect, lapply(seurat_list, function(x) rownames(x)))
    
    # Subset each object to include only the common genes
    for (i in 1:length(list)) {
  list[[i]] <- subset(list[[i]], features = common_genes)
}

##Preform SCT transform----
options(future.globals.maxSize = 600 * 1024^2)
list <- lapply(
  list,
  function(obj) {
    SCTransform(obj, new.assay.name = "SCT", verbose = FALSE)
  }
)
rm( "GSE135310_SS",
    "GSE135310_D1",
    "GSE135310_D3",
    "GSE135310_D5",
    "GSE135310_D7",
    "GSE157244_D4", 
    "GSE157244_D0",
    "GSE163465_sham",
    "GSE163465_D3",
    "GSE163465_D7",
    "GSE163465_D14")
features<- SelectIntegrationFeatures(object.list = list, nfeatures = 3000)

list<- PrepSCTIntegration(object.list = list, anchor.features = features)


# Find integration anchors----
anchors<- FindIntegrationAnchors(object.list = list,
                                  normalization.method = "SCT",
                                  anchor.features      = features)

# Integrate
integrated<- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# At this point, `integrated` has a combined metadata. 
# You can check: head(integrated@meta.data)

#run reduction and plot----
integrated <- RunPCA(integrated)# defult npcs is 50 
# Plot the elbow plot
ElbowPlot(object = integrated, 
          ndims = 50)
#npcs<- 1:30
integrated<- RunUMAP(integrated, dims = 1:30)
DimPlot(integrated)
DimPlot(integrated, reduction = "umap", split.by ='timepoint')
DimPlot(integrated, reduction = "umap", split.by ='orig.ident')
# Explore heatmap of PCs
DimHeatmap(integrated, 
           dims = 16:30, 
           cells = 500, 
           balanced = TRUE)
print(x = integrated[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)


integrated <- FindNeighbors(integrated, dims = 1:30)
integrated <- FindClusters(integrated, resolution = 0.5)
#integrated <- FindClusters(integrated, resolution = seq(0.5, 2, 0.1))
sapply(grep("^integrated_snn_res",colnames(integrated@meta.data),value = TRUE), function(x) length(unique(integrated@meta.data[,x])))

DefaultAssay(integrated)<- "RNA"

integrated<- NormalizeData(integrated, verbose = FALSE)
all.genes <- rownames(integrated)
integrated <- ScaleData(integrated, features = all.genes, verbose = TRUE)
integrated <- JoinLayers(integrated)

# Explore different resolutions
# seurat_integrated <- FindClusters(object = integrated, resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
#seurat_integrated@meta.data %>% View()
cluster_markers <- FindAllMarkers(integrated, only.pos = TRUE)
head(cluster_markers)

cluster_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
pdf("heatmap_merge_integrated_0.5.pdf", width = 12, height = 20)
DoHeatmap(integrated, features = top10$gene) + NoLegend()+
  theme(axis.text.y = element_text(size = 8)) 
dev.off()


#SingleR to label cells----
 # For example, using ImmGenData for mouse or BlueprintEncodeData for human
ref <- ImmGenData()
head(ref$label.main)
ref2 <-  MouseRNAseqData(ensembl = FALSE)
# Extract normalized expression data (ensure it is in a suitable format)
query <- GetAssayData(integrated, assay = "SCT", layer = "data")

pred<- SingleR(test = query, ref = ref, labels = ref$label.main) 
pred_mouse<- SingleR(test = query, ref = ref2, labels = ref2$label.main) 
head(ref2)
# SingleR will compare your query data with the reference. 
#SingleR will compute similarity scores and assign a label to each cell
head(pred$labels)
integrated$SingleR.labels <- pred$labels#now we can add the predictions as labels in the integrated seurat
integrated[["SingleR.labels2"]]<-pred_mouse$labels
DimPlot(integrated)
DimPlot(integrated, group.by = "SingleR.labels", repel = TRUE)
DimPlot(integrated, group.by = c("SingleR.labels","SingleR.labels2"), repel = TRUE)
DimPlot(integrated, group.by = "SingleR.labels", repel = TRUE, split.by = 'timepoint')
DimPlot(integrated,split.by = 'timepoint')
DimPlot(integrated,split.by = 'orig.ident')
write.csv(integrated@meta.data[["integrated_snn_res.0.5"]], file = "cluster numbers_snn0.5.csv")
write.csv(integrated@meta.data[["SingleR.labels"]], file = "cluster cell type.csv")
colnames(pred)
plotScoreHeatmap(pred)

tab <- table(cluster=integrated$integrated_snn_res.0.5, label=pred$labels) 
pheatmap::pheatmap(log10(tab+10)) # using a larger pseudo-count for smoothing. 

# replace with your features
FeaturePlot(integrated, features = Neutrophil_features)
summary(pruneScores(pred))
pruneScores(pred, get.thresholds=TRUE)

library(scuttle)
mrsd.se <- celldex::ImmGenData()#ref
MI.integrated.se <- as.SingleCellExperiment(integrated, assay = "RNA")#test
mrsd.common <- intersect(rownames(MI.integrated.se), rownames(mrsd.se))
mrsd.se <- mrsd.se[mrsd.common,]
MI.integrated.se <- MI.integrated.se[mrsd.common,]
MI.integrated.se <- logNormCounts(MI.integrated.se)
MI.mrsd.pred <- SingleR(test = MI.integrated.se, ref = mrsd.se, labels = mrsd.se$label.main)
#MI.integrated[["mrsd.fine"]] <- DC.mrsd.pred$labels
MI.integrated[["mrsd.fine"]] <- MI.mrsd.pred$labels
head(MI.mrsd.pred$labels)
integrated$SingleR.labelsB <- MI.mrsd.pred$labels
colnames(MI.mrsd.pred)
plotScoreHeatmap(MI.mrsd.pred)

tab2 <- table(cluster=integrated$integrated_snn_res.0.5, label=MI.mrsd.pred$labels) 
pheatmap::pheatmap(log10(tab2+10))
DimPlot(integrated, group.by = "SingleR.labelsB", repel = TRUE)
group.by = "seurat_annotations"
summary(pruneScores(MI.mrsd.pred))
pruneScores(MI.mrsd.pred, get.thresholds=TRUE)

plotDeltaDistribution(pred)
plotDeltaDistribution(MI.mrsd.pred)

new.cluster.ids <- c("Macrophages",	"Monocytes",	"Macrophages",	"Macrophages",	"Neutrophils",	"T cells/NK",	"B cells",	"Neutrophils",	"Fibroblsts",	"Mono/Mac",	"Macrophages",	"Mono/Mac",	"noise",	"Monocytes",	"DC",	"Macrophages",	"Macrophages")
names(x = new.cluster.ids) <- levels(x = integrated)
integrated <- RenameIdents(object = integrated, new.cluster.ids)

Neutrophil_features <- c("Csf3r","Cxcr1","Fcgr3b","Fpr1","Ceacam3","Tnfrsf10c","Plxnc1", "Fcar","Abtb1","Slc25a37","Bcl6","Ncf2","Cxcr2"
                         ,"Tnfrsf1a","Cd97","Mgam","Fam65b","Thbd","Ttnfsf14","XpO6","Lilra2","Camp","Sepx1") 
#"Icam1","Itga4","Ly6g","S100a8",



# Cell Cycle Scoring ----
# A list of cell cycle markers, from Kowalczyk et al, 2015, is loaded with Seurat (updated 2019).  
# We can segregate this list into markers of G2/M phase and markers of S phase
library(gprofiler2)
#Get the list of mouse cell cycle genes
m.s.genes = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
m.g2m.genes = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name

cellcyclet<-(integrated)
DefaultAssay(cellcyclet) <- "RNA"
cellcyclet <- JoinLayers(cellcyclet)
cellcyclet <- CellCycleScoring(cellcyclet, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)
DimPlot(cellcyclet, reduction = "umap", group.by = c("Phase"),combine = FALSE, label.size = 2)

# After further analysis, you might decide to regress out the cell cycle signature if it overwhelms cell type specific signatures. We will not do that here for the sake of time. 
cellcyclet <- ScaleData(cellcyclet, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(marrow))
#start custer specific ----
#neutrophils <- subset(integrated, subset = SingleR.labels == "Neutrophils")
# Subset the neutrophil cluster
neutrophil_cells <- WhichCells(integrated, idents = "Neutrophils")
neutrophil_data <- subset(integrated, cells = neutrophil_cells)
DimPlot(neutrophil_data)

# Re-run normalization, variable feature detection, and SCTransform on just neutrophils
DefaultAssay(neutrophil_data) <- "RNA"  # Use RNA counts for SCTransform
neutrophil_data <- SCTransform(neutrophil_data, new.assay.name = "SCT", verbose = FALSE)

# (Optional) If you have multiple batches within neutrophils, split by batch
neutro.list <- SplitObject(neutrophil_data, split.by = "batchid")
neutro.list <- lapply(neutro.list, function(x) {
  SCTransform(x, new.assay.name = "SCT", verbose = FALSE)
})

# Identify integration features specific to neutrophils
features_neutro <- SelectIntegrationFeatures(object.list = neutro.list, nfeatures = 2000)

# Prepare for integration
neutro.list <- PrepSCTIntegration(object.list = neutro.list, anchor.features = features_neutro)

# Find anchors and integrate
anchors_neutro <- FindIntegrationAnchors(object.list = neutro.list, normalization.method = "SCT", anchor.features = features_neutro)
neutrophils_integrated <- IntegrateData(anchorset = anchors_neutro, normalization.method = "SCT",  k.weight = 20)
DefaultAssay(neutrophil_data) <- "integrated" 
# Run PCA and UMAP to check integration
neutrophils_integrated <- RunPCA(neutrophil_data)
ElbowPlot(neutrophils_integrated, ndims= 50)
neutrophils_integrated <- RunUMAP(neutrophils_integrated, dims = 1:10 )
neutrophils_integrated <- FindNeighbors(neutrophils_integrated, dims = 1:10)
neutrophils_integrated <- FindClusters(neutrophils_integrated,resolution = seq(0.5, 2, 0.1))
# Display the number of clusters at one resolution (example: 1.4)
sapply(grep("^integrated_snn_res", colnames(neutrophils_integrated@meta.data), value = TRUE), 
       function(x) length(unique(neutrophils_integrated@meta.data[, x])))
Idents(neutrophils_integrated)<-"integrated_snn_res.0.6"
DimPlot(neutrophils_integrated, reduction = 'umap')
DimPlot(neutrophils_integrated, split.by = "batchid")
DimPlot(neutrophils_integrated)
DimPlot(neutrophils_integrated)
DimPlot(neutrophils_integrated,reduction = "umap", group.by = "seurat_clusters")
DimPlot(neutrophils_integrated, group.by = "seurat_clusters",split.by = 'timepoint')
DimPlot(neutrophils_integrated, split.by = 'timepoint')

neu_cluster_markers <- FindAllMarkers(neutrophils_integrated, only.pos = TRUE)
write.csv(neu_cluster_markers, file = "neutro_cluster_marker.csv")

# Explore heatmap of PCs
DimHeatmap(neutrophils_integrated, 
           dims = 1:10, 
           cells = 500, 
           balanced = TRUE)
print(x = neutrophils_integrated[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)
# Plot the elbow plot
ElbowPlot(object = neutrophils_integrated, 
          ndims = 40)

DefaultAssay(neutrophils_integrated)<- "RNA"

neutrophils_integrated<- NormalizeData(neutrophils_integrated, verbose = FALSE)
all.genes.neu <- rownames(neutrophils_integrated)
neutrophils_integrated <- ScaleData(neutrophils_integrated, features = all.genes.neu, verbose = TRUE)
neutrophils_integrated <- JoinLayers(neutrophils_integrated)

Idents(neutrophils_integrated) <- neutrophils_integrated$timepoint

deg_neutrophils_D0_D1 <- FindMarkers(neutrophils_integrated, 
                               ident.1 = "D1", 
                              ident.2 =  "D0", 
                               assay = "RNA", 
                               test.use = "wilcox")


head(deg_neutrophils_D0_D1)
write.csv(deg_neutrophils_D0_D1, file = "DEG_neutrophils_D0_D1.csv")


EnhancedVolcano(deg_neutrophils_D0_D1 ,
                lab = rownames(up_deg_neutrophils_D0_D1),
                x = "avg_log2FC",
                y = "p_val_adj",
                title = "Neutrophils: D1 vs D0",
                pCutoff = 0.05,
                FCcutoff = 1)





up_deg_neutrophils_D0_D1 <-subset(deg_neutrophils_D0_D1 , avg_log2FC>0)


EnhancedVolcano(up_deg_neutrophils_D0_D1 ,
                lab = rownames(up_deg_neutrophils_D0_D1 ),
                x = "avg_log2FC",
                y = "p_val_adj",
                title = "Neutrophils: D1 vs D0",
                pCutoff = 0.05,
                FCcutoff = 1)

write.csv(up_deg_neutrophils_D0_D1, file = "Upregulated_DEG_neutrophils_D0_D1.csv")

FeaturePlot(neutrophils_integrated, feature= c("Retnlg", "Lilrb4a", "Saa3","Chil3")) #lilrb4a , saa3 chil3 NGP,FPR1
deg_neutrophils_D0_D3 <- FindMarkers(neutrophils_integrated, 
                                     ident.1 = "D3", 
                                     ident.2 =  "D0", 
                                     assay = "RNA", 
                                     test.use = "wilcox")


head(deg_neutrophils_D0_D3)

EnhancedVolcano(deg_neutrophils_D0_D3 ,
                lab = rownames(deg_neutrophils_D0_D3 ),
                x = "avg_log2FC",
                y = "p_val_adj",
                title = "Neutrophils: D3 vs D0",
                pCutoff = 0.05,
                FCcutoff = 1)
write.csv(deg_neutrophils_D0_D3, file = "DEG_neutrophils_D0_D3.csv")
neutrophils_integrated$timegroup <- ifelse(neutrophils_integrated$timepoint == "D0", "Control",
                                ifelse(neutrophils_integrated$timepoint %in% c("D1", "D3", "D4"), "Early",
                                       ifelse(neutrophils_integrated$timepoint %in% c("D5", "D7", "D14"), "Late", NA)))

# Make sure "timegroup" is a factor with the desired order (optional)
neutrophils_integrated$timegroup <- factor(neutrophils_integrated$timegroup, levels = c("Control", "Early", "Late"))

Idents(neutrophils_integrated) <- neutrophils_integrated$timegroup

deg_neutrophils_early_ctrl <- FindMarkers(neutrophils_integrated, 
                                     ident.1 = "Early", 
                                     ident.2 =  "Control", 
                                     assay = "RNA", 
                                     test.use = "wilcox")


EnhancedVolcano(deg_neutrophils_early_ctrl ,
                lab = rownames(deg_neutrophils_early_ctrl ),
                x = "avg_log2FC",
                y = "p_val_adj",
                title = "Neutrophils: Early vs Control",
                pCutoff = 0.05,
                FCcutoff = 1)
write.csv(deg_neutrophils_early_ctrl, file = "deg_neutrophils_early_ctrl.csv")
#genes_of_interest <- c("Retnlg", "Lcn2", "Anxa3","Cd177","Itgam","Itgb2","Lypd10","Lypd11","Myo1f","Pikfyve","Pram1","Ptafr","Spi1","Stx11","Stxbp3","Syk","S100a8", "S100a9", "Vamp7")  # Replace with your gene names


#Psudeobulk analysis----
# Prepare Data for DESeq2 Analysis by Extract Raw Counts and Cell Metadata
# We use the raw counts from the "RNA" assay of your neutrophils_integrated object
# (assuming it contains the untransformed counts)

raw_counts <- as.matrix(GetAssayData(neutrophils_integrated, assay = "RNA", layer = "counts"))
raw_counts <-raw_counts+1 

# Extract condition metadata
metadata <- neutrophils_integrated@meta.data
metadata <- data.frame(condition = metadata$timepoint)
rownames(metadata) <- colnames(raw_counts)

# Ensure conditions are factors with explicit levels
metadata$condition <- factor(metadata$condition, levels = c("D0","D1", "D3","D4","D7","D14"))

# Align raw_counts and metadata if necessary
common_barcodes <- intersect(rownames(metadata), colnames(raw_counts))
print(common_barcodes)

metadata <- metadata[common_barcodes, , drop = FALSE]
raw_counts <- raw_counts[, common_barcodes]

metadata <- metadata[!is.na(metadata$condition), , drop = FALSE]
raw_counts <- raw_counts[, rownames(metadata)]
# Optionally, verify dimensions:
dim(raw_counts)    # Should be: (number of genes) x (number of cells)
dim(metadata)      # Should be: (number of cells) x (number of metadata columns)

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = raw_counts, 
  colData = metadata, 
  design = ~ condition  # Specify the design formula
)

# Pre-filter genes with low counts
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq2 differential expression analysis
dds <- DESeq(dds)

# Extract results for MI_1D vs MI_CTR comparison
res_1v0 <- results(dds, contrast = c("condition", "D1", "D0"))

res_1v0$significance <- "Not Significant"
res_1v0$significance[res_1v0$padj < 0.05 & res_1v0$log2FoldChange > 0] <- "Upregulated"
res_1v0$significance[res_1v0$padj < 0.05 & res_1v0$log2FoldChange < 0] <- "Downregulated"

# Preview top genes
head(res_1v0[order(res_1v0$padj), ])

# Convert results to data frame and add gene names
res_1v0 <- as.data.frame(res_1v0)
res_1v0$gene <- rownames(res_1v0)
rownames(res_1v0) <- NULL
write.csv(res_1v0,file="ddseq2_res_D1vsD0.csv")
# Clean data: remove rows with NA in padj or log2FoldChange, and adjust extreme values
res_1v0_clean <- res_1v0[!is.na(res_1v0$padj) & !is.na(res_1v0$log2FoldChange), ]
dim(res_1v0_clean)
res_1v0_clean$padj <- pmax(res_1v0_clean$padj, 1e-300)
res_1v0_clean$log2FoldChange <- pmin(pmax(res_1v0_clean$log2FoldChange, -10), 10)
write.csv(res_1v0_clean,file="ddseq2_res_filtered_D1vsD0.csv")

library(ggrepel)
library(ggplot2)
library(dplyr)

# Filter top upregulated genes (example: top 21 genes)
top_20_upregulated <- res_1v0_clean %>%
  filter(significance == "Upregulated") %>%
  arrange(desc(log2FoldChange)) %>%
  slice(1:21)

ggplot(res_1v0, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_text_repel(data = top_20_upregulated, aes(label = gene), size = 3.5, color = "black", max.overlaps = 21) +
  labs(title = "Volcano Plot: Top 20 Upregulated Genes",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value",
       color = "Gene Regulation") +
  theme_minimal() +
  theme(text = element_text(size = 14), legend.position = "top")

neutrophils_integrated <- RunUMAP(neutrophils_integrated, dims = 1:30 )
DimPlot(neutrophils_integrated, group.by = "seurat_clusters")

FeaturePlot(neutrophils_integrated1,features = c("Retnlg","Wfdc21","Lcn2","S100a8","Spp1","Thbs1"))
Idents(neutrophils_integrated1)<-neutrophils_integrated$seurat_clusters
VlnPlot(neutrophils_integrated1, features = c("Retnlg","Wfdc21","Lcn2","S100a8","Spp1","Thbs1"), pt.size = 0, assay = "RNA") 
+ NoLegend() + geom_boxplot(width=0.1, fill="white")
p2<-VlnPlot(neutrophils_integrated, features = c("Retnlg"), pt.size = 0, assay = "RNA") 
+ NoLegend() + geom_boxplot(width=0.1, fill="white")

#cols = c("0" = "#e68613", "1" = "#f8766d", "2" = "#7cae00", "3" = "#cd9600", "4" = "#aba300", "5" = "#0cb702", "6" = "#00be67", "7" = "#00c19a", "8" = "#c77cff", "9" = "#00bfc4", "10" = "#00b8e7"
# "11" = "#00a9ff", "12" = "#8494ff", "13" = "#ed68ed", "14" = "#ff61cc", "15" = "#ff68a1"






# EL----
# Subset the neutrophil cluster
subset_cells <- WhichCells(integrated, idents = "Macrophages")
subset_data <- subset(integrated, cells = subset_cells)
DimPlot(subset_data)

# Check batch distribution
DimPlot(subset_data, group.by = "batchid")

# Re-run normalization, variable feature detection, and SCTransform on just neutrophils
DefaultAssay(subset_data) <- "RNA"  # Use RNA counts for SCTransform
library(future)
options(future.globals.maxSize = 2 * 1024^3)  # 2 GB
subset_data <- SCTransform(subset_data, new.assay.name = "SCT", verbose = FALSE)

# (Optional) If you have multiple batches within neutrophils, split by batch
subset.list <- SplitObject(subset_data, split.by = "batchid")
subset.list <- lapply(subset.list, function(x) {
  SCTransform(x, new.assay.name = "SCT", verbose = FALSE)
})

# Identify integration features specific to neutrophils
features_subset <- SelectIntegrationFeatures(object.list = subset.list, nfeatures = 2000)

# Prepare for integration
subset.list <- PrepSCTIntegration(object.list = subset.list, anchor.features = features_subset)

# Find anchors and integrate
anchors_subset <- FindIntegrationAnchors(object.list = subset.list, normalization.method = "SCT", anchor.features = features_subset, dims=1:20)
cell_subset_integrated <- IntegrateData(anchorset = anchors_subset, normalization.method = "SCT",  k.weight = 20)
DefaultAssay(subset_data) <- "integrated" 

DefaultAssay(cell_subset_integrated)<- "RNA"

cell_subset_integrated<- NormalizeData(cell_subset_integrated, verbose = FALSE)
all.genes.subset <- rownames(cell_subset_integrated)
cell_subset_integrated <- ScaleData(cell_subset_integrated, features = all.genes.subset, verbose = TRUE)
cell_subset_integrated <- JoinLayers(cell_subset_integrated)

Idents(cell_subset_integrated) <- cell_subset_integrated$timepoint

deg_cell_subset_D0_D34 <- FindMarkers(cell_subset_integrated, 
                                     ident.1 = "D3/4", 
                                     ident.2 =  "D0", 
                                     assay = "RNA", 
                                     test.use = "wilcox")

head(deg_cell_subset_D0_D34)

# Make sure rownames are gene names
all_genes <- rownames(deg_cell_subset_D0_D34)

# Select top 5 upregulated genes based on log2FC (that also pass p < 0.05 and FC > 1)
filtered_up <- subset(deg_cell_subset_D0_D34, avg_log2FC > 1 & p_val_adj < 0.05)
top5_genes <- rownames(head(filtered_up[order(-filtered_up$avg_log2FC), ], 5))

# Now plot
EnhancedVolcano(deg_cell_subset_D0_D34,
                lab = all_genes,  # label all genes (internally matched to selectLab)
                selectLab = top5_genes,  # these are the only ones that will be shown
                x = "avg_log2FC",
                y = "p_val_adj",
                title = "Macrophages: D3/4 vs D0",
                pCutoff = 0.05,
                FCcutoff = 1,
                drawConnectors = TRUE)



# Fibroblasts
# Subset fibroblasts
subset_cells <- WhichCells(integrated, idents = "Fibroblsts")
subset_data <- subset(integrated, cells = subset_cells)

# Re-run SCTransform
DefaultAssay(subset_data) <- "RNA"
subset_data <- SCTransform(subset_data, new.assay.name = "SCT", verbose = FALSE)

# Use RNA assay for DE analysis
DefaultAssay(subset_data) <- "RNA"
subset_data <- NormalizeData(subset_data, verbose = FALSE)

# Set timepoint identities
Idents(subset_data) <- subset_data$timepoint

# Find DEGs between D3/4 and D0
deg_cell_subset_D0_D34 <- FindMarkers(subset_data, 
                                      ident.1 = "D3/4", 
                                      ident.2 = "D0", 
                                      assay = "RNA", 
                                      test.use = "wilcox")

# Get top 5 upregulated DEGs
filtered_up <- subset(deg_cell_subset_D0_D34, avg_log2FC > 1 & p_val_adj < 0.05)
top5_genes <- rownames(head(filtered_up[order(-filtered_up$avg_log2FC), ], 5))

# Volcano plot
EnhancedVolcano(deg_cell_subset_D0_D34,
                lab = rownames(deg_cell_subset_D0_D34),
                selectLab = top5_genes,
                x = "avg_log2FC",
                y = "p_val_adj",
                title = "Fibroblasts: D3/4 vs D0",
                pCutoff = 0.05,
                FCcutoff = 1,
                drawConnectors = TRUE)
