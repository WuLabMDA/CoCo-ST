rm(list = ls(all = TRUE))
graphics.off()

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)
library(kernlab)
library(monocle3)
library(SeuratWrappers)
library(ggridges)
library(CoCoST)
library(kBET)


options(future.globals.maxSize = 1024 * 1024 * 10000)

## 1. Load Data

normalTissue <- Load10X_Spatial(
  data.dir = "Y:/Projects/Colab/Colab Bo/Data/Old directory/Processed 129S4 Urethane model/53430/outs",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", 
  filter.matrix = TRUE, 
  to.upper = FALSE
)

abnormalTissue1 <- Load10X_Spatial("Y:/Projects/Colab/Colab Bo/Data/Old directory/Processed 129S4 Urethane model/26934/outs",
                                   filename = "filtered_feature_bc_matrix.h5",
                                   assay = "Spatial", 
                                   filter.matrix = TRUE, 
                                   to.upper = FALSE)
abnormalTissue2 <- Load10X_Spatial("Y:/Projects/Colab/Colab Bo/Data/Old directory/Processed 129S4 Urethane model/53432/outs",
                                   filename = "filtered_feature_bc_matrix.h5",
                                   assay = "Spatial", 
                                   filter.matrix = TRUE, 
                                   to.upper = FALSE)
abnormalTissue3 <- Load10X_Spatial("Y:/Projects/Colab/Colab Bo/Data/Old directory/Processed 129S4 Urethane model/53431/outs",
                                   filename = "filtered_feature_bc_matrix.h5",
                                   assay = "Spatial", 
                                   filter.matrix = TRUE, 
                                   to.upper = FALSE)

## 2. Merge & Normalize

combData <- merge(
  normalTissue,
  y = list(abnormalTissue1, abnormalTissue2, abnormalTissue3),
  add.cell.ids = c("normal","abnormal1","abnormal2","abnormal3"),
  project = "Spatial"
)

ngenes <- 4000

combData <- SCTransform(
  combData,
  assay = "Spatial",
  variable.features.n = ngenes,
  verbose = FALSE
)

## 3. Split Normal vs Abnormal

n_normal <- ncol(normalTissue)

bdata <- combData@assays[["SCT"]]@scale.data[1:3000, 1:n_normal]
fdata <- combData@assays[["SCT"]]@scale.data[1:3000, (n_normal+1):ncol(combData)]

## 4. Create Abnormal Seurat Object

combAbnormal <- CreateSeuratObject(
  counts = combData@assays[["SCT"]]@data[, (n_normal+1):ncol(combData)],
  assay = "Spatial"
)

combAbnormal <- SCTransform(
  combAbnormal,
  assay = "Spatial",
  variable.features.n = ngenes,
  verbose = FALSE
)

## 5. CoCo-ST

rbf <- laplacedot(sigma = 0.50)

Wf <- kernelMatrix(rbf, t(fdata))@.Data
Wb <- kernelMatrix(rbf, t(bdata))@.Data

CoCo_model <- CoCoST(
  t(fdata),
  Wf,
  t(bdata),
  Wb,
  para = 0.10,
  Dim = 50
)

CoCo <- CreateDimReducObject(
  embeddings = CoCo_model$fgComponents,
  loadings = CoCo_model$projMatrix,
  key = "CoCoST_",
  assay = "SCT"
)

rownames(CoCo@feature.loadings) <- rownames(fdata)

combAbnormal[["CoCoST"]] <- CoCo

## 6. Metadata

combAbnormal$batch <- factor(c(
  rep("MLP7", ncol(abnormalTissue1)),
  rep("MLP4", ncol(abnormalTissue2)),
  rep("MLP3", ncol(abnormalTissue3))
))

## 7. Clustering

combAbnormal <- RunUMAP(
  combAbnormal,
  reduction = "CoCoST",
  dims = 1:50,
  n.components = 50,
  reduction.name = "CoCo_UMAP"
)

combAbnormal <- FindNeighbors(combAbnormal, reduction = "CoCoST", dims = 1:50)
combAbnormal <- FindClusters(combAbnormal, resolution = 0.19)

combAbnormal$hcoco_clusters <- combAbnormal$SCT_snn_res.0.19

## 8. Visualization

cols <- c(
  "0"="#E1BD6D","1"="#ED0000FF","2"="#44B05B","3"="#a18e23",
  "4"="#0B775E","5"="#0e7ec4","6"="#ff00ff","7"="#7A3A9A","8"="#7B556C"
)

cols2 <- c("#cc440a","#1559b3","#20821d")

DimPlot(combAbnormal, reduction="CoCo_UMAP", group.by="batch", cols=cols2)
DimPlot(combAbnormal, reduction="CoCo_UMAP", cols=cols)

## 9. Map Clusters Back

meta <- data.frame(
  hcoco_clusters = combAbnormal$hcoco_clusters
)

## Function to map metadata safely

map_clusters <- function(seu_obj, prefix, meta_df) {
  
  meta_sub <- meta_df[grep(paste0("^", prefix, "_"), rownames(meta_df)), , drop = FALSE]
  
  rownames(meta_sub) <- sub(paste0("^", prefix, "_"), "", rownames(meta_sub))
  
  common_cells <- intersect(colnames(seu_obj), rownames(meta_sub))
  
  meta_sub <- meta_sub[common_cells, , drop = FALSE]
  
  seu_obj <- AddMetaData(seu_obj, metadata = meta_sub)
  
  return(seu_obj)
}

abnormalTissue1 <- map_clusters(abnormalTissue1, "abnormal1", meta)
abnormalTissue2 <- map_clusters(abnormalTissue2, "abnormal2", meta)
abnormalTissue3 <- map_clusters(abnormalTissue3, "abnormal3", meta)


## 10. Spatial Plots

SpatialDimPlot(abnormalTissue1, label = FALSE, label.size = 3, group.by = "hcoco_clusters", pt.size.factor = 3.2, stroke = 1, cols = cols)
SpatialDimPlot(abnormalTissue2, label = FALSE, label.size = 3, group.by = "hcoco_clusters", pt.size.factor = 3.2, stroke = 1, cols = cols)
SpatialDimPlot(abnormalTissue3, label = FALSE, label.size = 3, group.by = "hcoco_clusters", pt.size.factor = 3.2, stroke = 1, cols = cols)

## 11. kBET test for batch correction

umap_embeddings <- Embeddings(combAbnormal, "CoCo_UMAP")
batch_labels <- combAbnormal$batch

kbet_results <- kBET(umap_embeddings, batch_labels)

summary(kbet_results$results)

plot.data <- data.frame(
  class = rep(c('observed', 'expected'), each = length(kbet_results$stats$kBET.observed)), 
  data = c(kbet_results$stats$kBET.observed, kbet_results$stats$kBET.expected)
)

ggplot(plot.data, aes(class, data)) +
  geom_boxplot() +
  labs(x = 'Test', y = 'Rejection rate', title = 'kBET test') +
  theme_bw() +
  scale_y_continuous(limits = c(0, 1)) +
  theme(
    axis.text = element_text(face = "bold", size = 12),  
    axis.title = element_text(face = "bold", size = 12), 
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),  
    panel.border = element_blank(),  
    panel.grid = element_line(color = "grey90"),  
    axis.line = element_line(color = "black", size = 0.8),  
    axis.ticks = element_line(color = "black")  
  )
