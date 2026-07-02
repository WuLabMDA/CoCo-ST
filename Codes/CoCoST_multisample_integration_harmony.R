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
library(harmony)


options(future.globals.maxSize = 1024 * 1024 * 10000)

## Load Data

normalTissue <- Load10X_Spatial(
  data.dir = ".../53430/outs",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", 
  filter.matrix = TRUE, 
  to.upper = FALSE
)

abnormalTissue1 <- Load10X_Spatial(".../26934/outs",
                                   filename = "filtered_feature_bc_matrix.h5",
                                   assay = "Spatial", 
                                   filter.matrix = TRUE, 
                                   to.upper = FALSE)
abnormalTissue2 <- Load10X_Spatial(".../53432/outs",
                                   filename = "filtered_feature_bc_matrix.h5",
                                   assay = "Spatial", 
                                   filter.matrix = TRUE, 
                                   to.upper = FALSE)
abnormalTissue3 <- Load10X_Spatial(".../53431/outs",
                                   filename = "filtered_feature_bc_matrix.h5",
                                   assay = "Spatial", 
                                   filter.matrix = TRUE, 
                                   to.upper = FALSE)


## User inputs

ngenes <- 4000
n_coco_genes <- 3000
coco_dim <- 50
sigma_use <- 0.50
coco_para <- 0.10
cluster_resolution <- 0.19

sample_names <- c(
  "normal" = "Normal",
  "abnormal1" = "MLP7",
  "abnormal2" = "MLP4",
  "abnormal3" = "MLP3"
)


## Merge normal + abnormal tissues

combData <- merge(
  normalTissue,
  y = list(
    abnormalTissue1,
    abnormalTissue2,
    abnormalTissue3
  ),
  add.cell.ids = c(
    "normal",
    "abnormal1",
    "abnormal2",
    "abnormal3"
  ),
  project = "Spatial"
)


## Add sample and condition metadata

combData$sample_id <- NA
combData$condition <- NA

combData$sample_id[grepl("^normal_", colnames(combData))] <- "Normal"
combData$sample_id[grepl("^abnormal1_", colnames(combData))] <- "MLP7"
combData$sample_id[grepl("^abnormal2_", colnames(combData))] <- "MLP4"
combData$sample_id[grepl("^abnormal3_", colnames(combData))] <- "MLP3"

combData$condition[grepl("^normal_", colnames(combData))] <- "Normal"
combData$condition[grepl("^abnormal", colnames(combData))] <- "Abnormal"

combData$sample_id <- factor(
  combData$sample_id,
  levels = c("Normal", "MLP7", "MLP4", "MLP3")
)

combData$condition <- factor(
  combData$condition,
  levels = c("Normal", "Abnormal")
)


## SCTransform normalization

combData <- SCTransform(
  combData,
  assay = "Spatial",
  variable.features.n = ngenes,
  verbose = FALSE
)


## Split normal/background and abnormal/foreground data

normal_cells <- colnames(combData)[combData$condition == "Normal"]
abnormal_cells <- colnames(combData)[combData$condition == "Abnormal"]

sct_scaled <- combData@assays[["SCT"]]@scale.data

genes_use <- rownames(sct_scaled)[
  seq_len(min(n_coco_genes, nrow(sct_scaled)))
]

bdata <- sct_scaled[
  genes_use,
  normal_cells,
  drop = FALSE
]

fdata <- sct_scaled[
  genes_use,
  abnormal_cells,
  drop = FALSE
]


## Run CoCo-ST

rbf <- laplacedot(
  sigma = sigma_use
)

Wf <- kernelMatrix(
  rbf,
  t(fdata)
)@.Data

Wb <- kernelMatrix(
  rbf,
  t(bdata)
)@.Data

CoCo_model <- CoCoST(
  t(fdata),
  Wf,
  t(bdata),
  Wb,
  para = coco_para,
  Dim = coco_dim
)


## Create CoCo-ST embeddings for abnormal cells

abnormal_coco_embeddings <- CoCo_model$fgComponents

rownames(abnormal_coco_embeddings) <- colnames(fdata)

colnames(abnormal_coco_embeddings) <- paste0(
  "CoCoST_",
  seq_len(ncol(abnormal_coco_embeddings))
)


## Project normal cells into same CoCo-ST space

normal_coco_embeddings <- t(bdata) %*% CoCo_model$projMatrix

rownames(normal_coco_embeddings) <- colnames(bdata)

colnames(normal_coco_embeddings) <- paste0(
  "CoCoST_",
  seq_len(ncol(normal_coco_embeddings))
)


## Combine normal + abnormal CoCo-ST embeddings

all_coco_embeddings <- rbind(
  normal_coco_embeddings,
  abnormal_coco_embeddings
)

all_coco_embeddings <- all_coco_embeddings[
  colnames(combData),
  ,
  drop = FALSE
]


# Add CoCo-ST reduction to combined Seurat object

CoCo_combined <- CreateDimReducObject(
  embeddings = all_coco_embeddings,
  loadings = CoCo_model$projMatrix,
  key = "CoCoST_",
  assay = "SCT"
)

rownames(CoCo_combined@feature.loadings) <- rownames(fdata)

combData[["CoCoST"]] <- CoCo_combined

## Harmony correction on CoCo-ST components

combData <- RunHarmony(
  object = combData,
  group.by.vars = "sample_id",
  reduction = "CoCoST",
  assay.use = "SCT",
  reduction.save = "CoCoST_harmony",
  dims.use = 1:coco_dim
)


## UMAP and clustering on Harmony-corrected CoCo-ST

combData <- RunUMAP(
  combData,
  reduction = "CoCoST_harmony",
  dims = 1:coco_dim,
  reduction.name = "CoCoST_Harmony_UMAP"
)

combData <- FindNeighbors(
  combData,
  reduction = "CoCoST_harmony",
  dims = 1:coco_dim
)

combData <- FindClusters(
  combData,
  resolution = cluster_resolution
)

combData$hcoco_harmony_domains <- combData$SCT_snn_res.0.19


sample_cols <- c(
  "Normal" = "#555555",
  "MLP7" = "#cc440a",
  "MLP4" = "#1559b3",
  "MLP3" = "#20821d"
)

condition_cols <- c(
  "Normal" = "#555555",
  "Abnormal" = "#D73027"
)

domain_cols <- c(
  "0" = "#E1BD6D",
  "1" = "#ED0000FF",
  "2" = "#44B05B",
  "3" = "#a18e23",
  "4" = "#0B775E",
  "5" = "#0e7ec4",
  "6" = "#ff00ff",
  "7" = "#7A3A9A",
  "8" = "#7B556C",
  "9" = "#A6CEE3",
  "10" = "#FB9A99",
  "11" = "#FDBF6F",
  "12" = "#CAB2D6",
  "13" = "#FFFF99",
  "14" = "#B15928",
  "15" = "#1F78B4"
)


## UMAP visualization

DimPlot(
  combData,
  reduction = "CoCoST_Harmony_UMAP",
  group.by = "sample_id",
  cols = sample_cols
)

DimPlot(
  combData,
  reduction = "CoCoST_Harmony_UMAP",
  group.by = "condition",
  cols = condition_cols
)

DimPlot(
  combData,
  reduction = "CoCoST_Harmony_UMAP",
  group.by = "hcoco_harmony_domains",
  label = TRUE,
  cols = domain_cols
)


## Map CoCoST-Harmony domains back to original objects

meta_all <- data.frame(
  hcoco_harmony_domains = combData$hcoco_harmony_domains,
  sample_id = combData$sample_id,
  condition = combData$condition
)

rownames(meta_all) <- colnames(combData)


## Function to map metadata

map_hcoco_harmony_domains <- function(
    seu_obj,
    prefix,
    meta_df
) {
  
  meta_sub <- meta_df[
    grep(
      paste0("^", prefix, "_"),
      rownames(meta_df)
    ),
    ,
    drop = FALSE
  ]
  
  rownames(meta_sub) <- sub(
    paste0("^", prefix, "_"),
    "",
    rownames(meta_sub)
  )
  
  common_cells <- intersect(
    colnames(seu_obj),
    rownames(meta_sub)
  )
  
  meta_sub <- meta_sub[
    common_cells,
    ,
    drop = FALSE
  ]
  
  seu_obj <- AddMetaData(
    seu_obj,
    metadata = meta_sub
  )
  
  return(seu_obj)
}


## Apply metadata to each original tissue

normalTissue <- map_hcoco_harmony_domains(
  seu_obj = normalTissue,
  prefix = "normal",
  meta_df = meta_all
)

abnormalTissue1 <- map_hcoco_harmony_domains(
  seu_obj = abnormalTissue1,
  prefix = "abnormal1",
  meta_df = meta_all
)

abnormalTissue2 <- map_hcoco_harmony_domains(
  seu_obj = abnormalTissue2,
  prefix = "abnormal2",
  meta_df = meta_all
)

abnormalTissue3 <- map_hcoco_harmony_domains(
  seu_obj = abnormalTissue3,
  prefix = "abnormal3",
  meta_df = meta_all
)


## Spatial plots: CoCoST-Harmony domains

SpatialDimPlot(
  normalTissue,
  group.by = "hcoco_harmony_domains",
  label = FALSE,
  label.size = 3,
  pt.size.factor = 3.2,
  stroke = 1,
  cols = domain_cols
) +
  ggtitle("Normal tissue: CoCoST-Harmony domains")

SpatialDimPlot(
  abnormalTissue1,
  group.by = "hcoco_harmony_domains",
  label = FALSE,
  label.size = 3,
  pt.size.factor = 3.2,
  stroke = 1,
  cols = domain_cols
) +
  ggtitle("Abnormal tissue 1: CoCoST-Harmony domains")

SpatialDimPlot(
  abnormalTissue2,
  group.by = "hcoco_harmony_domains",
  label = FALSE,
  label.size = 3,
  pt.size.factor = 3.2,
  stroke = 1,
  cols = domain_cols
) +
  ggtitle("Abnormal tissue 2: CoCoST-Harmony domains")

SpatialDimPlot(
  abnormalTissue3,
  group.by = "hcoco_harmony_domains",
  label = FALSE,
  label.size = 3,
  pt.size.factor = 3.2,
  stroke = 1,
  cols = domain_cols
) +
  ggtitle("Abnormal tissue 3: CoCoST-Harmony domains")


