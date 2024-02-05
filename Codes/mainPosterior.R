# delete work space
rm(list = ls(all = TRUE))
graphics.off()

library(Seurat)
library(SeuratData)
library(ggpubr)
library(ggplot2)
library(patchwork)
library(dplyr)
library(MOFA2)
library(tidyverse)
library(cowplot)
library(magrittr)
library(SpatialPCA)
library(spacexr)
library(kernlab)
library(MAST)
# library(fdrtool)
library(qvalue)
library(gprofiler2)
library(forcats)
library(monocle3)
library(SeuratWrappers)
library(ggridges)

source("helperFunctions.R")
source("CoCoST.R")

anterior <- LoadData("stxBrain", type = "anterior1")
anterior <- SCTransform(anterior, assay = "Spatial", verbose = FALSE)
SpatialFeaturePlot(anterior, features = "nCount_Spatial") + theme(legend.position = "right")

posterior <- LoadData("stxBrain", type = "posterior1")
posterior <- SCTransform(posterior, assay = "Spatial", verbose = FALSE)
SpatialFeaturePlot(posterior, features = "nCount_Spatial") + theme(legend.position = "right")

# CoCoST
fdata <- posterior@assays[["SCT"]]@scale.data
# flocation <- abnormalTissue@images[["slice1"]]@coordinates[,c(2,3)] 
flocation <- GetTissueCoordinates(posterior)

bdata <- anterior@assays[["SCT"]]@scale.data
# blocation <- normalTissue@images[["slice1"]]@coordinates[,c(2,3)]
blocation <- GetTissueCoordinates(anterior)

# construct affinity matrices
rbf <- laplacedot(sigma = 0.50)

fKernel <- kernelMatrix(rbf, t(fdata))
Wf <- fKernel@.Data

bKernel <- kernelMatrix(rbf, t(bdata))
Wb <- bKernel@.Data

# Extract contrastive features
para <- 0.10
Dim <- 50
CoCo <- CoCoST(t(fdata),Wf,t(bdata),Wb,para,Dim)

fCoCo <- CreateDimReducObject(
  embeddings = CoCo[["fgComponents"]],
  loadings = CoCo[["projMatrix"]],
  stdev = numeric(),
  key = "CoCoST_",
  assay = "SCT"
)

rownames(fCoCo@feature.loadings) <- rownames(fdata)
posterior@reductions[["CoCo_ST"]] <- fCoCo
posterior <- RunUMAP(posterior, reduction = "CoCo_ST", dims = 1:50, 
                    n.components = 5, reduction.name = "CoCo_UMAP")
posterior <- FindNeighbors(posterior, reduction = "CoCo_ST", dims = 1:10)
posterior <- FindClusters(posterior, verbose = TRUE, resolution = 0.8)
posterior@meta.data[["CoCo_clusters"]]  <- posterior@meta.data[["SCT_snn_res.0.8"]]

# cols <- c("#E1BD6D","deepskyblue1","#7A3A9A","#ED0000FF","#0B775E","#ff00ff","#7B556C","#44B05B")
DimPlot(posterior, reduction = "CoCo_UMAP", label = FALSE)
SpatialDimPlot(posterior, label = FALSE, label.size = 3, group.by = "CoCo_clusters", pt.size.factor = 2.2) 

SpatialFeaturePlot(posterior, features = c("CoCoST_1", "CoCoST_2", "CoCoST_3", "CoCoST_4", "CoCoST_5"), ncol = 5, 
                   pt.size.factor = 2.1)

# Plot gene loading for foreground CoCoST components
refComponents <- c(1:5)
loadings <- as.data.frame(posterior@reductions[["CoCo_ST"]]@feature.loadings[,refComponents])
nfea <- nrow(loadings)
whichComp <- c(1:5)
numGenes <- 20

topGenes <- getTopGenes(loadings,refComponents,numGenes)

geneList <- topGenes[["topGenes"]]
pltTopG <- plotLoadings(geneList)
pltTopG[[1]]

compTopG <- topGenes[["compTopFirstGene"]]
SpatialFeaturePlot(posterior, features = compTopG, ncol = 5, pt.size.factor = 2.2)


