# delete work space
rm(list = ls(all = TRUE))
graphics.off()

setwd("C:/Projects/Contrastive PCA/ST/Mouse Data (Bo)/Codes")

library(spacexr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(kernlab)
library(MAST)
# library(fdrtool)
library(qvalue)
library(gprofiler2)
library(tidyverse)
library(forcats)
library(monocle3)
library(Seurat)
library(dplyr)
library(monocle3)
library(SeuratWrappers)
library(ggridges)

source("helperFunctions.R")
source("CoCoST.R")

# load background data
normalTissue <- Load10X_Spatial(
  data.dir = "X:/maminu/Processed 129S4 Urethane model/53430/outs", 
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", # specify name of the initial assay
  filter.matrix = TRUE, 
  to.upper = FALSE
)
normalTissue <- SCTransform(normalTissue, assay = "Spatial", verbose = FALSE)
SpatialFeaturePlot(normalTissue, features = "nCount_Spatial") + theme(legend.position = "right")

# load foreground data
abnormalTissue <- Load10X_Spatial(
  data.dir = "X:/maminu/Processed 129S4 Urethane model/53433/outs", 
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", # specify name of the initial assay
  filter.matrix = TRUE, 
  to.upper = FALSE
)
abnormalTissue <- SCTransform(abnormalTissue, assay = "Spatial", verbose = FALSE)
SpatialFeaturePlot(abnormalTissue, features = "nCount_Spatial") + theme(legend.position = "right")

# graph contrastive PCA
fdata <- abnormalTissue@assays[["SCT"]]@scale.data
# flocation <- abnormalTissue@images[["slice1"]]@coordinates[,c(2,3)] 
flocation <- GetTissueCoordinates(abnormalTissue)

bdata <- normalTissue@assays[["SCT"]]@scale.data
# blocation <- normalTissue@images[["slice1"]]@coordinates[,c(2,3)]
blocation <- GetTissueCoordinates(normalTissue)

# construct affinity matrices
rbf <- laplacedot(sigma = 0.50)

fKernel <- kernelMatrix(rbf, t(fdata))
Wf <- fKernel@.Data

bKernel <- kernelMatrix(rbf, t(bdata))
Wb <- bKernel@.Data

# Extract contrastive features
para <- 0.10
Dim <- 20
gCPCA <- CoCoST(t(fdata),Wf,t(bdata),Wb,para,Dim)

CoCo <- CreateDimReducObject(
  embeddings = gCPCA[["fgComponents"]],
  loadings = gCPCA[["projMatrix"]],
  stdev = numeric(),
  key = "CoCoST_",
  assay = "SCT"
)

rownames(CoCo@feature.loadings) <- rownames(fdata)
abnormalTissue@reductions[["CoCoST"]] <- CoCo

abnormalTissue <- RunUMAP(abnormalTissue, reduction = "CoCoST", dims = 1:10, 
                          n.components = 5, reduction.name = "CoCo_UMAP")
abnormalTissue <- FindNeighbors(abnormalTissue, reduction = "CoCo_UMAP", dims = 1:5)
abnormalTissue <- FindClusters(abnormalTissue, verbose = TRUE, resolution = 0.05)
abnormalTissue@meta.data[["coco_clusters"]]  <- abnormalTissue@meta.data[["SCT_snn_res.0.05"]]

cols <- c("#E1BD6D","deepskyblue1","#7A3A9A","#ED0000FF","#0B775E","#ff00ff","#7B556C","#44B05B")
DimPlot(abnormalTissue, reduction = "CoCo_UMAP", label = FALSE, cols = cols)

SpatialDimPlot(abnormalTissue, label = FALSE, label.size = 3, group.by = "coco_clusters", pt.size.factor = 2) +  
  scale_fill_manual(values=c("#E1BD6D","deepskyblue1","#7A3A9A","#ED0000FF","#0B775E","#ff00ff","#7B556C","#44B05B"))

abnormalTissue <- RenameIdents(object = abnormalTissue, `0` = "SD 1", `1` = "SD 3", `2` = "SD 2", `3` = "SD 4",
                               `4` = "SD 5", `5` = "SD 6")

SpatialDimPlot(abnormalTissue, label = FALSE, label.size = 3, pt.size.factor = 2) +  
  scale_fill_manual(values=c("#E1BD6D","deepskyblue1","#7A3A9A","#ED0000FF","#0B775E","#ff00ff","#7B556C","#44B05B"))

# Find differentially expressed features
abnormalClust1.markers <- FindMarkers(abnormalTissue, ident.1 = "SD 1", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
abnormalClust2.markers <- FindMarkers(abnormalTissue, ident.1 = "SD 2", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
abnormalClust3.markers <- FindMarkers(abnormalTissue, ident.1 = "SD 3", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
abnormalClust4.markers <- FindMarkers(abnormalTissue, ident.1 = "SD 4", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
abnormalClust5.markers <- FindMarkers(abnormalTissue, ident.1 = "SD 5", logfc.threshold = 0.25, test.use = "LR", only.pos = FALSE)
abnormalClust6.markers <- FindMarkers(abnormalTissue, ident.1 = "SD 6", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

library(EnhancedVolcano)
EnhancedVolcano(abnormalClust5.markers,
                lab = rownames(abnormalClust5.markers),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                selectLab = c(rownames(abnormalClust5.markers)[c(1:10)]),
                xlab = bquote(~Log[2]~ 'fold change'),
                # pCutoff = 10e-14,
                # FCcutoff = 2.0,
                pointSize = 3.0,
                labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')

EnhancedVolcano(abnormalClust5.markers,
                lab = rownames(abnormalClust5.markers),
                selectLab = c(rownames(abnormalClust5.markers)[c(1:10)],"Trf"),
                x = 'avg_log2FC',
                y = 'p_val_adj')

my_levels <- c("SD 1", "SD 2", "SD 3", "SD 4", "SD 5", "SD 6")
Idents(abnormalTissue) <- factor(Idents(abnormalTissue), levels= my_levels)
VlnPlot(abnormalTissue, features = c(rownames(abnormalClust1.markers)[1],rownames(abnormalClust2.markers)[1],
                                     rownames(abnormalClust3.markers)[1],rownames(abnormalClust4.markers)[1],
                                     rownames(abnormalClust5.markers)[1],rownames(abnormalClust6.markers)[2]),
        pt.size = 0, cols = c("#E1BD6D","#7A3A9A","deepskyblue1","#ED0000FF","#0B775E","#ff00ff"), ncol = 6)

DimPlot(abnormalTissue, reduction = "CoCo_UMAP", label = FALSE, cols = c("#E1BD6D","#7A3A9A","deepskyblue1","#ED0000FF","#0B775E","#ff00ff"))

FeaturePlot(abnormalTissue, features = c(rownames(abnormalClust1.markers)[1],rownames(abnormalClust2.markers)[1],
                                         rownames(abnormalClust3.markers)[1],rownames(abnormalClust4.markers)[1],
                                         rownames(abnormalClust5.markers)[1],rownames(abnormalClust6.markers)[2]),
            cols = c("lightgray","darkgray", "red"))

SpatialDimPlot(abnormalTissue, label = FALSE, label.size = 3, group.by = "coco_clusters", pt.size.factor = 2) +  
  scale_fill_manual(values=c("#E1BD6D","deepskyblue1","#7A3A9A","#ED0000FF","#0B775E","#ff00ff","#7B556C","#44B05B"),
                    labels = c("Normal lung","Fibrotic/Scar tissue","Adjacent normal","Bronchus/Alveoli",
                               "Adenoma","Membrane"))

# plot refine tissue maps
plotCluster(location=flocation,clusterlabel=abnormalTissue@meta.data[["coco_clusters"]],pointsize=1.5,
            title_in=paste0("CoCo-ST"),color_in=cols)

# SpatialDimPlot(abnormalTissue, cells.highlight = CellsByIdentities(object = abnormalTissue, idents = c(0,1,2,3,4,5)), 
#                facet.highlight = TRUE, ncol = 3)

SpatialFeaturePlot(abnormalTissue, features = c("CoCoST_1", "CoCoST_2", "CoCoST_3", "CoCoST_4", "CoCoST_5"), ncol = 5, 
                   pt.size.factor = 2)

# Plot gene loading for foreground gCPCA components
refComponents <- c(1:5)
loadings <- as.data.frame(abnormalTissue@reductions[["CoCoST"]]@feature.loadings[,refComponents])
nfea <- nrow(loadings)
whichComp <- c(1:5)
numGenes <- 50

topGenes <- getTopGenes(loadings,refComponents,numGenes)

geneList <- topGenes[["topGenes"]]
pltTopG <- plotLoadings(geneList)
pltTopG[[1]]

compTopG <- topGenes[["compTopFirstGene"]]
SpatialFeaturePlot(abnormalTissue, features = compTopG, ncol = 5)
SpatialFeaturePlot(abnormalTissue, features = c("Ctsh","Cxcl15","Slc34a2","Trf"), ncol = 5, pt.size.factor = 2.5)
SpatialFeaturePlot(abnormalTissue, features = c("Trf","Ctsh"), ncol = 5, pt.size.factor = 2.5)

# perform pathway enrichment analysis of component top genes
bg <- rownames(fdata)
genePathways <- getPathways(geneList,bg)
pltGO <- plotPathways(genePathways)
# pltGO[[1]][["data"]][["term_name"]][10] <- "oxidoreductase activity, acting on paired donors"
pltGO[[5]]


# perform deconvolution
# set up reference
ref <- readRDS("X:/maminu/Bo/Single cell RNA Seq data/129S4/129S4_Major_CellTypes.rds")
ref <- UpdateSeuratObject(ref)
subref <- subset(ref, model == "Urethane")
Idents(subref) <- "Major_CellType"
subref@meta.data[["Major_CellType"]] <- factor(subref@meta.data[["Major_CellType"]], 
                                               levels = c("Endothelial cells","Epithelial cells","Fibroblasts",
                                                          "Macrophages","cDC","Prolif_Macro","B cells","T cells",
                                                          "Proliferating T cells","pDC","Neutrophils","Plasma cells",
                                                          "Monocytes","NK cells"))

col2 = c("plum1","tomato","#762A83","deepskyblue1","#C4961A","#ff00ff","#DC0000FF","#4E84C4","chartreuse3",
                "#D16103","#58593FFF","lightblue1","#068105","yellow","#4E2A1E","#C3D7A4","black")
                
DimPlot(subref, group.by = "Major_CellType", label = F, cols = col2)

# markerGenes <- c("Cd34","Cd44","Cd13","Ccr5","Cd8a","Ccr7","Cd11c")
# DotPlot(subref, features = markerGenes)

# extract information to pass to the RCTD Reference function
counts <- subref@assays[["RNA"]]@counts
cluster <- as.factor(subref@meta.data[["Major_CellType"]])
names(cluster) <- colnames(subref)
nUMI <- subref@meta.data[["nCount_RNA"]]
names(nUMI) <- colnames(subref)
reference <- Reference(counts, cluster, nUMI)

# set up query with the RCTD function SpatialRNA
STcounts <- abnormalTissue@assays[["Spatial"]]@counts
coords <- GetTissueCoordinates(abnormalTissue)
colnames(coords) <- c("x", "y")
coords[is.na(colnames(coords))] <- NULL
query <- SpatialRNA(coords, STcounts, colSums(STcounts))

RCTD <- create.RCTD(query, reference, max_cores = 1)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet") # previously used doublet mode
abnormalTissue <- AddMetaData(abnormalTissue, metadata = RCTD@results$results_df)
abnormalTissue@meta.data[["second_type"]] <- factor(abnormalTissue@meta.data[["second_type"]], 
                                                    levels = c("Endothelial cells","Epithelial cells","Fibroblasts",
                                                               "Macrophages","cDC","Prolif_Macro","B cells","T cells",
                                                               "Proliferating T cells","pDC","Neutrophils","Plasma cells",
                                                               "Monocytes","NK cells"))
abnormalTissue@meta.data[["first_type"]] <- factor(abnormalTissue@meta.data[["first_type"]], 
                                                    levels = c("Endothelial cells","Epithelial cells","Fibroblasts",
                                                               "Macrophages","cDC","Prolif_Macro","B cells","T cells",
                                                               "Proliferating T cells","pDC","Neutrophils","Plasma cells",
                                                               "Monocytes","NK cells"))

SpatialDimPlot(abnormalTissue, group.by = "first_type") + scale_fill_manual(values=col2)
SpatialDimPlot(abnormalTissue, group.by = "second_type") + scale_fill_manual(values=col2)


# add deconvolution result to st data as a new assay
tempabnormal <- abnormalTissue
tempabnormal@assays[["RCTD"]] <- CreateAssayObject(data = t(RCTD@results[["weights"]]))

if (length(tempabnormal@assays$RCTD@key) == 0) {
  tempabnormal@assays$RCTD@key = "rctd_"
}

DefaultAssay(tempabnormal) <- "RCTD"
int1 <- intersect(colnames(tempabnormal),rownames(tempabnormal@images[["slice1"]]@coordinates))
tempabnormal@images[["slice1"]]@coordinates <- tempabnormal@images[["slice1"]]@coordinates[int1,]
SpatialFeaturePlot(tempabnormal, features = rownames(tempabnormal)[1:4], pt.size.factor = 2, ncol = 4)



#--------- plot stack barplot for cell types inference at spatial domains ----------
int <- intersect(colnames(abnormalTissue),rownames(abnormalTissue@images[["slice1"]]@coordinates))
abnormalTissue@images[["slice1"]]@coordinates <- abnormalTissue@images[["slice1"]]@coordinates[int1,]

method="gCPCA"
metadata_RCTD = abnormalTissue@meta.data
metadata_RCTD$celltype <- abnormalTissue@meta.data[["second_type"]]
metadata_RCTD$celltype <- as.character(metadata_RCTD$celltype)
metadata_RCTD$second_type <-NULL
metadata_RCTD$first_class <- NULL
metadata_RCTD$second_class <- NULL
metadata_RCTD$min_score <- NULL
metadata_RCTD$singlet_score <- NULL
metadata_RCTD$conv_all <- NULL
metadata_RCTD$conv_doublet <- NULL
metadata_RCTD$cellid <- rownames(metadata_RCTD)
metadata_RCTD$celltype <- as.factor(metadata_RCTD$celltype)

percentage = matrix(0,length(unique(metadata_RCTD$gcpca_clusters)),
                    (length(unique(metadata_RCTD$celltype))-1))
for(k in 1:length(unique(metadata_RCTD$gcpca_clusters))){
  metadata_sub = metadata_RCTD[which(metadata_RCTD$gcpca_clusters==k-1 ),]
  match_type = metadata_sub$celltype
  percentage[k,] = round(unlist(table(match_type))/dim(metadata_RCTD)[1]*100,2)
}

celltypenames =  names(table(metadata_RCTD$celltype))
datt2=preparedata(percentage, celltypenames)
datt2$cluster_vec = as.character(datt2$cluster_vec)
datt2$cluster_vec = factor(datt2$cluster_vec, levels=c(paste0("Cluster",c(1,2,3,4,5,6))), 
                           labels = c(paste0("Cluster",c(1,2,3,4,5,6))),order=T)
datt2$CellType <- factor(datt2$CellType, 
                         levels = c("Endothelial cells","Epithelial cells","Fibroblasts",
                                    "Macrophages","cDC","Prolif_Macro","B cells","T cells",
                                    "Proliferating T cells","pDC","Neutrophils","Plasma cells",
                                    "Monocytes","NK cells"))
makefigure(datt2)+ggtitle(paste0(method))

# plot cell type proportion in each cluster
metatable <- table(abnormalTissue@meta.data[["gcpca_clusters"]],abnormalTissue@meta.data[["second_type"]])
metatable_proportion <- matrix(0 ,dim(metatable)[1],dim(metatable)[2])
for(celltypeid in 1:dim(metatable)[2]){
  metatable_proportion[,celltypeid] <- metatable[,celltypeid]/sum(metatable[,celltypeid])
}
colnames(metatable_proportion) <- colnames(metatable)
# rownames(metatable_proportion) <- paste(rep("Cluster",nrow(metatable_proportion)),c(1:nrow(metatable_proportion)))
rownames(metatable_proportion) <- c("Normal lung","Fibrotic/Scar tissue","Adjacent normal","Bronchus/Alveoli",
                                    "Adenoma","Membrane")
# regionname <- paste(rep("Cluster",nrow(metatable_proportion)),c(1:nrow(metatable_proportion)))
regionname <- rownames(metatable_proportion)
# cbp_spatialpca <- c("#FFD92F", "skyblue1", "#E78AC3" ,"lightcyan2", "#66C2A5" ,"coral" ,"cornflowerblue", "lightyellow2")
# pdf(paste0("mouse cell type proportion in spatial domains gCPCA.pdf"),width=6,height=6)
for(celltypeid in 1:dim(metatable)[2]){
  dat <- data.frame("Regions"=regionname,"Proportion"=metatable_proportion[,celltypeid])
  # dat$Regions <- factor(dat$Regions, 
  #                       levels=paste(rep("Cluster",nrow(metatable_proportion)),c(1:nrow(metatable_proportion))),order=T)
  dat$Regions <- factor(dat$Regions, 
                        levels=c("Normal lung","Fibrotic/Scar tissue","Adjacent normal","Bronchus/Alveoli",
                                 "Adenoma","Membrane"),order=T)
  pl2 <- ggplot(dat, aes(x=Regions, y=Proportion,fill=Proportion)) +
    geom_bar(stat="identity", position=position_dodge(),color="black",width=0.8)+
    #scale_fill_brewer(palette="Paired")+
    #scale_fill_manual(values = cbp_spatialpca)+
    scale_fill_distiller(palette = "Greens")+
    labs(title=paste0("gCPCA ",colnames(metatable)[celltypeid]),x="Spatial domains", y = "Proprotions")+
    theme(legend.position="right") +
    theme_classic()+
    #scale_fill_manual(values = method_color)+
    theme(axis.text.x = element_text(angle = 60,  hjust=1))+
    theme(plot.title = element_text(size = 22),
          text = element_text(size = 22),
          #axis.title = element_text(face="bold"),
          #axis.text.x=element_text(size = 20) ,
          legend.position = "right")
  
  print(pl2)
}


# perform trajectory analysis (tissue wise)
library(slingshot)
sim <- SingleCellExperiment(assays = STcounts)
reducedDims(sim) <- SimpleList(DRM = gCPCA[["fgComponents"]])
colData(sim)$clusterlabel <- factor(abnormalTissue@meta.data[["gcpca_clusters"]])    
sim  <-slingshot(sim, clusterLabels = 'clusterlabel', reducedDim = 'DRM', start.clus="0") 
# in this data we set normal epithelial region as start cluster 

summary(sim@colData@listData)
pseudotime_traj1 = sim@colData@listData$slingPseudotime_1 # in this data only one trajectory was inferred
gridnum = 10
# color_in = c("#9C9EDE" ,"#5CB85C" ,"#E377C2", "#4DBBD5" ,"#FED439" ,"#FF9896", "#FFDC91","black")
color_in = c("#E1BD6D","deepskyblue1","#7A3A9A","#ED0000FF","#0B775E","#ff00ff","black")
p_traj1 = plotTrajectory(pseudotime_traj1, flocation,abnormalTissue@meta.data[["gcpca_clusters"]],gridnum,
                         color_in,pointsize=1.5 ,arrowlength=0.2,arrowsize=1,textsize=15 )
p_traj1$Arrowoverlay1
p_traj1$Pseudotime

# cell-ell communication
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

sub_abnormalTissue <- subset(abnormalTissue, second_type != "NA")
data.input = GetAssayData(sub_abnormalTissue, slot = "data", assay = "SCT") # normalized data matrix
# meta = data.frame(labels = as.numeric(Idents(abnormalTissue)), row.names = colnames(abnormalTissue))
meta = data.frame(labels = sub_abnormalTissue@meta.data[["second_type"]], row.names = colnames(sub_abnormalTissue))
# meta$labels <- as.factor(meta$labels)
# meta$labels[is.na(meta$labels)] <- "Endothelial cells"
# meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))
# meta$labels <- factor(meta$labels, levels = c("Endothelial cells","Epithelial cells","Fibroblasts",
#                                               "Macrophages","cDC","Prolif_Macro","B cells","T cells",
#                                               "Proliferating T cells","pDC","Neutrophils"))
meta$labels <- factor(meta$labels, levels = unique(meta$labels))
scale.factors = jsonlite::fromJSON(txt = file.path("X:/maminu/Processed 129S4 Urethane model/53433/outs/spatial", 
                                                   'scalefactors_json.json'))
scale.factors = list(spot.diameter = 65, spot = scale.factors$spot_diameter_fullres, # these two information are required
                     fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef, 
                     lowres = scale.factors$tissue_lowres_scalef # these three information are not required
)

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = GetTissueCoordinates(sub_abnormalTissue), 
                           scale.factors = scale.factors)

CellChatDB <- CellChatDB.mouse
cellchat@DB <- CellChatDB

cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 1)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- projectData(cellchat, PPI.mouse)

cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
                              distance.use = TRUE, interaction.length = 200, scale.distance = 1)

cellchat <- filterCommunication(cellchat, min.cells = 0)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(meta$labels))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions",
                 color.use = c("plum1","tomato","#762A83","deepskyblue1","#C4961A","#ff00ff","#DC0000FF","#4E84C4",
                                      "chartreuse3","#D16103","#58593FFF"))
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength",
                 color.use = c("plum1","tomato","#762A83","deepskyblue1","#C4961A","#ff00ff","#DC0000FF","#4E84C4",
                                      "chartreuse3","#D16103","#58593FFF"))

mat1 <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat1)) {
  mat11 <- matrix(0, nrow = nrow(mat1), ncol = ncol(mat1), dimnames = dimnames(mat1))
  mat11[i, ] <- mat1[i, ]
  netVisual_circle(mat11, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat1), title.name = rownames(mat1)[i],
                   color.use = c("plum1","tomato","#762A83","deepskyblue1","#C4961A","#ff00ff","#DC0000FF","#4E84C4",
                                        "chartreuse3","#D16103","#58593FFF"))
}

cellchat@netP$pathways

pathways.show <- c("WNT") 
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle",
                    color.use = c("plum1","tomato","#762A83","deepskyblue1","#C4961A","#ff00ff","#DC0000FF","#4E84C4",
                                         "chartreuse3","#D16103","#58593FFF"))

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord", vertex.label.cex = 1.5,
                    color.use = c("plum1","tomato","#762A83","deepskyblue1","#C4961A","#ff00ff","#DC0000FF","#4E84C4",
                                         "chartreuse3","#D16103","#58593FFF"))

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds",
                  color.use = c("plum1","tomato","#762A83","deepskyblue1","#C4961A","#ff00ff","#DC0000FF","#4E84C4",
                                       "chartreuse3","#D16103","#58593FFF"))

# Spatial plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", edge.width.max = 2, 
                    vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5,
                    color.use = c("plum1","tomato","#762A83","deepskyblue1","#C4961A","#ff00ff","#DC0000FF","#4E84C4",
                                         "chartreuse3","#D16103","#58593FFF"))

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10,
                                  color.use = c("plum1","tomato","#762A83","deepskyblue1","#C4961A","#ff00ff","#DC0000FF","#4E84C4",
                                                       "chartreuse3","#D16103","#58593FFF"))



############### Apply transformation on unseen data  #####################
# load second foreground data 
abnormalTissue2 <- Load10X_Spatial(
  data.dir = "X:/maminu/Processed 129S4 Urethane model/26932/outs", 
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", # specify name of the initial assay
  filter.matrix = TRUE, 
  to.upper = FALSE
)

abnormalTissue2 <- SCTransform(abnormalTissue2, assay = "Spatial", verbose = FALSE)
SpatialFeaturePlot(abnormalTissue2, features = "nCount_Spatial") + theme(legend.position = "right")

fdata2 <- abnormalTissue2@assays[["SCT"]]@scale.data
# flocation2 <- abnormalTissue2@images[["slice1"]]@coordinates[,c(2,3)] 
flocation2 <- GetTissueCoordinates(abnormalTissue2)

fgComponents2 <- t(fdata2) %*% gCPCA[["projMatrix"]]

fgraphCPCA2 <- CreateDimReducObject(
  embeddings = fgComponents2,
  loadings = gCPCA[["projMatrix"]],
  stdev = numeric(),
  key = "gCPCA_",
  assay = "SCT"
)

rownames(fgraphCPCA2@feature.loadings) <- rownames(fdata2)
abnormalTissue2@reductions[["graphCPCA"]] <- fgraphCPCA2

abnormalTissue2 <- RunUMAP(abnormalTissue2, reduction = "graphCPCA", dims = 1:10, reduction.name = "gCPCA_UMAP")
abnormalTissue2 <- FindNeighbors(abnormalTissue2, reduction = "graphCPCA", dims = 1:10)
abnormalTissue2 <- FindClusters(abnormalTissue2, verbose = TRUE, resolution = 0.44)
abnormalTissue2@meta.data[["gcpca_clusters"]]  <- abnormalTissue2@meta.data[["SCT_snn_res.0.44"]]

# DimPlot(abnormalTissue2, reduction = "gCPCA_UMAP", label = FALSE, cols = c("#E1BD6D","#ED0000FF","#44B05B"))
SpatialDimPlot(abnormalTissue2, label = FALSE, label.size = 3, group.by = "gcpca_clusters", pt.size.factor = 2) +  
  scale_fill_manual(values=c("#E1BD6D","#ED0000FF","#44B05B"))


# set up query with the RCTD function SpatialRNA
STcounts2 <- abnormalTissue2@assays[["Spatial"]]@counts
coords2 <- GetTissueCoordinates(abnormalTissue2)
colnames(coords2) <- c("x", "y")
coords2[is.na(colnames(coords2))] <- NULL
query2 <- SpatialRNA(coords2, STcounts2, colSums(STcounts2))

RCTD2 <- create.RCTD(query2, reference, max_cores = 1)
RCTD2 <- run.RCTD(RCTD2, doublet_mode = "doublet")
abnormalTissue2 <- AddMetaData(abnormalTissue2, metadata = RCTD2@results$results_df)

SpatialDimPlot(abnormalTissue2, group.by = "first_type")
SpatialDimPlot(abnormalTissue2, group.by = "second_type") +  scale_fill_manual(values=col2)
                                    
# perform trajectory analysis (tissue wise)
sim2 <- SingleCellExperiment(assays = STcounts2)
reducedDims(sim2) <- SimpleList(DRM = fgComponents2)
colData(sim2)$clusterlabel <- factor(abnormalTissue2@meta.data[["gcpca_clusters"]])    
sim2  <-slingshot(sim2, clusterLabels = 'clusterlabel', reducedDim = 'DRM', start.clus="0") 
# in this data we set normal epithelial region as start cluster 

summary(sim2@colData@listData)
pseudotime_traj2 = sim2@colData@listData$slingPseudotime_1 # in this data only one trajectory was inferred
# gridnum = 10
# # color_in = c("#9C9EDE" ,"#5CB85C" ,"#E377C2", "#4DBBD5" ,"#FED439" ,"#FF9896", "#FFDC91","black")
# color_in = c("#E1BD6D","deepskyblue1","#7A3A9A","#ED0000FF","#0B775E","#ff00ff","black")
p_traj2 = plotTrajectory(pseudotime_traj2, coords2,abnormalTissue2@meta.data[["gcpca_clusters"]],gridnum,
                         color_in,pointsize=1.5 ,arrowlength=0.2,arrowsize=1,textsize=15 )
p_traj2$Arrowoverlay1
p_traj2$Pseudotime


# cell-ell communication
# sub_abnormalTissue <- subset(abnormalTissue, second_type != "NA")
data.input2 = GetAssayData(abnormalTissue2, slot = "data", assay = "SCT") # normalized data matrix
# meta = data.frame(labels = as.numeric(Idents(abnormalTissue)), row.names = colnames(abnormalTissue))
meta2 = data.frame(labels = abnormalTissue2@meta.data[["second_type"]], row.names = colnames(abnormalTissue2))
# meta2$labels <- as.factor(meta2$labels)
# meta$labels[is.na(meta$labels)] <- "Endothelial cells"
# meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))
# meta2$labels <- factor(meta2$labels, levels = c("Endothelial cells","Epithelial cells","Fibroblasts",
#                                                 "Macrophages","cDC","B cells"))
meta2$labels <- factor(meta2$labels, levels = unique(meta2$labels))
scale.factors2 = jsonlite::fromJSON(txt = file.path("X:/maminu/Processed 129S4 Urethane model/26932/outs/spatial", 
                                                    'scalefactors_json.json'))
scale.factors2 = list(spot.diameter = 65, spot = scale.factors2$spot_diameter_fullres, # these two information are required
                      fiducial = scale.factors2$fiducial_diameter_fullres, hires = scale.factors2$tissue_hires_scalef, 
                      lowres = scale.factors2$tissue_lowres_scalef # these three information are not required
)

cellchat2 <- createCellChat(object = data.input2, meta = meta2, group.by = "labels",
                            datatype = "spatial", coordinates = GetTissueCoordinates(abnormalTissue2), 
                            scale.factors = scale.factors2)

# CellChatDB <- CellChatDB.mouse
cellchat2@DB <- CellChatDB

cellchat2 <- subsetData(cellchat2) # This step is necessary even if using the whole database
future::plan("multisession", workers = 1)

cellchat2 <- identifyOverExpressedGenes(cellchat2)
cellchat2 <- identifyOverExpressedInteractions(cellchat2)

cellchat2 <- projectData(cellchat2, PPI.mouse)

cellchat2 <- computeCommunProb(cellchat2, type = "truncatedMean", trim = 0.1, 
                               distance.use = TRUE, interaction.length = 200, scale.distance = 1)

cellchat2 <- filterCommunication(cellchat2, min.cells = 0)
cellchat2 <- computeCommunProbPathway(cellchat2)
cellchat2 <- aggregateNet(cellchat2)

groupSize2 <- as.numeric(table(meta2$labels))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat2@net$count, vertex.weight = rowSums(cellchat2@net$count), weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions", 
                 color.use = c("plum1","#762A83","tomato","deepskyblue1","#DC0000FF","#C4961A"))
netVisual_circle(cellchat2@net$weight, vertex.weight = rowSums(cellchat2@net$weight), weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength",
                 color.use = c("plum1","#762A83","tomato","deepskyblue1","#DC0000FF","#C4961A"))

mat2 <- cellchat2@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat2)) {
  mat22 <- matrix(0, nrow = nrow(mat2), ncol = ncol(mat2), dimnames = dimnames(mat2))
  mat22[i, ] <- mat2[i, ]
  netVisual_circle(mat22, vertex.weight = groupSize2, weight.scale = T, edge.weight.max = max(mat2), title.name = rownames(mat2)[i],
                   color.use = c("plum1","#762A83","tomato","deepskyblue1","#DC0000FF","#C4961A"))
}

cellchat2@netP$pathways

pathways.show <- c("WNT") 
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat2, signaling = pathways.show, layout = "circle",
                    color.use = c("plum1","#762A83","tomato","deepskyblue1","#DC0000FF","#C4961A"))

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat2, signaling = pathways.show, layout = "chord", vertex.label.cex = 0.8,
                    color.use = c("plum1","#762A83","tomato","deepskyblue1","#DC0000FF","#C4961A"))

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat2, signaling = pathways.show, color.heatmap = "Reds",
                  color.use = c("plum1","#762A83","tomato","deepskyblue1","#DC0000FF","#C4961A"))

# Spatial plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat2, signaling = pathways.show, layout = "spatial", edge.width.max = 2, 
                    vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5,
                    color.use = c("plum1","#762A83","tomato","deepskyblue1","#DC0000FF","#C4961A"))

# Compute the network centrality scores
cellchat2 <- netAnalysis_computeCentrality(cellchat2, slot.name = "netP") 
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat2, signaling = pathways.show, width = 8, height = 2.5, font.size = 10,
                                  color.use = c("plum1","#762A83","tomato","deepskyblue1","#DC0000FF","#C4961A"))


# load third foreground data 
abnormalTissue3 <- Load10X_Spatial(
  data.dir = "X:/maminu/Processed 129S4 Urethane model/53431/outs", 
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", # specify name of the initial assay
  filter.matrix = TRUE, 
  to.upper = FALSE
)

abnormalTissue3 <- SCTransform(abnormalTissue3, assay = "Spatial", verbose = FALSE)
SpatialFeaturePlot(abnormalTissue3, features = "nCount_Spatial") + theme(legend.position = "right")

fdata3 <- abnormalTissue3@assays[["SCT"]]@scale.data
# flocation2 <- abnormalTissue2@images[["slice1"]]@coordinates[,c(2,3)] 
flocation3 <- GetTissueCoordinates(abnormalTissue3)

fgComponents3 <- t(fdata3) %*% gCPCA[["projMatrix"]]

fgraphCPCA3 <- CreateDimReducObject(
  embeddings = fgComponents3,
  loadings = gCPCA[["projMatrix"]],
  stdev = numeric(),
  key = "gCPCA_",
  assay = "SCT"
)

rownames(fgraphCPCA3@feature.loadings) <- rownames(fdata3)
abnormalTissue3@reductions[["graphCPCA"]] <- fgraphCPCA3

abnormalTissue3 <- RunUMAP(abnormalTissue3, reduction = "graphCPCA", dims = 1:10, reduction.name = "gCPCA_UMAP")
abnormalTissue3 <- FindNeighbors(abnormalTissue3, reduction = "graphCPCA", dims = 1:20)
abnormalTissue3 <- FindClusters(abnormalTissue3, verbose = TRUE, resolution = 0.1)
abnormalTissue3@meta.data[["gcpca_clusters"]]  <- abnormalTissue3@meta.data[["SCT_snn_res.0.1"]]

# DimPlot(abnormalTissue3, reduction = "gCPCA_UMAP", label = FALSE, cols = c("#E1BD6D","#ED0000FF","#0B775E"))
SpatialDimPlot(abnormalTissue3, label = FALSE, label.size = 3, group.by = "gcpca_clusters", pt.size.factor = 2) +  
  scale_fill_manual(values=c("#E1BD6D","#ED0000FF","#0B775E"))

# set up query with the RCTD function SpatialRNA
STcounts3 <- abnormalTissue3@assays[["Spatial"]]@counts
coords3 <- GetTissueCoordinates(abnormalTissue3)
colnames(coords3) <- c("x", "y")
coords3[is.na(colnames(coords3))] <- NULL
query3 <- SpatialRNA(coords3, STcounts3, colSums(STcounts3))

RCTD3 <- create.RCTD(query3, reference, max_cores = 1)
RCTD3 <- run.RCTD(RCTD3, doublet_mode = "doublet")
abnormalTissue3 <- AddMetaData(abnormalTissue3, metadata = RCTD3@results$results_df)

SpatialDimPlot(abnormalTissue3, group.by = "first_type")
SpatialDimPlot(abnormalTissue3, group.by = "second_type") +  scale_fill_manual(values=col2)
                                    

# perform trajectory analysis (tissue wise)
sim3 <- SingleCellExperiment(assays = STcounts3)
reducedDims(sim3) <- SimpleList(DRM = fgComponents3)
colData(sim3)$clusterlabel <- factor(abnormalTissue3@meta.data[["gcpca_clusters"]])    
sim3  <-slingshot(sim3, clusterLabels = 'clusterlabel', reducedDim = 'DRM', start.clus="0") 
# in this data we set normal epithelial region as start cluster 

summary(sim3@colData@listData)
pseudotime_traj3 = sim3@colData@listData$slingPseudotime_1 # in this data only one trajectory was inferred
# gridnum = 10
# # color_in = c("#9C9EDE" ,"#5CB85C" ,"#E377C2", "#4DBBD5" ,"#FED439" ,"#FF9896", "#FFDC91","black")
# color_in = c("#E1BD6D","deepskyblue1","#7A3A9A","#ED0000FF","#0B775E","#ff00ff","black")
p_traj3 = plotTrajectory(pseudotime_traj3, coords3,abnormalTissue3@meta.data[["gcpca_clusters"]],gridnum,
                         color_in,pointsize=1.5 ,arrowlength=0.2,arrowsize=1,textsize=15 )
p_traj3$Arrowoverlay1
p_traj3$Pseudotime


# cell-ell communication
# sub_abnormalTissue3 <- subset(abnormalTissue3, second_type != "NA")
data.input3 = GetAssayData(abnormalTissue3, slot = "data", assay = "SCT") # normalized data matrix
# meta = data.frame(labels = as.numeric(Idents(abnormalTissue)), row.names = colnames(abnormalTissue))
meta3 = data.frame(labels = abnormalTissue3@meta.data[["second_type"]], row.names = colnames(abnormalTissue3))
# meta3$labels <- as.factor(meta3$labels)
# meta$labels[is.na(meta$labels)] <- "Endothelial cells"
# meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))
# meta3$labels <- factor(meta3$labels, levels = c("Endothelial cells","Epithelial cells","Fibroblasts",
#                                                 "Macrophages","cDC","Prolif_Macro","B cells","T cells",
#                                                 "Proliferating T cells","pDC","Neutrophils","Plasma cells",
#                                                 "Monocytes","NK cells"))
meta3$labels <- factor(meta3$labels, levels = unique(meta3$labels))
scale.factors3 = jsonlite::fromJSON(txt = file.path("X:/maminu/Processed 129S4 Urethane model/53431/outs/spatial", 
                                                    'scalefactors_json.json'))
scale.factors3 = list(spot.diameter = 65, spot = scale.factors3$spot_diameter_fullres, # these two information are required
                      fiducial = scale.factors3$fiducial_diameter_fullres, hires = scale.factors3$tissue_hires_scalef, 
                      lowres = scale.factors3$tissue_lowres_scalef # these three information are not required
)

cellchat3 <- createCellChat(object = data.input3, meta = meta3, group.by = "labels",
                            datatype = "spatial", coordinates = GetTissueCoordinates(abnormalTissue3), 
                            scale.factors = scale.factors3)

# CellChatDB <- CellChatDB.mouse
cellchat3@DB <- CellChatDB

cellchat3 <- subsetData(cellchat3) # This step is necessary even if using the whole database
future::plan("multisession", workers = 1)

cellchat3 <- identifyOverExpressedGenes(cellchat3)
cellchat3 <- identifyOverExpressedInteractions(cellchat3)

cellchat3 <- projectData(cellchat3, PPI.mouse)

cellchat3 <- computeCommunProb(cellchat3, type = "truncatedMean", trim = 0.1, 
                               distance.use = TRUE, interaction.length = 200, scale.distance = 1)

cellchat3 <- filterCommunication(cellchat3, min.cells = 0)
cellchat3 <- computeCommunProbPathway(cellchat3)
cellchat3 <- aggregateNet(cellchat3)

groupSize3 <- as.numeric(table(meta3$labels))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat3@net$count, vertex.weight = rowSums(cellchat3@net$count), weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions",
                 color.use = c("tomato","#762A83","plum1","#DC0000FF","deepskyblue1","chartreuse3","#C4961A","#4E84C4",
                                       "#ff00ff"))
netVisual_circle(cellchat3@net$weight, vertex.weight = rowSums(cellchat3@net$weight), weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength",
                 color.use = c("tomato","#762A83","plum1","#DC0000FF","deepskyblue1","chartreuse3","#C4961A","#4E84C4",
                                       "#ff00ff"))
                     
mat3 <- cellchat3@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat3)) {
  mat33 <- matrix(0, nrow = nrow(mat3), ncol = ncol(mat3), dimnames = dimnames(mat3))
  mat33[i, ] <- mat3[i, ]
  netVisual_circle(mat33, vertex.weight = groupSize3, weight.scale = T, edge.weight.max = max(mat3), title.name = rownames(mat3)[i],
                   color.use = c("tomato","#762A83","plum1","#DC0000FF","deepskyblue1","chartreuse3","#C4961A","#4E84C4",
                                         "#ff00ff"))
}

cellchat3@netP$pathways

pathways.show <- c("WNT") 
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat3, signaling = pathways.show, layout = "circle",
                    color.use = c("tomato","#762A83","plum1","#DC0000FF","deepskyblue1","chartreuse3","#C4961A","#4E84C4",
                                          "#ff00ff"))

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat3, signaling = pathways.show, layout = "chord", vertex.label.cex = 0.8,
                    color.use = c("tomato","#762A83","plum1","#DC0000FF","deepskyblue1","chartreuse3","#C4961A","#4E84C4",
                                          "#ff00ff"))

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat3, signaling = pathways.show, color.heatmap = "Reds",
                  color.use = c("tomato","#762A83","plum1","#DC0000FF","deepskyblue1","chartreuse3","#C4961A","#4E84C4",
                                        "#ff00ff"))

# Spatial plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat3, signaling = pathways.show, layout = "spatial", edge.width.max = 2, 
                    vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5,
                    color.use = c("tomato","#762A83","plum1","#DC0000FF","deepskyblue1","chartreuse3","#C4961A","#4E84C4",
                                          "#ff00ff"))

# Compute the network centrality scores
cellchat3 <- netAnalysis_computeCentrality(cellchat3, slot.name = "netP") 
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat3, signaling = pathways.show, width = 8, height = 2.5, font.size = 10,
                                  color.use = c("tomato","#762A83","plum1","#DC0000FF","deepskyblue1","chartreuse3","#C4961A","#4E84C4",
                                                        "#ff00ff"))




# load fourth foreground data 
abnormalTissue4 <- Load10X_Spatial(
  data.dir = "X:/maminu/Processed 129S4 Urethane model/53432/outs", 
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", # specify name of the initial assay
  filter.matrix = TRUE, 
  to.upper = FALSE
)

abnormalTissue4 <- SCTransform(abnormalTissue4, assay = "Spatial", verbose = FALSE)
SpatialFeaturePlot(abnormalTissue4, features = "nCount_Spatial") + theme(legend.position = "right")

fdata4 <- abnormalTissue4@assays[["SCT"]]@scale.data
# flocation2 <- abnormalTissue2@images[["slice1"]]@coordinates[,c(2,3)] 
flocation4 <- GetTissueCoordinates(abnormalTissue4)

fgComponents4 <- t(fdata4) %*% gCPCA[["projMatrix"]]

fgraphCPCA4 <- CreateDimReducObject(
  embeddings = fgComponents4,
  loadings = gCPCA[["projMatrix"]],
  stdev = numeric(),
  key = "gCPCA_",
  assay = "SCT"
)

rownames(fgraphCPCA4@feature.loadings) <- rownames(fdata4)
abnormalTissue4@reductions[["graphCPCA"]] <- fgraphCPCA4

abnormalTissue4 <- RunUMAP(abnormalTissue4, reduction = "graphCPCA", dims = 1:10, reduction.name = "gCPCA_UMAP")
abnormalTissue4 <- FindNeighbors(abnormalTissue4, reduction = "graphCPCA", dims = 1:20)
abnormalTissue4 <- FindClusters(abnormalTissue4, verbose = TRUE, resolution = 0.4)
abnormalTissue4@meta.data[["gcpca_clusters"]]  <- abnormalTissue4@meta.data[["SCT_snn_res.0.4"]]

# DimPlot(abnormalTissue4, reduction = "gCPCA_UMAP", label = FALSE, cols = c("#E1BD6D","#ED0000FF","#0B775E"))
SpatialDimPlot(abnormalTissue4, label = FALSE, label.size = 3, group.by = "gcpca_clusters", pt.size.factor = 2) +  
  scale_fill_manual(values=c("#E1BD6D","#ED0000FF","#0B775E"))


# set up query with the RCTD function SpatialRNA
STcounts4 <- abnormalTissue4@assays[["Spatial"]]@counts
coords4 <- GetTissueCoordinates(abnormalTissue4)
colnames(coords4) <- c("x", "y")
coords4[is.na(colnames(coords4))] <- NULL
query4 <- SpatialRNA(coords4, STcounts4, colSums(STcounts4))

RCTD4 <- create.RCTD(query4, reference, max_cores = 1)
RCTD4 <- run.RCTD(RCTD4, doublet_mode = "doublet")
abnormalTissue4 <- AddMetaData(abnormalTissue4, metadata = RCTD4@results$results_df)

SpatialDimPlot(abnormalTissue4, group.by = "first_type")
SpatialDimPlot(abnormalTissue4, group.by = "second_type") +  scale_fill_manual(values=col2)
                                    
# perform trajectory analysis (tissue wise)
sim4 <- SingleCellExperiment(assays = STcounts4)
reducedDims(sim4) <- SimpleList(DRM = fgComponents4)
colData(sim4)$clusterlabel <- factor(abnormalTissue4@meta.data[["gcpca_clusters"]])    
sim4  <-slingshot(sim4, clusterLabels = 'clusterlabel', reducedDim = 'DRM', start.clus="0") 
# in this data we set normal epithelial region as start cluster 

summary(sim4@colData@listData)
pseudotime_traj4 = sim4@colData@listData$slingPseudotime_1 # in this data only one trajectory was inferred
# gridnum = 10
# # color_in = c("#9C9EDE" ,"#5CB85C" ,"#E377C2", "#4DBBD5" ,"#FED439" ,"#FF9896", "#FFDC91","black")
# color_in = c("#E1BD6D","deepskyblue1","#7A3A9A","#ED0000FF","#0B775E","#ff00ff","black")
p_traj4 = plotTrajectory(pseudotime_traj4, coords4,abnormalTissue4@meta.data[["gcpca_clusters"]],gridnum,
                         color_in,pointsize=1.5 ,arrowlength=0.2,arrowsize=1,textsize=15 )
p_traj4$Arrowoverlay1
p_traj4$Pseudotime


# cell-ell communication
# sub_abnormalTissue3 <- subset(abnormalTissue3, second_type != "NA")
data.input4 = GetAssayData(abnormalTissue4, slot = "data", assay = "SCT") # normalized data matrix
# meta = data.frame(labels = as.numeric(Idents(abnormalTissue)), row.names = colnames(abnormalTissue))
meta4 = data.frame(labels = abnormalTissue4@meta.data[["second_type"]], row.names = colnames(abnormalTissue4))
meta4$labels <- as.factor(meta4$labels)
# meta$labels[is.na(meta$labels)] <- "Endothelial cells"
# meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))
# meta4$labels <- factor(meta4$labels, levels = c("Endothelial cells","Epithelial cells","Fibroblasts",
#                                                 "Macrophages","cDC","Prolif_Macro","B cells","T cells",
#                                                 "Proliferating T cells","pDC","Neutrophils","Plasma cells",
#                                                 "Monocytes","NK cells"))
meta4$labels <- factor(meta4$labels, levels = unique(meta4$labels))
scale.factors4 = jsonlite::fromJSON(txt = file.path("X:/maminu/Processed 129S4 Urethane model/53432/outs/spatial", 
                                                    'scalefactors_json.json'))
scale.factors4 = list(spot.diameter = 65, spot = scale.factors4$spot_diameter_fullres, # these two information are required
                      fiducial = scale.factors4$fiducial_diameter_fullres, hires = scale.factors4$tissue_hires_scalef, 
                      lowres = scale.factors4$tissue_lowres_scalef # these three information are not required
)

cellchat4 <- createCellChat(object = data.input4, meta = meta4, group.by = "labels",
                            datatype = "spatial", coordinates = GetTissueCoordinates(abnormalTissue4), 
                            scale.factors = scale.factors4)

# CellChatDB <- CellChatDB.mouse
cellchat4@DB <- CellChatDB

cellchat4 <- subsetData(cellchat4) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4)

cellchat4 <- identifyOverExpressedGenes(cellchat4)
cellchat4 <- identifyOverExpressedInteractions(cellchat4)

cellchat4 <- projectData(cellchat4, PPI.mouse)

cellchat4 <- computeCommunProb(cellchat4, type = "truncatedMean", trim = 0.1, 
                               distance.use = TRUE, interaction.length = 200, scale.distance = 1)

cellchat4 <- filterCommunication(cellchat4, min.cells = 0)
cellchat4 <- computeCommunProbPathway(cellchat4)
cellchat4 <- aggregateNet(cellchat4)

groupSize4 <- as.numeric(table(meta4$labels))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat4@net$count, vertex.weight = rowSums(cellchat4@net$count), weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions",
                 color.use = c("plum1","tomato","#762A83","chartreuse3","deepskyblue1","#D16103","#4E84C4","#C4961A",
                                      "#DC0000FF","#ff00ff","lightblue1"))
netVisual_circle(cellchat4@net$weight, vertex.weight = rowSums(cellchat4@net$weight), weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength",
                 color.use = c("plum1","tomato","#762A83","chartreuse3","deepskyblue1","#D16103","#4E84C4","#C4961A",
                                      "#DC0000FF","#ff00ff","lightblue1"))

mat4 <- cellchat4@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat4)) {
  mat44 <- matrix(0, nrow = nrow(mat4), ncol = ncol(mat4), dimnames = dimnames(mat4))
  mat44[i, ] <- mat4[i, ]
  netVisual_circle(mat44, vertex.weight = groupSize4, weight.scale = T, edge.weight.max = max(mat4), title.name = rownames(mat4)[i],
                   color.use = c("plum1","tomato","#762A83","chartreuse3","deepskyblue1","#D16103","#4E84C4","#C4961A",
                                        "#DC0000FF","#ff00ff","lightblue1"))
}

cellchat4@netP$pathways

pathways.show <- c("WNT") 
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat4, signaling = pathways.show, layout = "circle",
                    color.use = c("plum1","tomato","#762A83","chartreuse3","deepskyblue1","#D16103","#4E84C4","#C4961A",
                                         "#DC0000FF","#ff00ff","lightblue1"))

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat4, signaling = pathways.show, layout = "chord", vertex.label.cex = 0.8,
                    color.use = c("plum1","tomato","#762A83","chartreuse3","deepskyblue1","#D16103","#4E84C4","#C4961A",
                                         "#DC0000FF","#ff00ff","lightblue1"))

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat4, signaling = pathways.show, color.heatmap = "Reds",
                  color.use = c("plum1","tomato","#762A83","chartreuse3","deepskyblue1","#D16103","#4E84C4","#C4961A",
                                       "#DC0000FF","#ff00ff","lightblue1"))

# Spatial plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat4, signaling = pathways.show, layout = "spatial", edge.width.max = 2, 
                    vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5,
                    color.use = c("plum1","tomato","#762A83","chartreuse3","deepskyblue1","#D16103","#4E84C4","#C4961A",
                                         "#DC0000FF","#ff00ff","lightblue1"))

# Compute the network centrality scores
cellchat4 <- netAnalysis_computeCentrality(cellchat4, slot.name = "netP") 
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat4, signaling = pathways.show, width = 8, height = 2.5, font.size = 10,
                                  color.use = c("plum1","tomato","#762A83","chartreuse3","deepskyblue1","#D16103","#4E84C4","#C4961A",
                                                       "#DC0000FF","#ff00ff","lightblue1"))




# load fifth foreground data 
abnormalTissue5 <- Load10X_Spatial(
  data.dir = "X:/maminu/Processed 129S4 Urethane model/26933/outs", 
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", # specify name of the initial assay
  filter.matrix = TRUE, 
  to.upper = FALSE
)

abnormalTissue5 <- SCTransform(abnormalTissue5, assay = "Spatial", verbose = FALSE)
SpatialFeaturePlot(abnormalTissue5, features = "nCount_Spatial") + theme(legend.position = "right")

fdata5 <- abnormalTissue5@assays[["SCT"]]@scale.data
# flocation2 <- abnormalTissue2@images[["slice1"]]@coordinates[,c(2,3)] 
flocation5 <- GetTissueCoordinates(abnormalTissue5)

fgComponents5 <- t(fdata5) %*% gCPCA[["projMatrix"]]

fgraphCPCA5 <- CreateDimReducObject(
  embeddings = fgComponents5,
  loadings = gCPCA[["projMatrix"]],
  stdev = numeric(),
  key = "gCPCA_",
  assay = "SCT"
)

rownames(fgraphCPCA5@feature.loadings) <- rownames(fdata5)
abnormalTissue5@reductions[["graphCPCA"]] <- fgraphCPCA5

abnormalTissue5 <- RunUMAP(abnormalTissue5, reduction = "graphCPCA", dims = 1:10, reduction.name = "gCPCA_UMAP")
abnormalTissue5 <- FindNeighbors(abnormalTissue5, reduction = "graphCPCA", dims = 1:50)
abnormalTissue5 <- FindClusters(abnormalTissue5, verbose = TRUE, resolution = 0.65)
abnormalTissue5@meta.data[["gcpca_clusters"]]  <- abnormalTissue5@meta.data[["SCT_snn_res.0.65"]]

# DimPlot(abnormalTissue5, reduction = "gCPCA_UMAP", label = FALSE, cols = c("#E1BD6D","#ED0000FF","#44B05B","#0B775E"))
SpatialDimPlot(abnormalTissue5, label = FALSE, label.size = 3, group.by = "gcpca_clusters", pt.size.factor = 2.2) +  
  scale_fill_manual(values=c("#E1BD6D","#9A8822","#ED0000FF","#0B775E","#ff00ff","#7B556C"))


# set up query with the RCTD function SpatialRNA
STcounts5 <- abnormalTissue5@assays[["Spatial"]]@counts
coords5 <- GetTissueCoordinates(abnormalTissue5)
colnames(coords5) <- c("x", "y")
coords5[is.na(colnames(coords5))] <- NULL
query5 <- SpatialRNA(coords5, STcounts5, colSums(STcounts5))

RCTD5 <- create.RCTD(query5, reference, max_cores = 1)
RCTD5 <- run.RCTD(RCTD5, doublet_mode = "doublet")
abnormalTissue5 <- AddMetaData(abnormalTissue5, metadata = RCTD5@results$results_df)

SpatialDimPlot(abnormalTissue5, group.by = "first_type")
SpatialDimPlot(abnormalTissue5, group.by = "second_type") +  scale_fill_manual(values=col2)
                                    
# perform trajectory analysis (tissue wise)
sim5 <- SingleCellExperiment(assays = STcounts5)
reducedDims(sim5) <- SimpleList(DRM = fgComponents5)
colData(sim5)$clusterlabel <- factor(abnormalTissue5@meta.data[["gcpca_clusters"]])    
sim5  <-slingshot(sim5, clusterLabels = 'clusterlabel', reducedDim = 'DRM', start.clus="0") 
# in this data we set normal epithelial region as start cluster 

summary(sim5@colData@listData)
pseudotime_traj5 = sim5@colData@listData$slingPseudotime_1 # in this data only one trajectory was inferred
# gridnum = 10
# # color_in = c("#9C9EDE" ,"#5CB85C" ,"#E377C2", "#4DBBD5" ,"#FED439" ,"#FF9896", "#FFDC91","black")
# color_in = c("#E1BD6D","deepskyblue1","#7A3A9A","#ED0000FF","#0B775E","#ff00ff","black")
p_traj5 = plotTrajectory(pseudotime_traj5, coords5,abnormalTissue5@meta.data[["gcpca_clusters"]],gridnum,
                         color_in,pointsize=1.5 ,arrowlength=0.2,arrowsize=1,textsize=15 )
p_traj5$Arrowoverlay1
p_traj5$Pseudotime


# cell-ell communication
# sub_abnormalTissue3 <- subset(abnormalTissue3, second_type != "NA")
data.input5 = GetAssayData(abnormalTissue5, slot = "data", assay = "SCT") # normalized data matrix
# meta = data.frame(labels = as.numeric(Idents(abnormalTissue)), row.names = colnames(abnormalTissue))
meta5 = data.frame(labels = abnormalTissue5@meta.data[["second_type"]], row.names = colnames(abnormalTissue5))
meta5$labels <- as.factor(meta5$labels)
# meta$labels[is.na(meta$labels)] <- "Endothelial cells"
# meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))
# meta5$labels <- factor(meta5$labels, levels = c("Endothelial cells","Epithelial cells","Fibroblasts",
#                                                 "Macrophages","cDC","Prolif_Macro","B cells","T cells",
#                                                 "Proliferating T cells","pDC","Neutrophils","Plasma cells",
#                                                 "Monocytes","NK cells"))
meta5$labels <- factor(meta5$labels, levels = unique(meta5$labels))
scale.factors5 = jsonlite::fromJSON(txt = file.path("X:/maminu/Processed 129S4 Urethane model/26933/outs/spatial", 
                                                    'scalefactors_json.json'))
scale.factors5 = list(spot.diameter = 65, spot = scale.factors5$spot_diameter_fullres, # these two information are required
                      fiducial = scale.factors5$fiducial_diameter_fullres, hires = scale.factors5$tissue_hires_scalef, 
                      lowres = scale.factors5$tissue_lowres_scalef # these three information are not required
)

cellchat5 <- createCellChat(object = data.input5, meta = meta5, group.by = "labels",
                            datatype = "spatial", coordinates = GetTissueCoordinates(abnormalTissue5), 
                            scale.factors = scale.factors5)

# CellChatDB <- CellChatDB.mouse
cellchat5@DB <- CellChatDB

cellchat5 <- subsetData(cellchat5) # This step is necessary even if using the whole database
future::plan("multisession", workers = 1)

cellchat5 <- identifyOverExpressedGenes(cellchat5)
cellchat5 <- identifyOverExpressedInteractions(cellchat5)

cellchat5 <- projectData(cellchat5, PPI.mouse)

cellchat5 <- computeCommunProb(cellchat5, type = "truncatedMean", trim = 0.1, 
                               distance.use = TRUE, interaction.length = 200, scale.distance = 1)

cellchat5 <- filterCommunication(cellchat5, min.cells = 0)
cellchat5 <- computeCommunProbPathway(cellchat5)
cellchat5 <- aggregateNet(cellchat5)

groupSize5 <- as.numeric(table(meta5$labels))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat5@net$count, vertex.weight = rowSums(cellchat5@net$count), weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions",
                 color.use = c("plum1","#762A83","tomato","deepskyblue1","chartreuse3","#C4961A","#DC0000FF","#4E84C4",
                                      "#ff00ff"))
netVisual_circle(cellchat5@net$weight, vertex.weight = rowSums(cellchat5@net$weight), weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength",
                 color.use = c("plum1","#762A83","tomato","deepskyblue1","chartreuse3","#C4961A","#DC0000FF","#4E84C4",
                                      "#ff00ff"))

mat5 <- cellchat5@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat5)) {
  mat55 <- matrix(0, nrow = nrow(mat5), ncol = ncol(mat5), dimnames = dimnames(mat5))
  mat55[i, ] <- mat5[i, ]
  netVisual_circle(mat55, vertex.weight = groupSize5, weight.scale = T, edge.weight.max = max(mat5), title.name = rownames(mat5)[i],
                   color.use = c("plum1","#762A83","tomato","deepskyblue1","chartreuse3","#C4961A","#DC0000FF","#4E84C4",
                                        "#ff00ff"))
}

cellchat5@netP$pathways

pathways.show <- c("WNT") 
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat5, signaling = pathways.show, layout = "circle",
                    color.use = c("plum1","#762A83","tomato","deepskyblue1","chartreuse3","#C4961A","#DC0000FF","#4E84C4",
                                         "#ff00ff"))

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat5, signaling = pathways.show, layout = "chord", vertex.label.cex = 0.8,
                    color.use = c("plum1","#762A83","tomato","deepskyblue1","chartreuse3","#C4961A","#DC0000FF","#4E84C4",
                                         "#ff00ff"))

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat5, signaling = pathways.show, color.heatmap = "Reds",
                  color.use = c("plum1","#762A83","tomato","deepskyblue1","chartreuse3","#C4961A","#DC0000FF","#4E84C4",
                                       "#ff00ff"))

# Spatial plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat5, signaling = pathways.show, layout = "spatial", edge.width.max = 2, 
                    vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5,
                    color.use = c("plum1","#762A83","tomato","deepskyblue1","chartreuse3","#C4961A","#DC0000FF","#4E84C4",
                                         "#ff00ff"))

# Compute the network centrality scores
cellchat5 <- netAnalysis_computeCentrality(cellchat5, slot.name = "netP") 
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat5, signaling = pathways.show, width = 8, height = 2.5, font.size = 10,
                                  color.use = c("plum1","#762A83","tomato","deepskyblue1","chartreuse3","#C4961A","#DC0000FF","#4E84C4",
                                                       "#ff00ff"))





# load sixth foreground data 
abnormalTissue6 <- Load10X_Spatial(
  data.dir = "X:/maminu/Processed 129S4 Urethane model/26934/outs", 
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", # specify name of the initial assay
  filter.matrix = TRUE, 
  to.upper = FALSE
)

abnormalTissue6 <- SCTransform(abnormalTissue6, assay = "Spatial", verbose = FALSE)
SpatialFeaturePlot(abnormalTissue6, features = "nCount_Spatial") + theme(legend.position = "right")

fdata6 <- abnormalTissue6@assays[["SCT"]]@scale.data
# flocation2 <- abnormalTissue2@images[["slice1"]]@coordinates[,c(2,3)] 
flocation6 <- GetTissueCoordinates(abnormalTissue6)

fgComponents6 <- t(fdata6) %*% gCPCA[["projMatrix"]]

fgraphCPCA6 <- CreateDimReducObject(
  embeddings = fgComponents6,
  loadings = gCPCA[["projMatrix"]],
  stdev = numeric(),
  key = "gCPCA_",
  assay = "SCT"
)

rownames(fgraphCPCA6@feature.loadings) <- rownames(fdata6)
abnormalTissue6@reductions[["graphCPCA"]] <- fgraphCPCA6

abnormalTissue6 <- RunUMAP(abnormalTissue6, reduction = "graphCPCA", dims = 1:10, reduction.name = "gCPCA_UMAP")
abnormalTissue6 <- FindNeighbors(abnormalTissue6, reduction = "graphCPCA", dims = 1:6)
abnormalTissue6 <- FindClusters(abnormalTissue6, verbose = TRUE, resolution = 0.2)
abnormalTissue6@meta.data[["gcpca_clusters"]]  <- abnormalTissue6@meta.data[["SCT_snn_res.0.2"]]

# DimPlot(abnormalTissue6, reduction = "gCPCA_UMAP", label = FALSE, cols = c("#E1BD6D","#ED0000FF","#0B775E"))
SpatialDimPlot(abnormalTissue6, label = FALSE, label.size = 3, group.by = "gcpca_clusters", pt.size.factor = 2) +  
  scale_fill_manual(values=c("#E1BD6D","#ED0000FF","#0B775E"))


# set up query with the RCTD function SpatialRNA
STcounts6 <- abnormalTissue6@assays[["Spatial"]]@counts
coords6 <- GetTissueCoordinates(abnormalTissue6)
colnames(coords6) <- c("x", "y")
coords6[is.na(colnames(coords6))] <- NULL
query6 <- SpatialRNA(coords6, STcounts6, colSums(STcounts6))

RCTD6 <- create.RCTD(query6, reference, max_cores = 1)
RCTD6 <- run.RCTD(RCTD6, doublet_mode = "doublet")
abnormalTissue6 <- AddMetaData(abnormalTissue6, metadata = RCTD6@results$results_df)

SpatialDimPlot(abnormalTissue6, group.by = "first_type")
SpatialDimPlot(abnormalTissue6, group.by = "second_type") +  scale_fill_manual(values=col2)
                                    

# perform trajectory analysis (tissue wise)
sim6 <- SingleCellExperiment(assays = STcounts6)
reducedDims(sim6) <- SimpleList(DRM = fgComponents6)
colData(sim6)$clusterlabel <- factor(abnormalTissue6@meta.data[["gcpca_clusters"]])    
sim6  <-slingshot(sim6, clusterLabels = 'clusterlabel', reducedDim = 'DRM', start.clus="0") 
# in this data we set normal epithelial region as start cluster 

summary(sim6@colData@listData)
pseudotime_traj6 = sim6@colData@listData$slingPseudotime_1 # in this data only one trajectory was inferred
# gridnum = 10
# # color_in = c("#9C9EDE" ,"#5CB85C" ,"#E377C2", "#4DBBD5" ,"#FED439" ,"#FF9896", "#FFDC91","black")
# color_in = c("#E1BD6D","deepskyblue1","#7A3A9A","#ED0000FF","#0B775E","#ff00ff","black")
p_traj6 = plotTrajectory(pseudotime_traj6, coords6,abnormalTissue6@meta.data[["gcpca_clusters"]],gridnum,
                         color_in,pointsize=1.5 ,arrowlength=0.2,arrowsize=1,textsize=15 )
p_traj6$Arrowoverlay1
p_traj6$Pseudotime


# cell-ell communication
sub_abnormalTissue6 <- subset(abnormalTissue6, second_type != "NA")
data.input6 = GetAssayData(sub_abnormalTissue6, slot = "data", assay = "SCT") # normalized data matrix
# meta = data.frame(labels = as.numeric(Idents(abnormalTissue)), row.names = colnames(abnormalTissue))
meta6 = data.frame(labels = sub_abnormalTissue6@meta.data[["second_type"]], row.names = colnames(sub_abnormalTissue6))
meta6$labels <- as.factor(meta6$labels)
# meta$labels[is.na(meta$labels)] <- "Endothelial cells"
# meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))
# meta6$labels <- factor(meta6$labels, levels = c("Endothelial cells","Epithelial cells","Fibroblasts",
#                                                 "Macrophages","cDC","Prolif_Macro","B cells","T cells",
#                                                 "Proliferating T cells","pDC","Neutrophils","Plasma cells",
#                                                 "Monocytes","NK cells"))
meta6$labels <- factor(meta6$labels, levels = unique(meta6$labels))
scale.factors6 = jsonlite::fromJSON(txt = file.path("X:/maminu/Processed 129S4 Urethane model/26934/outs/spatial", 
                                                    'scalefactors_json.json'))
scale.factors6 = list(spot.diameter = 65, spot = scale.factors6$spot_diameter_fullres, # these two information are required
                      fiducial = scale.factors6$fiducial_diameter_fullres, hires = scale.factors6$tissue_hires_scalef, 
                      lowres = scale.factors6$tissue_lowres_scalef # these three information are not required
)

cellchat6 <- createCellChat(object = data.input6, meta = meta6, group.by = "labels",
                            datatype = "spatial", coordinates = GetTissueCoordinates(sub_abnormalTissue6), 
                            scale.factors = scale.factors6)

# CellChatDB <- CellChatDB.mouse
cellchat6@DB <- CellChatDB

cellchat6 <- subsetData(cellchat6) # This step is necessary even if using the whole database
future::plan("multisession", workers = 1)

cellchat6 <- identifyOverExpressedGenes(cellchat6)
cellchat6 <- identifyOverExpressedInteractions(cellchat6)

cellchat6 <- projectData(cellchat6, PPI.mouse)

cellchat6 <- computeCommunProb(cellchat6, type = "truncatedMean", trim = 0.1, 
                               distance.use = TRUE, interaction.length = 200, scale.distance = 1)

cellchat6 <- filterCommunication(cellchat6, min.cells = 0)
cellchat6 <- computeCommunProbPathway(cellchat6)
cellchat6 <- aggregateNet(cellchat6)

groupSize6 <- as.numeric(table(meta6$labels))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat6@net$count, vertex.weight = rowSums(cellchat6@net$count), weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions",
                 color.use = c("deepskyblue1","#ff00ff","plum1","#762A83","#DC0000FF","#4E84C4","tomato","lightblue1",
                                             "#C4961A","chartreuse3","#D16103"))
netVisual_circle(cellchat6@net$weight, vertex.weight = rowSums(cellchat6@net$weight), weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength",
                 color.use = c("deepskyblue1","#ff00ff","plum1","#762A83","#DC0000FF","#4E84C4","tomato","lightblue1",
                                             "#C4961A","chartreuse3","#D16103"))

mat6 <- cellchat6@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat6)) {
  mat66 <- matrix(0, nrow = nrow(mat6), ncol = ncol(mat6), dimnames = dimnames(mat6))
  mat66[i, ] <- mat6[i, ]
  netVisual_circle(mat66, vertex.weight = groupSize6, weight.scale = T, edge.weight.max = max(mat6), title.name = rownames(mat6)[i],
                   color.use = c("deepskyblue1","#ff00ff","plum1","#762A83","#DC0000FF","#4E84C4","tomato","lightblue1",
                                               "#C4961A","chartreuse3","#D16103"))
}

cellchat6@netP$pathways

pathways.show <- c("WNT") 
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat6, signaling = pathways.show, layout = "circle",
                    color.use = c("deepskyblue1","#ff00ff","plum1","#762A83","#DC0000FF","#4E84C4","tomato","lightblue1",
                                                "#C4961A","chartreuse3","#D16103"))

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat6, signaling = pathways.show, layout = "chord", vertex.label.cex = 0.8,
                    color.use = c("deepskyblue1","#ff00ff","plum1","#762A83","#DC0000FF","#4E84C4","tomato","lightblue1",
                                                "#C4961A","chartreuse3","#D16103"))

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat6, signaling = pathways.show, color.heatmap = "Reds",
                  color.use = c("deepskyblue1","#ff00ff","plum1","#762A83","#DC0000FF","#4E84C4","tomato","lightblue1",
                                              "#C4961A","chartreuse3","#D16103"))

# Spatial plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat6, signaling = pathways.show, layout = "spatial", edge.width.max = 2, 
                    vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5,
                    color.use = c("deepskyblue1","#ff00ff","plum1","#762A83","#DC0000FF","#4E84C4","tomato","lightblue1",
                                                "#C4961A","chartreuse3","#D16103"))

# Compute the network centrality scores
cellchat6 <- netAnalysis_computeCentrality(cellchat6, slot.name = "netP") 
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat6, signaling = pathways.show, width = 8, height = 2.5, font.size = 10,
                                  color.use = c("deepskyblue1","#ff00ff","plum1","#762A83","#DC0000FF","#4E84C4","tomato","lightblue1",
                                                              "#C4961A","chartreuse3","#D16103"))




# load seventh foreground data 
abnormalTissue7 <- Load10X_Spatial(
  data.dir = "X:/maminu/Processed 129S4 Urethane model/53434/outs", 
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", # specify name of the initial assay
  filter.matrix = TRUE, 
  to.upper = FALSE
)

abnormalTissue7 <- SCTransform(abnormalTissue7, assay = "Spatial", verbose = FALSE)
SpatialFeaturePlot(abnormalTissue7, features = "nCount_Spatial") + theme(legend.position = "right")

fdata7 <- abnormalTissue7@assays[["SCT"]]@scale.data
# flocation2 <- abnormalTissue2@images[["slice1"]]@coordinates[,c(2,3)] 
flocation7 <- GetTissueCoordinates(abnormalTissue7)

fgComponents7 <- t(fdata7) %*% gCPCA[["projMatrix"]]

fgraphCPCA7 <- CreateDimReducObject(
  embeddings = fgComponents7,
  loadings = gCPCA[["projMatrix"]],
  stdev = numeric(),
  key = "gCPCA_",
  assay = "SCT"
)

rownames(fgraphCPCA7@feature.loadings) <- rownames(fdata7)
abnormalTissue7@reductions[["graphCPCA"]] <- fgraphCPCA7

abnormalTissue7 <- RunUMAP(abnormalTissue7, reduction = "graphCPCA", dims = 1:50, reduction.name = "gCPCA_UMAP",
                           n.components = 30)
abnormalTissue7 <- FindNeighbors(abnormalTissue7, reduction = "graphCPCA", dims = 1:50)
abnormalTissue7 <- FindClusters(abnormalTissue7, verbose = TRUE, resolution = 0.5)
abnormalTissue7@meta.data[["gcpca_clusters"]]  <- abnormalTissue7@meta.data[["SCT_snn_res.0.5"]]

# DimPlot(abnormalTissue7, reduction = "gCPCA_UMAP", label = FALSE, 
#         cols = c("#E1BD6D","deepskyblue1","#7A3A9A","#ED0000FF","#0B775E","#ff00ff","#7B556C","#44B05B"))
SpatialDimPlot(abnormalTissue7, label = FALSE, label.size = 3, group.by = "gcpca_clusters", pt.size.factor = 2) +  
  scale_fill_manual(values=c("#9A8822","#E1BD6D","#ED0000FF","#0B775E","#4E2A1E","#44B05B","#ff00ff","#7B556C","deepskyblue1"),
                    labels = c("Lymphatic domain","Normal lung","Bronchus/Alveoli","Adenoma","Adenocarcinoma","Unknown domain"))
  


# set up query with the RCTD function SpatialRNA
STcounts7 <- abnormalTissue7@assays[["Spatial"]]@counts
coords7 <- GetTissueCoordinates(abnormalTissue7)
colnames(coords7) <- c("x", "y")
coords7[is.na(colnames(coords7))] <- NULL
query7 <- SpatialRNA(coords7, STcounts7, colSums(STcounts7))

RCTD7 <- create.RCTD(query7, reference, max_cores = 1)
RCTD7 <- run.RCTD(RCTD7, doublet_mode = "doublet")
abnormalTissue7 <- AddMetaData(abnormalTissue7, metadata = RCTD7@results$results_df)

SpatialDimPlot(abnormalTissue7, group.by = "first_type")
SpatialDimPlot(abnormalTissue7, group.by = "second_type") +  scale_fill_manual(values=col2)
                                    
# perform trajectory analysis (tissue wise)
sim7 <- SingleCellExperiment(assays = STcounts7)
reducedDims(sim7) <- SimpleList(DRM = fgComponents7)
colData(sim7)$clusterlabel <- factor(abnormalTissue7@meta.data[["gcpca_clusters"]])    
sim7  <-slingshot(sim7, clusterLabels = 'clusterlabel', reducedDim = 'DRM', start.clus="0") 
# in this data we set normal epithelial region as start cluster 

summary(sim7@colData@listData)
pseudotime_traj7 = sim7@colData@listData$slingPseudotime_1 # in this data only one trajectory was inferred
# gridnum = 10
# # color_in = c("#9C9EDE" ,"#5CB85C" ,"#E377C2", "#4DBBD5" ,"#FED439" ,"#FF9896", "#FFDC91","black")
# color_in = c("#E1BD6D","deepskyblue1","#7A3A9A","#ED0000FF","#0B775E","#ff00ff","black")
p_traj7 = plotTrajectory(pseudotime_traj7, coords7,abnormalTissue7@meta.data[["gcpca_clusters"]],gridnum,
                         color_in,pointsize=1.5 ,arrowlength=0.2,arrowsize=1,textsize=15 )
p_traj7$Arrowoverlay1
p_traj7$Pseudotime


# cell-ell communication
sub_abnormalTissue7 <- subset(abnormalTissue7, second_type != "NA")
data.input7 = GetAssayData(sub_abnormalTissue7, slot = "data", assay = "SCT") # normalized data matrix
# meta = data.frame(labels = as.numeric(Idents(abnormalTissue)), row.names = colnames(abnormalTissue))
meta7 = data.frame(labels = sub_abnormalTissue7@meta.data[["second_type"]], row.names = colnames(sub_abnormalTissue7))
meta7$labels <- as.factor(meta7$labels)
# meta$labels[is.na(meta$labels)] <- "Endothelial cells"
# meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))
# meta7$labels <- factor(meta7$labels, levels = c("Endothelial cells","Epithelial cells","Fibroblasts",
#                                                 "Macrophages","cDC","Prolif_Macro","B cells","T cells",
#                                                 "Proliferating T cells","pDC","Neutrophils","Plasma cells",
#                                                 "Monocytes","NK cells"))
meta7$labels <- factor(meta7$labels, levels = unique(meta7$labels))
scale.factors7 = jsonlite::fromJSON(txt = file.path("X:/maminu/Processed 129S4 Urethane model/53434/outs/spatial", 
                                                    'scalefactors_json.json'))
scale.factors7 = list(spot.diameter = 65, spot = scale.factors7$spot_diameter_fullres, # these two information are required
                      fiducial = scale.factors7$fiducial_diameter_fullres, hires = scale.factors7$tissue_hires_scalef, 
                      lowres = scale.factors7$tissue_lowres_scalef # these three information are not required
)

cellchat7 <- createCellChat(object = data.input7, meta = meta7, group.by = "labels",
                            datatype = "spatial", coordinates = GetTissueCoordinates(sub_abnormalTissue7), 
                            scale.factors = scale.factors7)

# CellChatDB <- CellChatDB.mouse
cellchat7@DB <- CellChatDB

cellchat7 <- subsetData(cellchat7) # This step is necessary even if using the whole database
future::plan("multisession", workers = 1)

cellchat7 <- identifyOverExpressedGenes(cellchat7)
cellchat7 <- identifyOverExpressedInteractions(cellchat7)

cellchat7 <- projectData(cellchat7, PPI.mouse)

cellchat7 <- computeCommunProb(cellchat7, type = "truncatedMean", trim = 0.1, 
                               distance.use = TRUE, interaction.length = 200, scale.distance = 1)

cellchat7 <- filterCommunication(cellchat7, min.cells = 0)
cellchat7 <- computeCommunProbPathway(cellchat7)
cellchat7 <- aggregateNet(cellchat7)

groupSize7 <- as.numeric(table(meta7$labels))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat7@net$count, vertex.weight = rowSums(cellchat7@net$count), weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions",
                 color.use = c("plum1","#762A83","#DC0000FF","deepskyblue1","#4E84C4","tomato","#ff00ff","lightblue1",
                                      "#C4961A","chartreuse3","#D16103","#068105"))
netVisual_circle(cellchat7@net$weight, vertex.weight = rowSums(cellchat7@net$weight), weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength",
                 color.use = c("plum1","#762A83","#DC0000FF","deepskyblue1","#4E84C4","tomato","#ff00ff","lightblue1",
                                      "#C4961A","chartreuse3","#D16103","#068105"))

mat7 <- cellchat7@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat7)) {
  mat77 <- matrix(0, nrow = nrow(mat7), ncol = ncol(mat7), dimnames = dimnames(mat7))
  mat77[i, ] <- mat7[i, ]
  netVisual_circle(mat77, vertex.weight = groupSize7, weight.scale = T, edge.weight.max = max(mat7), title.name = rownames(mat7)[i],
                   color.use = c("plum1","#762A83","#DC0000FF","deepskyblue1","#4E84C4","tomato","#ff00ff","lightblue1",
                                        "#C4961A","chartreuse3","#D16103","#068105"))
}

cellchat7@netP$pathways

pathways.show <- c("WNT") 
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat7, signaling = pathways.show, layout = "circle",
                    color.use = c("plum1","#762A83","#DC0000FF","deepskyblue1","#4E84C4","tomato","#ff00ff","lightblue1",
                                         "#C4961A","chartreuse3","#D16103","#068105"))

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat7, signaling = pathways.show, layout = "chord", vertex.label.cex = 0.8,
                    color.use = c("plum1","#762A83","#DC0000FF","deepskyblue1","#4E84C4","tomato","#ff00ff","lightblue1",
                                         "#C4961A","chartreuse3","#D16103","#068105"))

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat7, signaling = pathways.show, color.heatmap = "Reds",
                  color.use = c("plum1","#762A83","#DC0000FF","deepskyblue1","#4E84C4","tomato","#ff00ff","lightblue1",
                                       "#C4961A","chartreuse3","#D16103","#068105"))

# Spatial plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat7, signaling = pathways.show, layout = "spatial", edge.width.max = 2, 
                    vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5,
                    color.use = c("plum1","#762A83","#DC0000FF","deepskyblue1","#4E84C4","tomato","#ff00ff","lightblue1",
                                         "#C4961A","chartreuse3","#D16103","#068105"))

# Compute the network centrality scores
cellchat7 <- netAnalysis_computeCentrality(cellchat7, slot.name = "netP") 
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat7, signaling = pathways.show, width = 8, height = 2.5, font.size = 10,
                                  color.use = c("plum1","#762A83","#DC0000FF","deepskyblue1","#4E84C4","tomato","#ff00ff","lightblue1",
                                                       "#C4961A","chartreuse3","#D16103","#068105"))



# load eight foreground data 
abnormalTissue8 <- Load10X_Spatial(
  data.dir = "X:/maminu/Processed 129S4 Urethane model/26935/outs", 
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", # specify name of the initial assay
  filter.matrix = TRUE, 
  to.upper = FALSE
)

abnormalTissue8 <- SCTransform(abnormalTissue8, assay = "Spatial", verbose = FALSE)
SpatialFeaturePlot(abnormalTissue8, features = "nCount_Spatial") + theme(legend.position = "right")

fdata8 <- abnormalTissue8@assays[["SCT"]]@scale.data
# flocation2 <- abnormalTissue2@images[["slice1"]]@coordinates[,c(2,3)] 
flocation8 <- GetTissueCoordinates(abnormalTissue8)

fgComponents8 <- t(fdata8) %*% gCPCA[["projMatrix"]]

fgraphCPCA8 <- CreateDimReducObject(
  embeddings = fgComponents8,
  loadings = gCPCA[["projMatrix"]],
  stdev = numeric(),
  key = "gCPCA_",
  assay = "SCT"
)

rownames(fgraphCPCA8@feature.loadings) <- rownames(fdata8)
abnormalTissue8@reductions[["graphCPCA"]] <- fgraphCPCA8

abnormalTissue8 <- RunUMAP(abnormalTissue8, reduction = "graphCPCA", dims = 1:50, reduction.name = "gCPCA_UMAP")
abnormalTissue8 <- FindNeighbors(abnormalTissue8, reduction = "graphCPCA", dims = 1:20)
abnormalTissue8 <- FindClusters(abnormalTissue8, verbose = TRUE, resolution = 0.14)
abnormalTissue8@meta.data[["gcpca_clusters"]]  <- abnormalTissue8@meta.data[["SCT_snn_res.0.14"]]

# DimPlot(abnormalTissue8, reduction = "gCPCA_UMAP", label = FALSE, 
#         cols = c("#E1BD6D","deepskyblue1","#7A3A9A","#ED0000FF","#0B775E","#ff00ff","#7B556C","#44B05B"))
SpatialDimPlot(abnormalTissue8, label = FALSE, label.size = 3, group.by = "gcpca_clusters", pt.size.factor = 2.2) +  
  scale_fill_manual(values=c("#E1BD6D","#0B775E","#ED0000FF","deepskyblue1","#7A3A9A","#ff00ff","#7B556C","#44B05B"))


# set up query with the RCTD function SpatialRNA
STcounts8 <- abnormalTissue8@assays[["Spatial"]]@counts
coords8 <- GetTissueCoordinates(abnormalTissue8)
colnames(coords8) <- c("x", "y")
coords8[is.na(colnames(coords8))] <- NULL
query8 <- SpatialRNA(coords8, STcounts8, colSums(STcounts8))

RCTD8 <- create.RCTD(query8, reference, max_cores = 1)
RCTD8 <- run.RCTD(RCTD8, doublet_mode = "doublet")
abnormalTissue8 <- AddMetaData(abnormalTissue8, metadata = RCTD8@results$results_df)

SpatialDimPlot(abnormalTissue8, group.by = "first_type")
SpatialDimPlot(abnormalTissue8, group.by = "second_type") +  scale_fill_manual(values=col2)
                                    

# perform trajectory analysis (tissue wise)
sim8 <- SingleCellExperiment(assays = STcounts8)
reducedDims(sim8) <- SimpleList(DRM = fgComponents8)
colData(sim8)$clusterlabel <- factor(abnormalTissue8@meta.data[["gcpca_clusters"]])    
sim8  <-slingshot(sim8, clusterLabels = 'clusterlabel', reducedDim = 'DRM', start.clus="0") 
# in this data we set normal epithelial region as start cluster 

summary(sim8@colData@listData)
pseudotime_traj8 = sim8@colData@listData$slingPseudotime_1 # in this data only one trajectory was inferred
# gridnum = 10
# # color_in = c("#9C9EDE" ,"#5CB85C" ,"#E377C2", "#4DBBD5" ,"#FED439" ,"#FF9896", "#FFDC91","black")
# color_in = c("#E1BD6D","deepskyblue1","#7A3A9A","#ED0000FF","#0B775E","#ff00ff","black")
p_traj8 = plotTrajectory(pseudotime_traj8, coords8,abnormalTissue8@meta.data[["gcpca_clusters"]],gridnum,
                         color_in,pointsize=1.5 ,arrowlength=0.2,arrowsize=1,textsize=15 )
p_traj8$Arrowoverlay1
p_traj8$Pseudotime


# cell-ell communication
# sub_abnormalTissue8 <- subset(abnormalTissue8, second_type != "NA")
data.input8 = GetAssayData(abnormalTissue8, slot = "data", assay = "SCT") # normalized data matrix
# meta = data.frame(labels = as.numeric(Idents(abnormalTissue)), row.names = colnames(abnormalTissue))
meta8 = data.frame(labels = abnormalTissue8@meta.data[["second_type"]], row.names = colnames(abnormalTissue8))
meta8$labels <- as.factor(meta8$labels)
# meta$labels[is.na(meta$labels)] <- "Endothelial cells"
# meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))
# meta8$labels <- factor(meta8$labels, levels = c("Endothelial cells","Epithelial cells","Fibroblasts",
#                                                 "Macrophages","cDC","Prolif_Macro","B cells","T cells",
#                                                 "Proliferating T cells","pDC","Neutrophils","Plasma cells",
#                                                 "Monocytes","NK cells"))
meta8$labels <- factor(meta8$labels, levels = unique(meta8$labels))
scale.factors8 = jsonlite::fromJSON(txt = file.path("X:/maminu/Processed 129S4 Urethane model/26935/outs/spatial", 
                                                    'scalefactors_json.json'))
scale.factors8 = list(spot.diameter = 65, spot = scale.factors8$spot_diameter_fullres, # these two information are required
                      fiducial = scale.factors8$fiducial_diameter_fullres, hires = scale.factors8$tissue_hires_scalef, 
                      lowres = scale.factors8$tissue_lowres_scalef # these three information are not required
)

cellchat8 <- createCellChat(object = data.input8, meta = meta8, group.by = "labels",
                            datatype = "spatial", coordinates = GetTissueCoordinates(abnormalTissue8), 
                            scale.factors = scale.factors8)

# CellChatDB <- CellChatDB.mouse
cellchat8@DB <- CellChatDB

cellchat8 <- subsetData(cellchat8) # This step is necessary even if using the whole database
future::plan("multisession", workers = 1)

cellchat8 <- identifyOverExpressedGenes(cellchat8)
cellchat8 <- identifyOverExpressedInteractions(cellchat8)

cellchat8 <- projectData(cellchat8, PPI.mouse)

cellchat8 <- computeCommunProb(cellchat8, type = "truncatedMean", trim = 0.1, 
                               distance.use = TRUE, interaction.length = 200, scale.distance = 1)

cellchat8 <- filterCommunication(cellchat8, min.cells = 0)
cellchat8 <- computeCommunProbPathway(cellchat8)
cellchat8 <- aggregateNet(cellchat8)

groupSize8 <- as.numeric(table(meta8$labels))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat8@net$count, vertex.weight = rowSums(cellchat8@net$count), weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions",
                 color.use = c("#762A83","plum1","deepskyblue1","#DC0000FF","lightblue1","tomato","#ff00ff","chartreuse3",
                                        "#4E84C4","#D16103","#C4961A"))
netVisual_circle(cellchat8@net$weight, vertex.weight = rowSums(cellchat8@net$weight), weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength",
                 color.use = c("#762A83","plum1","deepskyblue1","#DC0000FF","lightblue1","tomato","#ff00ff","chartreuse3",
                                        "#4E84C4","#D16103","#C4961A"))

mat8 <- cellchat8@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat8)) {
  mat88 <- matrix(0, nrow = nrow(mat8), ncol = ncol(mat8), dimnames = dimnames(mat8))
  mat88[i, ] <- mat8[i, ]
  netVisual_circle(mat88, vertex.weight = groupSize8, weight.scale = T, edge.weight.max = max(mat8), title.name = rownames(mat8)[i],
                   color.use = c("#762A83","plum1","deepskyblue1","#DC0000FF","lightblue1","tomato","#ff00ff","chartreuse3",
                                          "#4E84C4","#D16103","#C4961A"))
}

cellchat8@netP$pathways

pathways.show <- c("WNT") 
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat8, signaling = pathways.show, layout = "circle",
                    color.use = c("#762A83","plum1","deepskyblue1","#DC0000FF","lightblue1","tomato","#ff00ff","chartreuse3",
                                           "#4E84C4","#D16103","#C4961A"))

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat8, signaling = pathways.show, layout = "chord", vertex.label.cex = 0.8,
                    color.use = c("#762A83","plum1","deepskyblue1","#DC0000FF","lightblue1","tomato","#ff00ff","chartreuse3",
                                           "#4E84C4","#D16103","#C4961A"))

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat8, signaling = pathways.show, color.heatmap = "Reds",
                  color.use = c("#762A83","plum1","deepskyblue1","#DC0000FF","lightblue1","tomato","#ff00ff","chartreuse3",
                                         "#4E84C4","#D16103","#C4961A"))

# Spatial plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat8, signaling = pathways.show, layout = "spatial", edge.width.max = 2, 
                    vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5,
                    color.use = c("#762A83","plum1","deepskyblue1","#DC0000FF","lightblue1","tomato","#ff00ff","chartreuse3",
                                           "#4E84C4","#D16103","#C4961A"))

# Compute the network centrality scores
cellchat8 <- netAnalysis_computeCentrality(cellchat8, slot.name = "netP") 
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat8, signaling = pathways.show, width = 8, height = 2.5, font.size = 10,
                                  color.use = c("#762A83","plum1","deepskyblue1","#DC0000FF","lightblue1","tomato","#ff00ff","chartreuse3",
                                                         "#4E84C4","#D16103","#C4961A"))



# load ninth foreground data 
abnormalTissue9 <- Load10X_Spatial(
  data.dir = "X:/maminu/Processed 129S4 Urethane model/53435/outs", 
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial", # specify name of the initial assay
  filter.matrix = TRUE, 
  to.upper = FALSE
)

abnormalTissue9 <- SCTransform(abnormalTissue9, assay = "Spatial", verbose = FALSE)
SpatialFeaturePlot(abnormalTissue9, features = "nCount_Spatial") + theme(legend.position = "right")

fdata9 <- abnormalTissue9@assays[["SCT"]]@scale.data
# flocation2 <- abnormalTissue2@images[["slice1"]]@coordinates[,c(2,3)] 
flocation9 <- GetTissueCoordinates(abnormalTissue9)

fgComponents9 <- t(fdata9) %*% gCPCA[["projMatrix"]]

fgraphCPCA9 <- CreateDimReducObject(
  embeddings = fgComponents9,
  loadings = gCPCA[["projMatrix"]],
  stdev = numeric(),
  key = "gCPCA_",
  assay = "SCT"
)

rownames(fgraphCPCA9@feature.loadings) <- rownames(fdata9)
abnormalTissue9@reductions[["graphCPCA"]] <- fgraphCPCA9

abnormalTissue9 <- RunUMAP(abnormalTissue9, reduction = "graphCPCA", dims = 1:10, reduction.name = "gCPCA_UMAP")
abnormalTissue9 <- FindNeighbors(abnormalTissue9, reduction = "graphCPCA", dims = 1:10)
abnormalTissue9 <- FindClusters(abnormalTissue9, verbose = TRUE, resolution = 0.3)
abnormalTissue9@meta.data[["gcpca_clusters"]]  <- abnormalTissue9@meta.data[["SCT_snn_res.0.3"]]

# DimPlot(abnormalTissue9, reduction = "gCPCA_UMAP", label = FALSE, 
#         cols = c("#E1BD6D","deepskyblue1","#7A3A9A","#ED0000FF","#0B775E","#ff00ff","#7B556C","#44B05B"))
SpatialDimPlot(abnormalTissue9, label = FALSE, label.size = 3, group.by = "gcpca_clusters", pt.size.factor = 2) +  
  scale_fill_manual(values=c("#E1BD6D","#4E2A1E","#0B775E","deepskyblue1","#7A3A9A","#ED0000FF","#ff00ff","#7B556C","#44B05B"))


# set up query with the RCTD function SpatialRNA
STcounts9 <- abnormalTissue9@assays[["Spatial"]]@counts
coords9 <- GetTissueCoordinates(abnormalTissue9)
colnames(coords9) <- c("x", "y")
coords9[is.na(colnames(coords9))] <- NULL
query9 <- SpatialRNA(coords9, STcounts9, colSums(STcounts9))

RCTD9 <- create.RCTD(query9, reference, max_cores = 1)
RCTD9 <- run.RCTD(RCTD9, doublet_mode = "doublet")
abnormalTissue9 <- AddMetaData(abnormalTissue9, metadata = RCTD9@results$results_df)

SpatialDimPlot(abnormalTissue9, group.by = "first_type")
SpatialDimPlot(abnormalTissue9, group.by = "second_type") +  scale_fill_manual(values=col2)
                                    

# perform trajectory analysis (tissue wise)
sim9 <- SingleCellExperiment(assays = STcounts9)
reducedDims(sim9) <- SimpleList(DRM = fgComponents9)
colData(sim9)$clusterlabel <- factor(abnormalTissue9@meta.data[["gcpca_clusters"]])    
sim9  <-slingshot(sim9, clusterLabels = 'clusterlabel', reducedDim = 'DRM', start.clus="0") 
# in this data we set normal epithelial region as start cluster 

summary(sim9@colData@listData)
pseudotime_traj9 = sim9@colData@listData$slingPseudotime_1 # in this data only one trajectory was inferred
# gridnum = 10
# # color_in = c("#9C9EDE" ,"#5CB85C" ,"#E377C2", "#4DBBD5" ,"#FED439" ,"#FF9896", "#FFDC91","black")
# color_in = c("#E1BD6D","deepskyblue1","#7A3A9A","#ED0000FF","#0B775E","#ff00ff","black")
p_traj9 = plotTrajectory(pseudotime_traj9, coords9,abnormalTissue9@meta.data[["gcpca_clusters"]],gridnum,
                         color_in,pointsize=1.5 ,arrowlength=0.2,arrowsize=1,textsize=15 )
p_traj9$Arrowoverlay1
p_traj9$Pseudotime


# cell-ell communication
# sub_abnormalTissue9 <- subset(abnormalTissue9, second_type != "NA")
data.input9 = GetAssayData(abnormalTissue9, slot = "data", assay = "SCT") # normalized data matrix
# meta = data.frame(labels = as.numeric(Idents(abnormalTissue)), row.names = colnames(abnormalTissue))
meta9 = data.frame(labels = abnormalTissue9@meta.data[["second_type"]], row.names = colnames(abnormalTissue9))
meta9$labels <- as.factor(meta9$labels)
# meta$labels[is.na(meta$labels)] <- "Endothelial cells"
# meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))
# meta9$labels <- factor(meta9$labels, levels = c("Endothelial cells","Epithelial cells","Fibroblasts",
#                                                 "Macrophages","cDC","Prolif_Macro","B cells","T cells",
#                                                 "Proliferating T cells","pDC","Neutrophils","Plasma cells",
#                                                 "Monocytes","NK cells"))
meta9$labels <- factor(meta9$labels, levels = unique(meta9$labels))
scale.factors9 = jsonlite::fromJSON(txt = file.path("X:/maminu/Processed 129S4 Urethane model/53435/outs/spatial", 
                                                    'scalefactors_json.json'))
scale.factors9 = list(spot.diameter = 65, spot = scale.factors9$spot_diameter_fullres, # these two information are required
                      fiducial = scale.factors9$fiducial_diameter_fullres, hires = scale.factors9$tissue_hires_scalef, 
                      lowres = scale.factors9$tissue_lowres_scalef # these three information are not required
)

cellchat9 <- createCellChat(object = data.input9, meta = meta9, group.by = "labels",
                            datatype = "spatial", coordinates = GetTissueCoordinates(abnormalTissue9), 
                            scale.factors = scale.factors9)

# CellChatDB <- CellChatDB.mouse
cellchat9@DB <- CellChatDB

cellchat9 <- subsetData(cellchat9) # This step is necessary even if using the whole database
future::plan("multisession", workers = 1)

cellchat9 <- identifyOverExpressedGenes(cellchat9)
cellchat9 <- identifyOverExpressedInteractions(cellchat9)

cellchat9 <- projectData(cellchat9, PPI.mouse)

cellchat9 <- computeCommunProb(cellchat9, type = "truncatedMean", trim = 0.1, 
                               distance.use = TRUE, interaction.length = 200, scale.distance = 1)

cellchat9 <- filterCommunication(cellchat9, min.cells = 0)
cellchat9 <- computeCommunProbPathway(cellchat9)
cellchat9 <- aggregateNet(cellchat9)

groupSize9 <- as.numeric(table(meta9$labels))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat9@net$count, vertex.weight = rowSums(cellchat9@net$count), weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions",
                 color.use = c("plum1","deepskyblue1","#762A83","lightblue1","#ff00ff","#C4961A","tomato","chartreuse3",
                                      "#4E84C4","#DC0000FF","#D16103"))
netVisual_circle(cellchat9@net$weight, vertex.weight = rowSums(cellchat9@net$weight), weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength",
                 color.use = c("plum1","deepskyblue1","#762A83","lightblue1","#ff00ff","#C4961A","tomato","chartreuse3",
                                      "#4E84C4","#DC0000FF","#D16103"))

mat9 <- cellchat9@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat9)) {
  mat99 <- matrix(0, nrow = nrow(mat9), ncol = ncol(mat9), dimnames = dimnames(mat9))
  mat99[i, ] <- mat9[i, ]
  netVisual_circle(mat99, vertex.weight = groupSize9, weight.scale = T, edge.weight.max = max(mat9), title.name = rownames(mat9)[i],
                   color.use = c("plum1","deepskyblue1","#762A83","lightblue1","#ff00ff","#C4961A","tomato","chartreuse3",
                                        "#4E84C4","#DC0000FF","#D16103"))
}

cellchat9@netP$pathways

pathways.show <- c("WNT") 
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat9, signaling = pathways.show, layout = "circle",
                    color.use = c("plum1","deepskyblue1","#762A83","lightblue1","#ff00ff","#C4961A","tomato","chartreuse3",
                                         "#4E84C4","#DC0000FF","#D16103"))

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat9, signaling = pathways.show, layout = "chord", vertex.label.cex = 0.8,
                    color.use = c("plum1","deepskyblue1","#762A83","lightblue1","#ff00ff","#C4961A","tomato","chartreuse3",
                                         "#4E84C4","#DC0000FF","#D16103"))

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat9, signaling = pathways.show, color.heatmap = "Reds",
                  color.use = c("plum1","deepskyblue1","#762A83","lightblue1","#ff00ff","#C4961A","tomato","chartreuse3",
                                       "#4E84C4","#DC0000FF","#D16103"))

# Spatial plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat9, signaling = pathways.show, layout = "spatial", edge.width.max = 2, 
                    vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5,
                    color.use = c("plum1","deepskyblue1","#762A83","lightblue1","#ff00ff","#C4961A","tomato","chartreuse3",
                                         "#4E84C4","#DC0000FF","#D16103"))

# Compute the network centrality scores
cellchat9 <- netAnalysis_computeCentrality(cellchat9, slot.name = "netP") 
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat9, signaling = pathways.show, width = 8, height = 2.5, font.size = 10,
                                  color.use = c("plum1","deepskyblue1","#762A83","lightblue1","#ff00ff","#C4961A","tomato","chartreuse3",
                                                       "#4E84C4","#DC0000FF","#D16103"))

# Aggregate interaction weights form all samples

library(reshape2)
nmat1 <- as.data.frame(mat1)
nmat1$celltypes <- rownames(mat1)
nmat2 <- as.data.frame(mat2)
nmat2$celltypes <- rownames(mat2)
nmat3 <- as.data.frame(mat3)
nmat3$celltypes <- rownames(mat3)
nmat4 <- as.data.frame(mat4)
nmat4$celltypes <- rownames(mat4)
nmat5 <- as.data.frame(mat5)
nmat5$celltypes <- rownames(mat5)
nmat6 <- as.data.frame(mat6)
nmat6$celltypes <- rownames(mat6)
nmat7 <- as.data.frame(mat7)
nmat7$celltypes <- rownames(mat7)
nmat8 <- as.data.frame(mat8)
nmat8$celltypes <- rownames(mat8)
nmat9 <- as.data.frame(mat9)
nmat9$celltypes <- rownames(mat9)

finalMat <- rbind(melt(nmat1,id="celltypes"),melt(nmat2,id="celltypes"),melt(nmat3,id="celltypes"),melt(nmat4,id="celltypes"),
                  melt(nmat5,id="celltypes"),melt(nmat6,id="celltypes"),melt(nmat7,id="celltypes"),melt(nmat8,id="celltypes"),
                  melt(nmat9,id="celltypes"))
finalMat <- dcast(data=finalMat,formula=celltypes ~ variable,fun.aggregate=sum)
rownames(finalMat) <- finalMat$celltypes
finalMat <- finalMat[,-c(1)]
finalMat <- as.matrix(finalMat[colnames(finalMat),])

grpSize <- as.numeric(table(rbind(meta,meta2,meta3,meta4,meta5,meta6,meta7,meta8,meta9)))

par(mfrow=c(1,1))
netVisual_circle(finalMat, vertex.weight = rowSums(finalMat), weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength",
                 color.use = c("plum1","tomato","#762A83","deepskyblue1","#C4961A","#ff00ff","#DC0000FF","#4E84C4","chartreuse3",
                                      "#D16103","#58593FFF","lightblue1","#068105"))
                                         
par(mfrow = c(4,4), xpd=TRUE)
for (i in 1:nrow(finalMat)) {
  finalMatt <- matrix(0, nrow = nrow(finalMat), ncol = ncol(finalMat), dimnames = dimnames(finalMat))
  finalMatt[i, ] <- finalMat[i, ]
  netVisual_circle(finalMatt, vertex.weight = grpSize, weight.scale = T, edge.weight.max = max(finalMat), title.name = rownames(finalMat)[i],
                   color.use = c("plum1","tomato","#762A83","deepskyblue1","#C4961A","#ff00ff","#DC0000FF","#4E84C4","chartreuse3",
                                        "#D16103","#58593FFF","lightblue1","#068105"))
}

# Adenoma tissues interaction weights
adenomaMat <- rbind(melt(nmat1,id="celltypes"),melt(nmat3,id="celltypes"),melt(nmat4,id="celltypes"),
                  melt(nmat5,id="celltypes"),melt(nmat6,id="celltypes"),melt(nmat8,id="celltypes"))
adenomaMat <- dcast(data=adenomaMat,formula=celltypes ~ variable,fun.aggregate=sum)
rownames(adenomaMat) <- adenomaMat$celltypes
adenomaMat <- adenomaMat[,-c(1)]
adenomaMat <- as.matrix(adenomaMat[colnames(adenomaMat),])

adenomaSize <- as.numeric(table(rbind(meta,meta3,meta4,meta5,meta6,meta8)))

par(mfrow=c(1,1))
netVisual_circle(adenomaMat, vertex.weight = rowSums(adenomaMat), weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength",
                 color.use = c("plum1","tomato","#762A83","deepskyblue1","#C4961A","#ff00ff","#DC0000FF","#4E84C4","chartreuse3",
                                      "#D16103","#58593FFF","lightblue1","#068105"))
                                      
par(mfrow = c(4,4), xpd=TRUE)
for (i in 1:nrow(adenomaMat)) {
  adenomaMatt <- matrix(0, nrow = nrow(adenomaMat), ncol = ncol(adenomaMat), dimnames = dimnames(adenomaMat))
  adenomaMatt[i, ] <- adenomaMat[i, ]
  netVisual_circle(adenomaMatt, vertex.weight = adenomaSize, weight.scale = T, edge.weight.max = max(adenomaMat), title.name = rownames(adenomaMat)[i],
                   color.use = c("plum1","tomato","#762A83","deepskyblue1","#C4961A","#ff00ff","#DC0000FF","#4E84C4","chartreuse3",
                                        "#D16103","#58593FFF","lightblue1","#068105"))
}

# Adenocarcinoma tissues interaction weights
adcMat <- rbind(melt(nmat7,id="celltypes"),melt(nmat9,id="celltypes"))
adcMat <- dcast(data=adcMat,formula=celltypes ~ variable,fun.aggregate=sum)
rownames(adcMat) <- adcMat$celltypes
adcMat <- adcMat[,-c(1)]
adcMat <- as.matrix(adcMat[colnames(adcMat),])

adcSize <- as.numeric(table(rbind(meta7,meta9)))

par(mfrow=c(1,1))
netVisual_circle(adcMat, vertex.weight = rowSums(adcMat), weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength",
                 color.use = c("plum1","#762A83","#DC0000FF","deepskyblue1","#4E84C4","tomato","#ff00ff","lightblue1","#C4961A",
                                      "chartreuse3","#D16103","#068105","#58593FFF"))
                                      
par(mfrow = c(4,4), xpd=TRUE)
for (i in 1:nrow(adcMat)) {
  adcMatt <- matrix(0, nrow = nrow(adcMat), ncol = ncol(adcMat), dimnames = dimnames(adcMat))
  adcMatt[i, ] <- adcMat[i, ]
  netVisual_circle(adcMatt, vertex.weight = adcSize, weight.scale = T, edge.weight.max = max(adcMat), title.name = rownames(adcMat)[i],
                   color.use = c("plum1","#762A83","#DC0000FF","deepskyblue1","#4E84C4","tomato","#ff00ff","lightblue1","#C4961A",
                                        "chartreuse3","#D16103","#068105","#58593FFF"))
}


# subset tumor regions for trajectory inference
normal1 <- subset(abnormalTissue, gcpca_clusters == "0")
tumor1 <- subset(abnormalTissue, gcpca_clusters == "4")
tumor3 <- subset(abnormalTissue3, gcpca_clusters == "2")
tumor4 <- subset(abnormalTissue4, gcpca_clusters == "2")
tumor5 <- subset(abnormalTissue5, gcpca_clusters == "3")
tumor6 <- subset(abnormalTissue6, gcpca_clusters == "2")
tumor7 <- subset(abnormalTissue7, gcpca_clusters == "4")
tumor8 <- subset(abnormalTissue8, gcpca_clusters == "1")
tumor9 <- subset(abnormalTissue9, gcpca_clusters == "1")

ttrajData <- merge(x = normal1, y = list(tumor1, tumor3, tumor4, tumor5, tumor6, tumor7, tumor8, tumor9))
Adenoma <- merge(x = tumor1, y = list(tumor3, tumor4, tumor5, tumor6, tumor8))
Adenocarcinoma <- merge(x = tumor7, y = tumor9)


rdata <- ttrajData@assays[["SCT"]]@counts
# Now all set to run trajectory analysis
trajData <- CreateSeuratObject(counts = rdata,
                               project = "PreCancer",
                               assay = "PreCancer")

trajData@meta.data[["labels"]] <- c(rep("Normal",ncol(normal1)),rep("Adenoma",ncol(tumor1)),rep("Adenoma",ncol(tumor3)),
                                    rep("Adenoma",ncol(tumor4)),rep("Adenoma",ncol(tumor5)),rep("Adenoma",ncol(tumor6)),
                                    rep("Adenocarcinoma",ncol(tumor7)),rep("Adenoma",ncol(tumor8)),
                                    rep("Adenocarcinoma",ncol(tumor9)))
trajData <- NormalizeData(trajData)

trajData <- FindVariableFeatures(trajData, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(trajData), 10)
top10

plot1 <- VariableFeaturePlot(trajData)
LabelPoints(plot = plot1, points = top10, repel = T)#+ylim(c(0,600000))

all.genes <- rownames(trajData)
trajData <- ScaleData(trajData, features = all.genes)

trajData <- RunPCA(trajData, features = VariableFeatures(object = trajData))

DimHeatmap(trajData, dims = 1:6, cells = 500, balanced = T)

ElbowPlot(trajData)

# Idents(trajData) <- gnd2
# trajData@meta.data[["labels"]] <- gnd2
# fdata@meta.data[["Slabels"]] <- gnd3

trajData <- RunUMAP(trajData, dims = 1:10, n.components = 5)
# trajData <- RunUMAP(trajData, dims = 20, n.neighbors = 30L)
# fdata <- RunTSNE(fdata, dims = 1:10)

DimPlot(trajData, reduction = "umap", group.by = "labels", cols = c("#4E2A1E","#0B775E","#E1BD6D"), label = F)
# DimPlot(fdata, reduction = "tsne")

# trajectory analysis with monocle
cds <- as.cell_data_set(trajData)
head(colData(cds))

fData(cds)
rownames(fData(cds))[1:10]

fData(cds)$gene_short_name <- rownames(fData(cds))
head(fData(cds))
head(counts(cds))

cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds, preprocess_method = 'PCA')

recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

list.cluster <- as.factor(trajData@meta.data[["labels"]])
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- trajData@reductions$umap@cell.embeddings

cluster.before.traj <-plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,
                                 group_label_size = 5) + theme(legend.position = "right")
# cluster.before.traj

cds <- learn_graph(cds,learn_graph_control=list(ncenter=180))

plot_cells(cds, color_cells_by = "labels", label_groups_by_cluster = F,
           label_branch_points = F, label_roots = F, label_leaves = F,
           group_label_size = 5)

cds <- order_cells(cds, reduction_method = "UMAP", 
                   root_cells = colnames(cds[, cds@clusters@listData[["UMAP"]][["clusters"]] == "Normal"]))
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = F, label_roots = F, label_leaves = F) & scale_color_viridis_c()

head(pseudotime(cds), 10)

cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

ggplot(data.pseudo, aes(monocle3_pseudotime, labels, fill = labels)) + geom_boxplot()

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(labels, monocle3_pseudotime), fill = labels)) + xlab("Pseudotime") + 
  geom_boxplot()

deg <- graph_test(cds, neighbor_graph = "principal_graph")
deg %>% arrange(q_value) %>% dplyr::filter(status == "OK") %>% head()

FeaturePlot(trajData, features = rownames(deg)[c(1:6)]) & scale_color_viridis_c()

trajData$pseudotime <- pseudotime(cds)
FeaturePlot(trajData, features = "pseudotime", pt.size = 0.7) & scale_color_viridis_c()

ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))

plot_cells(cds, genes=pr_deg_ids[c(1:4)],
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)

gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=0.00001)

cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$labels)

agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat, 
                   scale="column", clustering_method="ward.D2")

efea <- pr_deg_ids[c(41:50)]
lineage_cds <- cds[rowData(cds)$gene_short_name %in% efea,]

plot_genes_in_pseudotime(lineage_cds,
                         color_cells_by="labels",
                         min_expr = 1e-5)


# # cell-ell communication
# library(CellChat)
# library(patchwork)
# options(stringsAsFactors = FALSE)
# 
# dataNormal = GetAssayData(normal1, slot = "data", assay = "SCT") # normalized data matrix
# # meta = data.frame(labels = as.numeric(Idents(abnormalTissue)), row.names = colnames(abnormalTissue))
# metaNormal = data.frame(labels = normal1@meta.data[["second_type"]], row.names = colnames(normal1))
# metaNormal$labels <- as.factor(metaNormal$labels)
# # meta$labels[is.na(meta$labels)] <- "Endothelial cells"
# # meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))
# meta$labels <- factor(meta$labels, levels = c("B cells","cDC","Endothelial cells","Epithelial cells","Fibroblasts","Macrophages",
#                                               "Prolif_Macro","Proliferating T cells","T cells"))
# scale.factors = jsonlite::fromJSON(txt = file.path("X:/maminu/Processed 129S4 Urethane model/26935/outs/spatial", 
#                                                    'scalefactors_json.json'))
# scale.factors = list(spot.diameter = 65, spot = scale.factors$spot_diameter_fullres, # these two information are required
#                      fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef, 
#                      lowres = scale.factors$tissue_lowres_scalef # these three information are not required
# )
# 
# cellchat <- createCellChat(object = dataNormal, meta = meta, group.by = "labels",
#                            datatype = "spatial", coordinates = flocation, scale.factors = scale.factors)
# 
# CellChatDB <- CellChatDB.mouse
# cellchat@DB <- CellChatDB
# 
# cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
# future::plan("multisession", workers = 4)
# 
# cellchat <- identifyOverExpressedGenes(cellchat)
# cellchat <- identifyOverExpressedInteractions(cellchat)
# 
# cellchat <- projectData(cellchat, PPI.mouse)
# 
# cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
#                               distance.use = TRUE, interaction.length = 200, scale.distance = 1)
# 
# cellchat <- filterCommunication(cellchat, min.cells = 0)
# cellchat <- computeCommunProbPathway(cellchat)
# cellchat <- aggregateNet(cellchat)
# 
# # groupSize <- as.numeric(table(cellchat@idents))
# groupSize <- as.numeric(table(meta$labels))
# par(mfrow = c(1,2), xpd=TRUE)
# netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, 
#                  label.edge= F, title.name = "Number of interactions")
# netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, 
#                  label.edge= F, title.name = "Interaction weights/strength")
# 
# cellchat@netP$pathways
# 
# pathways.show <- c("IL1") 
# # Circle plot
# par(mfrow=c(1,1))
# netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
# 
# # Spatial plot
# par(mfrow=c(1,1))
# netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", edge.width.max = 2, 
#                     vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5)
# 
# # Compute the network centrality scores
# cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
# par(mfrow=c(1,1))
# netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
# 
# 
# save.image("C:/Projects/Contrastive PCA/ST/Mouse Data (Bo)/Final Figures/Workspace_AllSamples_Final2.rds")
