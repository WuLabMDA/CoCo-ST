###########################
preparedata = function(percentage,celltypenames){
  rownames(percentage) = paste0("Cluster",1:nrow(percentage))
  celltype = celltypenames
  colnames(percentage) = celltype
  percentage_vec = c(percentage)
  cluster_vec = c(rep(c(paste0("Cluster",1:nrow(percentage))),length(celltype)))
  CellType = c(rep(celltype,each=nrow(percentage)))
  datt = data.frame(cluster_vec, percentage_vec,CellType)
  datt$cluster_vec = factor(cluster_vec, level=paste0("Cluster",1:nrow(percentage)))
  return(datt)
}

makefigure = function(datt){
  # cbp2 = c("#FFDB6D", "#C4961A", "#F4EDCA", "tomato","#C3D7A4",  "#4E84C4","#52854C",
  #          "#D16103", "deepskyblue1", "cadetblue3","lightblue1","plum1","chartreuse3")
  # cbp2 = c("#FB786E", "#ED8F00", "#CFA200", "#99A800","#56BB00",  "#00BC56","#00C094",
  #                   "#06A6FF", "#DF70F8", "#FB61D7","#FF69AD")
  cbp2 <- c("plum1","tomato","#762A83","deepskyblue1","#C4961A","#ff00ff","#DC0000FF","#4E84C4","chartreuse3",
                   "#D16103","#58593FFF","lightblue1","#068105","yellow","#4E2A1E",
                   "#C3D7A4","black")
  p=ggplot(datt, aes(y = percentage_vec,
                     x = factor(cluster_vec ), fill = CellType)) +        ## global aes
    scale_fill_manual(values=cbp2)+
    geom_bar(position="stack", stat="identity",width=0.7,color="black") +
    theme_bw()+xlab("")+ylab("")+
    theme(plot.title = element_text(size = 20),
          text = element_text(size = 20),
          #axis.title = element_text(face="bold"),
          axis.text.x=element_text(size = 12,angle = 60,hjust = 1),
          #axis.text.x=element_blank(),
          legend.position = "right")# +
  return(p)
}

plotCluster = function(location, clusterlabel, pointsize=3,text_size=15 ,title_in,color_in,legend="none"){
  cluster = clusterlabel
  loc_x=location[,1]
  loc_y=location[,2]
  datt = data.frame(cluster, loc_x, loc_y)
  p = ggplot(datt, aes(x = location[,1], y = location[,2], color = cluster)) +
    geom_point( alpha = 1,size=pointsize) +
    scale_color_manual(values = color_in)+
    ggtitle(paste0(title_in))+
    theme_void()+
    theme(plot.title = element_text(size = text_size,  face = "bold"),
          text = element_text(size = text_size),
          #axis.title = element_text(face="bold"),
          #axis.text.x=element_text(size = 15) ,
          legend.position =legend) 
  p
}

plotTrajectory = function(pseudotime, location,clusterlabels,gridnum,color_in,pointsize=5 ,arrowlength=0.2,arrowsize=1,textsize=22){
  
  pseudotime_use=pseudotime
  info = as.data.frame(location)
  colnames(info) = c("sdimx","sdimy")
  grids = gridnum
  
  min_x = min(info$sdimx)
  min_y = min(info$sdimy)
  max_x = max(info$sdimx)
  max_y = max(info$sdimy)
  
  x_anchor = c()
  for(x_i in 1:(grids+1)){
    space_x = (max_x - min_x)/grids
    x_anchor[x_i] = min_x+(x_i-1)*space_x
  }
  y_anchor = c()
  for(y_i in 1:(grids+1)){
    space_y = (max_y - min_y)/grids
    y_anchor[y_i] = min_y+(y_i-1)*space_y
  }
  
  # label square by num_x, num_y
  count = 0
  squares = list()
  direction_pseudotime_point = list()
  start_x_dat = c()
  start_y_dat = c()
  end_x_dat = c()
  end_y_dat = c()
  for(num_x in 1:grids){
    for(num_y in 1:grids){
      
      filter_x = which(info$sdimx >= x_anchor[num_x] & info$sdimx <= x_anchor[num_x+1])
      filter_y = which(info$sdimy >= y_anchor[num_y] & info$sdimy <= y_anchor[num_y+1])
      # find points in each grid
      points_in_grid = intersect(filter_x, filter_y)
      
      
      # find min pseudotime and max pseudotime in each grid
      if(length(points_in_grid)>1 & sum(which.min(pseudotime_use[points_in_grid]))>0){
        count = count + 1
        squares[[count]]= intersect(filter_x, filter_y)
        direction_pseudotime_point[[count]] = list()
        direction_pseudotime_point[[count]]$min_point = info[squares[[count]][which.min(pseudotime_use[squares[[count]]])],]
        direction_pseudotime_point[[count]]$max_point = info[squares[[count]][which.max(pseudotime_use[squares[[count]]])],]
        start_x_dat[count] = unlist(direction_pseudotime_point[[count]]$min_point$sdimx)
        start_y_dat[count] = unlist(direction_pseudotime_point[[count]]$min_point$sdimy)
        end_x_dat[count] = unlist(direction_pseudotime_point[[count]]$max_point$sdimx)
        end_y_dat[count] = unlist(direction_pseudotime_point[[count]]$max_point$sdimy)
      }
    }
  }
  
  
  loc1 = info$sdimx
  loc2 = info$sdimy
  
  time = pseudotime_use+0.01
  datt = data.frame(time, loc1, loc2)
  datt2 = data.frame(start_x_dat, start_y_dat, end_x_dat, end_y_dat)
  p01=ggplot(datt, aes(x = loc1, y = loc2, color = time)) +
    geom_point( alpha = 1,size=pointsize) +
    scale_color_gradientn(colours = c("red", "green")) +
    theme_void()+
    theme(plot.title = element_text(size = textsize,  face = "bold"),
          text = element_text(size = textsize),
          legend.position = "bottom")
  p02= ggplot()+
    geom_segment(aes(x = start_x_dat, y = start_y_dat, xend = end_x_dat, yend = end_y_dat,colour = "black"),
                 arrow = arrow(length = unit(arrowlength,"cm")),size=arrowsize,color="black",data = datt2) +
    theme_void()+
    theme(plot.title = element_text(size = textsize,  face = "bold"),
          text = element_text(size = textsize),
          legend.position = "bottom")
  
  p03= ggplot()+
    geom_segment(aes(x = end_x_dat, y = end_y_dat, xend = start_x_dat, yend = start_y_dat,colour = "black"),
                 arrow = arrow(length = unit(arrowlength,"cm")),size=arrowsize,color="black",data = datt2) +
    theme_void()+
    theme(plot.title = element_text(size = textsize,  face = "bold"),
          text = element_text(size = textsize),
          legend.position = "bottom")
  
  time = pseudotime_use+0.1
  datt1 = data.frame(time, loc1, loc2)
  p1 = ggplot(datt1, aes(x = loc1, y = loc2)) +
    geom_point( alpha =1,size=pointsize,aes(color=clusterlabels)) +
    theme_void()+
    scale_colour_manual(values=color_in)+
    theme(plot.title = element_text(size = textsize,  face = "bold"),
          text = element_text(size = textsize),
          legend.position = "bottom")
  datt2 = data.frame(start_x_dat, start_y_dat, end_x_dat, end_y_dat)
  p2= geom_segment(aes(x = start_x_dat, y = start_y_dat, xend = end_x_dat, yend = end_y_dat,colour = "segment"),
                   arrow = arrow(length = unit(arrowlength,"cm")),size=arrowsize,arrow.fill="black",data = datt2)
  p22=p1+p2
  
  
  time = pseudotime_use+0.1
  datt1 = data.frame(time, loc1, loc2)
  p1 = ggplot(datt1, aes(x = loc1, y = loc2)) +
    geom_point( alpha =1,size=pointsize,aes(color=clusterlabels)) +
    theme_void()+
    scale_colour_manual(values=color_in)+
    theme(plot.title = element_text(size = textsize,  face = "bold"),
          text = element_text(size = textsize),
          legend.position = "bottom")
  datt2 = data.frame(start_x_dat, start_y_dat, end_x_dat, end_y_dat)
  p2= geom_segment(aes(x = end_x_dat, y = end_y_dat, xend = start_x_dat, yend = start_y_dat,colour = "segment"),
                   arrow = arrow(length = unit(arrowlength,"cm")),size=arrowsize,arrow.fill="black",data = datt2)
  p33=p1+p2
  
  
  return(list("Pseudotime"=p01,"Arrowplot1"=p02,"Arrowplot2"=p03,"Arrowoverlay1"=p22,"Arrowoverlay2"=p33))
  
}

strsplit_func = function(input){
  strsplit(input, split = "_")[[1]][1]
}

# get top gene loading for foreground gCPCA components
getTopGenes = function (loadings,whichComp,numGenes){
  refComponents <- whichComp
  nfea <- nrow(loadings)
  
  tempfea <- data.frame(matrix(ncol = length(refComponents), nrow = nfea))
  colnames(tempfea) <- refComponents
  
  genes <- data.frame(matrix(ncol = length(refComponents), nrow = nfea))
  colnames(genes) <- refComponents
  
  sign <- data.frame(matrix(ncol = length(refComponents), nrow = nfea))
  colnames(sign) <- refComponents
  
  for (i in 1:length(refComponents)) {
    ind <- sort(abs(loadings[,i]), decreasing = TRUE, index.return=TRUE)$ix
    tempfea[,i] <- loadings[ind[c(1:nfea)],i]
    genes[,i] <- rownames(loadings)[ind[c(1:nfea)]]
    for (j in 1:nrow(tempfea)){
      if(tempfea[j,i] > 0){
        sign[j,i] <- "+"
      }
      else{
        sign[j,i] <- "-"
      }
    }
  }
  
  ffea <- data.frame(matrix(ncol = 3, nrow = nfea))
  numComp <- length(whichComp)
  topG = list()
  compTopG <- c()
  
  for (g in 1:numComp){
    component <- g
    ffea[,1] <- genes[,component]
    ffea[,2] <- abs(tempfea[,component])
    ffea[,3] <- sign[,component]
    colnames(ffea) <- c("genes","weights","sign")
    rownames(ffea)  <- ffea$genes
    topG[[g]] <- ffea[c(1:numGenes),]
    compTopG[g] <- topG[[g]][["genes"]][1]
  }
  names(topG) <- paste(rep("component",numComp),c(1:numComp))
  res <- list(topG,compTopG)
  names(res) <- c("topGenes","compTopFirstGene")
  return(res)
}


# Plot top gene loading for foreground gCPCA components
require(ggpubr)
plotLoadings = function (geneList){
  pltTopG <- list()
  for (g in 1:length(geneList)){
    topGenes <- geneList[[g]]
    pltTopG[[g]] <- ggdotchart(topGenes, x = "genes", y = "weights",
                               sorting = "descending",
                               add = "segments",
                               add.params = list(color = "#999999", size = 4),
                               rotate = TRUE,
                               dot.size = 8,
                               label = topGenes$sign,
                               font.label = list(face = "bold", color = "white", size = 18, 
                                                 vjust = 0.4), 
                               ggtheme = theme_pubr()
    )+
      font("x.text", size = 18, color = "black", face = "bold.italic")+
      font("y.text", size = 18, color = "black", face = "bold.italic")+
      font("xy", size = 18, color = "black", face = "bold.italic")
    
  }
  names(pltTopG) <- paste(rep("component",length(geneList)),c(1:length(geneList)))
  return(pltTopG)
}

# perform GSEA of top genes
require(gprofiler2)
getPathways = function (geneList,bg){
  datt = list()
  for (g in 1:length(geneList)){
    tGenes <- geneList[[g]]
    gostres <- gost(query = tGenes$genes, 
                    organism = "mmusculus", ordered_query = FALSE, 
                    multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = FALSE, 
                    user_threshold = 0.05, correction_method = "fdr", 
                    domain_scope = "custom", custom_bg = bg, 
                    numeric_ns = "", 
                    sources = c("GO:MF","GO:BP", "GO:CC","KEGG"), 
                    as_short_link = FALSE)
    
    gostres$result = gostres$result[order(gostres$result$p_value),]
    gostres$result$Source = unlist(lapply(gostres$result$term_id,strsplit_func))
    datt[[g]] = gostres$result[1:10,]
    datt[[g]]$log10p = -log10(datt[[g]]$p_value)
  }
  names(datt) <- paste(rep("component",length(geneList)),c(1:length(geneList)))
  return(datt)
}

# plot GSEA results
plotPathways = function(genePathways){
  titleComp <- paste(rep("Component",length(genePathways)),c(1:length(genePathways)))
  pltGO=list()
  for(num in 1:length(genePathways)){
    pltGO[[num]]=ggplot(data=genePathways[[num]], aes(x=fct_reorder(tolower(term_name),log10p), y=log10p, 
                                                      fill=source,label = ifelse(significant ==TRUE, "*",""))) +
      geom_bar(stat="identity", position=position_dodge(),color="black",width=0.8)+
      #scale_fill_continuous(low='#F0E442', high='red', guide=guide_colorbar(reverse=TRUE)) + 
      scale_fill_manual(values = c("#FC4E07","#E7B800", "#56B4E9", "#009E73","#00AFBB","#F0E442", 
                                            "#0072B2", "#D55E00", "#CC79A7","#09622A"))+
                                              labs(title=paste0(titleComp[num]),x="Biological terms", y = "-log10(p value)")+
      coord_flip()+
      theme_classic()+
      #geom_text(vjust = 0.5) +
      geom_text(vjust = 1, nudge_y = 0.5)+
      #ylim(0,1)+
      theme(plot.title = element_text(size = 30,color="black",face="bold"),
            text = element_text(size = 30,color="black",face="bold"),
            #axis.title = element_text(size = 25,color="black",face="bold"),
            axis.text.x=element_text(size = 30,color="black",face="bold") ,
            legend.position = "right")# +
  }
  return(pltGO)
}
