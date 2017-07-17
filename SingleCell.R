
testGeneSets_sykes <- function (dataNorm, geneSet, sampleColors, keyword= "") {
  
  sd.data <- apply(dataNorm,1, sd) 
  dataNorm <- subset(dataNorm, sd.data > 0) ## keep those genes that have variance > 0
  
  pdf(paste0(keyword, "_checkpreDCPattern.pdf"))
  par(mar = c(5, 5, 4, 2) + 0.1)
  sc.MDS(dataNorm[rownames(dataNorm) %in% geneSet, ], colors = sampleColors[colnames(dataNorm)])
  legend("bottomleft", legend = c("Trm", "Tem", "bloodTem", "singleTcell"),
         pch = 16, col = c("red", "black", "deepskyblue", "blue"), cex=1)
  # sc.PlotGeneHeatmap(dataNorm, hk.genes = geneSet,
  #                    col.colors = cellType.color[colnames(dataNorm)], log = T, do.zscore = F)
  # sc.MDS(preDC.normData[rownames(dataNorm) %in% geneSet, ], colors = sampleColors[colnames(preDC.normData)] )
  dev.off()
  
  pca(dataNorm[rownames(dataNorm) %in% geneSet, ], colors = sampleColors[colnames(dataNorm)], bound = 1, prefix = keyword )
  
  #   distance.matrix <- 1 -cor(log10(dataNorm[preDC.bimodal_genes, ] + 1), method = "spearman")
  #   run_tsne_with_counts_distance(distance.matrix, keyword = "/Volumes/big/buffer/Kang/singleCell/tsne_result/corIRF8TFs",
  #                                 colors = batch.colors[colnames(dataNorm)], isDistance = TRUE)
  #   
  #   run_tsne_with_counts_distance(distance.matrix, keyword = "/Volumes/big/buffer/Kang/singleCell/tsne_result/corIRF8TFs",
  #                                 colors = cellType.color[colnames(dataNorm)], isDistance = TRUE)
  
  system(paste0("mkdir /Volumes/big/buffer/sykes/tsne_result/", keyword))
  distance.matrix <- 1 - cor(log10(dataNorm[rownames(dataNorm) %in% geneSet, ] + 1), method = "spearman")
  run_tsne_with_counts_distance(distance.matrix, keyword = paste0("/Volumes/big/buffer/sykes/tsne_result/", keyword, "/pipe_distance"),
                                colors = sampleColors[colnames(dataNorm)], isDistance = TRUE, legend = F)
  
  # run_tsne_with_counts_distance(dataNorm[rownames(dataNorm) %in% geneSet, ], keyword = paste0("/Volumes/big/buffer/Kang/singleCell/tsne_result/", keyword, "/pipe_counts"),
  #                               colors = cellType.color[colnames(dataNorm)], isDistance = FALSE)
  
}










check.twofold2 <- function(data){
  subdata <- abs(as.numeric(data[1:6]))
  max.data <- max(subdata)
  return(max.data > 1 )
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


MytTest <- function(m, n, data, outfile, cutoff = 0, log = T ) {
  
  c1 <- as.integer(m) + 1
  c2 <- as.integer(m) + as.integer(n)
  
  if(log){
    data <- log2(data + 1)
  }
  
  rowmeans <- apply(data, 1, mean)
  # baseMean <- numeric(length=nrow(data))
  meanA <- numeric(length=nrow(data))
  meanB <- numeric(length=nrow(data))
  # foldChange <- numeric(length=nrow(data))
  log2FoldChange <- numeric(length=nrow(data))
  pvalue <- numeric(length=nrow(data))
  
  outdata <- data.frame(rownames(data))
  for (i in 1:nrow(data)) {
    
    x <- data[i, 1:m]
    y <- data[i, c1:c2]
    
    # baseMean[i] <- mean(as.numeric(data[i,]))
    meanA[i] <- mean(as.numeric(x))
    meanB[i] <- mean(as.numeric(y))
    
    obj <- try(Ttest <- t.test(x, y), silent=TRUE)
    if (is(obj, "try-error")) 
      pvalue[i]<-1
    else
      pvalue[i] <- Ttest$p.value
    
  }
  
  # foldChange <- (meanA)/(meanB)
  # log2FoldChange <- log2(foldChange)
  log2FoldChange <- meanA - meanB
  
  
  bh = p.adjust(pvalue, method="BH")
  
  
  outdata[,"baseMean"] <- rowmeans
  outdata[,"meanA"] <- meanA
  outdata[,"meanB"] <- meanB
  # outdata[,"foldChange"] <- foldChange
  outdata[,"log2FoldChange"] <- log2FoldChange
  outdata[,"pvalue"] <- pvalue
  outdata[,"padj"] <- bh
  
  pfilter <- which(outdata[,"padj"] <= 0.05 & outdata$baseMean > cutoff & abs(outdata$log2FoldChange) > 1)
  outdata.sig <- outdata[pfilter,]
  sortedData <- outdata.sig[with(outdata.sig,order(padj)),]
  write.csv(outdata, paste0("allGenesList_", outfile), na="", quote = FALSE, col.names = T, row.names = T)
  write.csv(sortedData, outfile, na="", quote = FALSE, col.names = T,row.names = T)
}


## rank is small to big, the larger the x is, the larger percentile is
computePercentile <- function(x) trunc(rank(x))/length(x)

## make color transparent 
transparentColor <- function(color, alpha = 0.8){
  trans_colors <- color
  for(i in 1: length(trans_colors)){
    trans_colors[i] <- colorRampAlpha(trans_colors[i], n = 1, alpha = alpha)
  }
  return(trans_colors)
}


filterTFSimpleVersion <- function(TFFile, inputSize,  prefix){
  
  # TFFile <- "ILC1_NK_ChEA_2016_table.txt", the output file from enrichr/ChEA
  # inputSize <- 242, the number of DE genes 
  # prefix <- "ILC1", prefix of the files
  TFs_result  <- read.table(TFFile, header = T, sep = "\t")
  TFs_result <- subset(TFs_result, TFs_result$Adjusted.P.value < 0.05 & TFs_result$Combined.Score > 0) ## zscore should be at least negative, to the left side of the distribution
  
  ## compute fold change 
  overlapSta <- as.character(TFs_result$Overlap) 
  overlappingGeneNum <- as.numeric(matrix(unlist(strsplit(overlapSta, '/')), ncol = 2, byrow = T)[,1])
  TFTargetNum<- as.numeric((matrix(unlist(strsplit(overlapSta, '/')), ncol = 2, byrow = T)[,2]))
  
  #   TFTargetNum[is.na(TFTargetNum)] <- 1
  log2FC <- log2((overlappingGeneNum/inputSize)/(TFTargetNum / 20000))
  out <- data.frame(TFs_result[, 1:2], log2FC, TFs_result[, 3:4])
  
  out_file <- subset(out, abs(out$log2FC) > 1 & out$Adjusted.P.value < 0.05)
  
  TopTFs <- matrix(unlist(strsplit(as.character(out_file[, 1]), '_')), ncol = 5, byrow = T)[,1]
  
  
  write.csv(out_file, paste0(prefix, "_filteredTFResult.csv"))
  
  # return(data.frame(TopTFs, out)) ## TF result with FC computed
}



testGeneSets_breton <- function (dataNorm, cellColors = "blue",  shapes = 16, geneSet, keyword= "") {
  library(matrixStats)
  dataNorm <- subset(dataNorm, rowSds(as.matrix(dataNorm)) > 0)
  # dataNorm <- dataNorm[, apply(dataNorm[geneSet, ], 2, sum) > 0] ### remove samples that do not express any of the TFs
  # pdf(paste0(keyword, "_checkpreDCPattern.pdf"))
  # # sc.MDS(dataNorm[rownames(dataNorm) %in% geneSet, ], log = F)
  # neHeatmap(dataNorm, hk.genes = geneSet,
  #                    col.colors = cellColors, log = F, do.zscore = F)
  # # sc.MDS(preDC.normData[rownames(dataNorm) %in% geneSet, ], colors = cellType.color[colnames(preDC.normData)], log = F )
  # dev.off()
  # 
  # pca(dataNorm[rownames(dataNorm) %in% geneSet, ], colors = cellColors, bound = 1, prefix = keyword )
  
  #   distance.matrix <- 1 -cor(log10(dataNorm[preDC.bimodal_genes, ] + 1), method = "spearman")
  #   run_tsne_with_counts_distance(distance.matrix, keyword = "/Volumes/big/buffer/Kang/singleCell/tsne_result/corIRF8TFs",
  #                                 colors = batch.colors[colnames(dataNorm)], isDistance = TRUE)
  #   
  #   run_tsne_with_counts_distance(distance.matrix, keyword = "/Volumes/big/buffer/Kang/singleCell/tsne_result/corIRF8TFs",
  #                                 colors = cellType.color[colnames(dataNorm)], isDistance = TRUE)
  
  system(paste0("mkdir /Volumes/big/buffer/Kang/singleCell/tsne_result/", keyword))
  distance.matrix <- 1 - cor(dataNorm[rownames(dataNorm) %in% geneSet, ], method = "spearman")
  
  run_tsne_with_counts_distance(distance.matrix, colors = cellColors, shapes = shapes,  keyword = paste0("/Volumes/big/buffer/Kang/singleCell/tsne_result/", keyword, "/pipe_distance"),
                                isDistance = TRUE)
  
  # run_tsne_with_counts_distance(distance.matrix, keyword = paste0("/Volumes/big/buffer/Kang/singleCell/tsne_result/", keyword, "/pipe_counts"),
  #                                isDistance = FALSE)
}



post_process_ChEAResult <- function (filename) {
  ChEA_result <- read.csv(filename, header=F, stringsAsFactors = FALSE)
  colnames(ChEA_result) <- c("TFName","NuTargets","FracBackground","NuTargetsInput", "FracInput","FracDiffer", "Pvalue","Genes")
  head(ChEA_result)
  ChEA_result$p.adj <- p.adjust(ChEA_result$Pvalue, method="BH")
  ChEA_result$FoldChange <-  log2(ChEA_result$FracInput/ChEA_result$FracBackground)
  ChEA_result <- subset(ChEA_result, ChEA_result$Pvalue < 0.05 & ChEA_result$FoldChange > 0.5 )
  
  ## extract the genes name
  genes <- character(length=length(ChEA_result$TFName))
  for (i in 1:length(ChEA_result$TFName)){
    genes[i] <- unlist(strsplit(ChEA_result$TFName[i], "-"))[1]
  }
  ChEA_result <- ChEA_result[order(ChEA_result$Pvalue), ]
  write.csv(ChEA_result, paste0(filename, ".ChEA.sig.csv"))
  return(unique(as.character(genes)))
}

converHuman2mouse <- function(humanGenes, mapping = T) {
## basically convert the second to end character to lower case
  
  mouse_genes <- humanGenes
  if(mapping){
    mouse2human <- read.csv("~/Dropbox (CGC)/resources/mouse2human.csv", stringsAsFactors = F)
    mouse_genes <- mouse2human[,1][mouse2human[,2] %in% humanGenes]
  } else {
    for(i in 1:length(humanGenes)){
      temp <- as.character(humanGenes[i])
      mouse_genes[i] <- paste0(substr(temp, 1, 1 ), tolower(substr(temp, 2, nchar(temp))))
    }
  }
  return(mouse_genes)
}

convertMouse2Human <- function(mouseGenes, mapping = T) {
  ## basically convert the second to end character to lower case
  
  humanGenes <- mouseGenes 
  if(mapping){
    mouse2human <- read.csv("~/Dropbox (CGC)/resources/mouse2human.csv", stringsAsFactors = F)
    humanGenes <- mouse2human[,2][mouse2human[,1] %in% mouseGenes]
  } else {
    for(i in 1:length(mouseGenes)){
      humanGenes[i] <- toupper(mouseGenes[i])
    }
  }
  return(humanGenes)
}



vioplot_show_preDCSubV2 <- function(value_vec) {
  library(vioplot)
  par(cex.lab=2.3, cex.axis=2, mar = c(7.5, 7, 5, 5), srt = 45) 
  vioplot(value_vec[cells.right_preDC], value_vec[cells.left_preDC], 
          value_vec[cellType.color[colnames(solid.sc_data.norm)] == "red"],
          value_vec[cellType.color[colnames(solid.sc_data.norm)] == "blue"],
          names=c("", "", "", ""),
          col="lightcoral", border="black")
  
  vioplot( value_vec[cells.left_preDC],  at = 2,
           names="right",
           col="steelblue1", border="black", add=TRUE)
  vioplot( value_vec[cellType.color[colnames(solid.sc_data.norm)] == "red"],  at = 3,
           names="right",
           col="red", border="black", add=TRUE)
  
  vioplot( value_vec[cellType.color[colnames(solid.sc_data.norm)] == "blue"],  at = 4,
           names="right",
           col="blue", border="black", add=TRUE)
  text(x=1:4, y=par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]),##y=min(value_vec),
       labels=c("Pre-DC1", "Pre-DC2", "cDC1", "cDC2"), srt=45, adj=1, xpd=TRUE, cex = 2.3)
}


vioplot_show_preDCSub <- function(value_vec) {
  library(vioplot)
  par(cex.lab=2.3, cex.axis=2, mar = c(7.5, 7, 5, 5), srt = 45) 
  vioplot(value_vec[cells.left_preDC], value_vec[cells.right_preDC], 
          value_vec[cellType.color[colnames(solid.sc_data.norm)] == "red"],
          value_vec[cellType.color[colnames(solid.sc_data.norm)] == "blue"],
          names=c("", "", "", ""),
          col="steelblue1", border="black")
  
  vioplot( value_vec[cells.right_preDC],  at = 2,
           names="right",
           col="lightcoral", border="black", add=TRUE)
  vioplot( value_vec[cellType.color[colnames(solid.sc_data.norm)] == "red"],  at = 3,
           names="right",
           col="red", border="black", add=TRUE)
  
  vioplot( value_vec[cellType.color[colnames(solid.sc_data.norm)] == "blue"],  at = 4,
           names="right",
           col="blue", border="black", add=TRUE)
  text(x=1:4, y=par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]),##y=min(value_vec),
       labels=c("Pre-DC2", "Pre-DC1", "cDC2", "cDC1"), srt=45, adj=1, xpd=TRUE, cex = 2.3)
}


plotAQCHeatmap <- function(data,  facVec, labels = colnames(data), prefix="", margin=c(8,8)){
  ## example: plotAQCHeatmap(mouse_data, facVec = factor(mouse_condition$Trm), prefix = "mouse", margin=c(10,10), labels = with(mouse_condition, paste(celltype, Trm, batch, sep = ":")))
  
  
  condition <- data.frame(facVec)
  dds <- DESeqDataSetFromMatrix(countData = data, colData=condition, design = ~ facVec)
  tophat.rld <- rlog(dds)
  distsRL <- dist(t(assay(tophat.rld)))
  mat <- as.matrix(distsRL)
  rownames(mat) <- colnames(mat) <- labels
  pdf(paste0(prefix, "heatmap.pdf"))
  library(gplots)
  library(RColorBrewer)
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  heatmap.2(mat, trace="none", col=rev(hmcol), margin = margin, cexCol=0.8, cexRow = 0.8)
  dev.off()
}


heatmap_col <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                      "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                      "#4393C3", "#2166AC", "#053061")))(200)
heatmap_col2 <- c(colorRampPalette(c("#104E8B", "#1E90FF"))(50), colorRampPalette(c("#1E90FF",'#E5E5E5','#FF6347'))(30),colorRampPalette(c("#FF6347","#B22222"))(50))




home_cooked_colorgradient_onGene_Exp <- function (coloringGradient, backgroundIndex, 
                                                  shapes = 20, matrix, xlab="", ylab="", breaks = 30,
                                                  title = "", cex= 1.5, legend = TRUE) {
  #   coloringGradient <- log2(as.numeric(solid.sc_data.norm[gene, ]) + 1)
  library(plotrix)
  if ( sum(coloringGradient) < 0 ){
    color <- "blue"
  } else {
  my_palette <- colorRampAlpha(c( "blue","cyan", "yellow", "red"), n = breaks,  alpha = 0.8)
  color <- my_palette[as.numeric(cut(coloringGradient, breaks = breaks))]
  names(color) <- rownames(matrix) 
  
  solid_my_palette <- colorRampAlpha(c( "blue","cyan", "yellow", "red"), n = breaks,  alpha = 1)
  solid_color <- solid_my_palette[as.numeric(cut(coloringGradient, breaks = breaks))]
  names(solid_color) <- rownames(matrix) 
  
  }
#   plot(matrix, col = color, pch = shapes, cex.lab=1.25, cex.axis=1.25, xlab="t-SNE Dim 1", ylab="t-SNE Dim 2", cex=cex, main = title)
  
  
  color[backgroundIndex] <- colorRampAlpha("gray", n = 1, alpha = 0.7)
  plot(matrix, col = color, pch = shapes, cex.lab=1.25, cex.axis=1.25, xlab=xlab, ylab=ylab, cex=cex, main = title, asp =1)
  
  solid_color[backgroundIndex] <- colorRampAlpha("gray", n = 1, alpha = 0.7)
  points(matrix, col = solid_color, pch = shapes - 15, cex.lab=1.25, cex.axis=1.25, xlab=xlab, ylab=ylab, cex=cex, main = title, asp =1)
  
  
  x_length <- 10
  y_length <-  2
  left_x <- max(matrix[,1] - 11)
  left_y <- min(matrix[,2] - 0.7)
  #   text(left_x + 8, left_y+6, labels = "CD45+ output in log10", cex = 0.8 )
  if(legend){
  color.legend(left_x,left_y,left_x + x_length ,left_y + y_length , c("min", "max"), my_palette, align="rb",gradient="x")
  }
}




home_cooked_colorgradient_onGene_Exp_breton <- function (matrix, coloringGradient, backgroundIndex, 
                                                  shapes = 20,  xlab="", ylab="", breaks = 50,
                                                  cols.use = c( "blue","cyan", "yellow", "red"),
                                                  title = "", cex= 1.5, legend = TRUE, legend.pos = "left") {
  #   coloringGradient <- log2(as.numeric(solid.sc_data.norm[gene, ]) + 1)
  library(plotrix)
  if (abs(sum(coloringGradient)) == sum(backgroundIndex) ) {
    print("gradient vector is 0")
    color <- rep("blue", nrow(matrix))
    names(color) <- rownames(matrix)
    solid_color <- color
    my_palette <- rep("blue", length = breaks)
  } else {
    my_palette <- colorRampAlpha(cols.use, n = breaks,  alpha = 0.8)
    # my_palette <- colorRampAlpha(c( "#EF00F2","black", "#F9FE0D"), n = breaks,  alpha = 0.8)
    color <- my_palette[as.numeric(cut(coloringGradient, breaks = breaks))]
    names(color) <- rownames(matrix) 
    
    solid_my_palette <- colorRampAlpha(cols.use, n = breaks,  alpha = 1)
    # solid_my_palette <- colorRampAlpha(c( "#EF00F2","black", "#F9FE0D"), n = breaks,  alpha = 1)
    solid_color <- solid_my_palette[as.numeric(cut(coloringGradient, breaks = breaks))]
    names(solid_color) <- rownames(matrix) 
    
  }
  #   plot(matrix, col = color, pch = shapes, cex.lab=1.25, cex.axis=1.25, xlab="t-SNE Dim 1", ylab="t-SNE Dim 2", cex=cex, main = title)
  
  
  color[backgroundIndex] <- colorRampAlpha("gray", n = 1, alpha = 0.7)
  plot(matrix, col = color, pch = shapes, cex.lab=1.25, cex.axis=1.25, xlab=xlab, ylab=ylab, cex=cex, main = title, asp =1)
  points(matrix[!backgroundIndex, ], col = color[!backgroundIndex], pch = shapes[!backgroundIndex], cex.lab=1.25, cex.axis=1.25, xlab=xlab, ylab=ylab, cex=cex, main = title, asp =1)
  
  solid_color[backgroundIndex] <- colorRampAlpha("gray", n = 1, alpha = 0.7)
  points(matrix[!backgroundIndex, ], col = solid_color[!backgroundIndex], pch = shapes[!backgroundIndex] - 15, cex.lab=1.25, cex.axis=1.25, xlab=xlab, ylab=ylab, cex=cex, main = title, asp =1)
  
  
  x_length <- 10
  y_length <-  2
  # left_x <- max(matrix[,1] -12)
  if(legend.pos == "left"){
    left_x <- min(matrix[, 1] + 2)
  } else {
    left_x <- max(matrix[, 1] - 12)
  }
  left_y <- min(matrix[, 2] -3 )
  #   text(left_x + 8, left_y+6, labels = "CD45+ output in log10", cex = 0.8 )
  if(legend){
    color.legend(left_x, left_y, left_x + x_length, left_y + y_length, c("min", "max"), my_palette, align="rb", cex = 1, gradient="x")
  }
}


 ## plot IRF8/IRF4 gradient on the tSNE plot
showGeneRatioInTsne <- function(gene1, gene2, exp.data, shapes, log = T, 
                                cell_type_color, legend = T, cols.use = c("blue", "black", "red"),
                                cex = 1.5, 
                                keyword = "", legend.pos = "left", plotLine= F, ...) {
  pdf(paste0(keyword, "_", gene1, "Over", gene2, "_Tsne.pdf"), height = 4.5, width = 13.5)
  par(mfrow=c(1, 3))
  for (cellColor in c("orange", "blue", "red")) {
    if (cellColor == "orange") {
      cell <- "preDC"
    } else if (cellColor == "blue") {
      cell <- "CD1c"
    } else if (cellColor == "red" ) {
      cell <- "CD141"
    }
    
    cex_vec <- rep(cex, length(cell_type_color))
    names(cex_vec) <- names(cell_type_color)
    
    if (log) {
      tempColors <- log2(as.numeric(exp.data["IRF8", ]) + 1) - log2(as.numeric(exp.data["IRF4", ]) + 1)
    } else {
      tempColors <- as.numeric(exp.data["IRF8", ]) - as.numeric(exp.data["IRF4", ])
    }
    
    names(tempColors) <- colnames(exp.data)
    backgroundIndex <- names(cell_type_color)[ cell_type_color!= cellColor & names(cell_type_color) %in% colnames(exp.data)]
    
    backgroundIndex <- cell_type_color!= cellColor & names(cell_type_color) %in% colnames(exp.data)
    
    # tempColors[backgroundIndex] <- -1
    cex_vec[backgroundIndex] <- cex - 0.2
    home_cooked_colorgradient_onGene_Exp_breton(coloringGradient = tempColors[colnames(exp.data)],  backgroundIndex = backgroundIndex,
                                                legend = legend, cols.use = cols.use,  legend.pos = legend.pos, 
                                                shapes = shapes[colnames(exp.data)], matrix = tsne_result, title = "", cex = cex_vec[colnames(exp.data)], ...)
    if (plotLine) {
      abline(v = 0.6, lty = "dotted")
    }
  }
  dev.off()
}

showMarkerInTsne <- function(gene, exp.data, shapes, log = F, 
                             cell_type_color, legend = T, 
                             cols.use = c( "blue","cyan", "yellow", "red"), keyword = "", 
                             legend.pos = "left",
                             cex = 1.8, 
                             plotLine = F) {
  pdf(paste0(keyword, "_", gene,".pdf"), height = 4.5, width = 13.5)
  par(mfrow=c(1, 3))
  for(cellColor in c("orange", "blue", "red")){
    if(cellColor == "orange"){
      cell <- "preDC"
    } else if (cellColor == "blue") {
      cell <- "CD1c"
    } else if (cellColor == "red" ) {
      cell <- "CD141"
    }
    cex_vec <- rep(cex, length(cell_type_color))
    names(cex_vec) <- names(cell_type_color)
    if (log) {
      tempColors <- log2(as.numeric(exp.data[gene, ])+1)
    } else {
      tempColors <- as.numeric(exp.data[gene, ])
    }
   
    names(tempColors) <- colnames(exp.data)
    backgroundIndex <- names(cell_type_color)[ cell_type_color!= cellColor & names(cell_type_color) %in% colnames(exp.data)]
    
    backgroundIndex <- cell_type_color!= cellColor & names(cell_type_color) %in% colnames(exp.data)
    
    # tempColors[backgroundIndex] <- -1
    cex_vec[backgroundIndex] <- cex - 0.2
    home_cooked_colorgradient_onGene_Exp_breton(coloringGradient = tempColors[colnames(exp.data)],  backgroundIndex = backgroundIndex,
                                                legend = legend, cols.use = cols.use, legend.pos = legend.pos,
                                                shapes = shapes[colnames(exp.data)], matrix = tsne_result, title = "",
                                                cex = cex_vec[colnames(exp.data)])
    if (plotLine) {
      abline(v = 0.6, lty = "dotted")
    }
  
    }
  dev.off()
}



Get.gene.id.to.name.dict <- function(){
  gene.name.dict <- read.table("/Users/wm2313/Dropbox (CGC)/SingleCell_DC/DCAnalysis/Data/humanEnsemblandERCC92_geneid_genename_dictionary.txt")
  gene.names <- gene.name.dict$V2
  names(gene.names) <- gene.name.dict$V1
  return(gene.names)
}

Get.gene.name.to.id.dict <- function(){
  gene.name.dict <- read.table("/Users/wm2313/Dropbox (CGC)/SingleCell_DC/DCAnalysis/Data/humanEnsemblandERCC92_geneid_genename_dictionary.txt")
  gene.names <- gene.name.dict$V1
  names(gene.names) <- gene.name.dict$V2
  return(gene.names)
}


run_tsne_with_counts_distance <- function (matrix, keyword, shapes = 16, colors = "blue", isDistance = FALSE, legend = T) {
  library(Rtsne)
  set.seed(42) # Set a seed if you want reproducible results
  
  # for(perplexity in seq(10, 22, by = 2)){
  for(perplexity in seq(20, 40, by = 5)){
    for(i in 1:4){
      for(theta in seq(0, 0.1, by = 0.05)){
        data.file.name <- paste0(keyword, '_','perplexity=', perplexity,'_try', i , '_theta=', theta, '_tsne_map.csv')
        
        tsne_out <- Rtsne(matrix, perplexity = perplexity, theta = theta, is.distance = isDistance, check_duplicates = FALSE) # Run TSNE
        
        # Show the objects in the 2D tsne representation
        pdf(paste0(keyword, '_','perplexity=', perplexity,'_try', i , '_theta=', theta, '_tsne_map.pdf'))
#         plot(tsne_out$Y,col=colors, pch = 16)
#         legend("bottomright", legend = c("Batch 1 (DCs)", "Batch 2 (preDCs)", "Batch 3 (mixed)"), pch = 16, col = c("orange", "brown", "black"), cex=1)
# #         dev.off()
        write.csv(tsne_out$Y, paste0(keyword, '_','perplexity=', perplexity,'_try', i , '_theta=', theta, '_tsne_map.csv'))
        
        map <- read.csv(data.file.name, row.names = 1)
#         pdf(paste0(keyword, '_','perplexity=', perplexity,'_try', i , '_theta=', theta, '_tsne_map_cellType.pdf'))
        plot(map,col=colors, cex.lab=1.25, cex.axis=1.25, xlab = "t-SNE dim1", ylab="t-SNE dim1", pch = shapes)
        if(legend) {
        legend('bottomright', pch = 16, col = c('red', 'blue', 'orange'),legend = c("CD141", "CD1c", "PreDC"))
        }
        dev.off()
      }
    }
  }
}



testGeneSets <- function (dataNorm, geneSet, keyword= "") {
  
  sd.data <- apply(dataNorm,1, sd) 
  dataNorm <- subset(dataNorm, sd.data > 0) ## keep those genes that have variance > 0
  
  pdf(paste0(keyword, "_checkpreDCPattern.pdf"))
  par(mar = c(5, 5, 4, 2) + 0.1)
  sc.MDS(dataNorm[rownames(dataNorm) %in% geneSet, ], colors = cellType.color[colnames(dataNorm)])
  sc.PlotGeneHeatmap(dataNorm, hk.genes = geneSet,
                     col.colors = cellType.color[colnames(dataNorm)], log = T, do.zscore = F)
  sc.MDS(preDC.normData[rownames(dataNorm) %in% geneSet, ], colors = cellType.color[colnames(preDC.normData)] )
  dev.off()
  
  pca(dataNorm[rownames(dataNorm) %in% geneSet, ], colors = cellType.color[colnames(dataNorm)], bound = 1, prefix = keyword )
  
  #   distance.matrix <- 1 -cor(log10(dataNorm[preDC.bimodal_genes, ] + 1), method = "spearman")
  #   run_tsne_with_counts_distance(distance.matrix, keyword = "/Volumes/big/buffer/Kang/singleCell/tsne_result/corIRF8TFs",
  #                                 colors = batch.colors[colnames(dataNorm)], isDistance = TRUE)
  #   
  #   run_tsne_with_counts_distance(distance.matrix, keyword = "/Volumes/big/buffer/Kang/singleCell/tsne_result/corIRF8TFs",
  #                                 colors = cellType.color[colnames(dataNorm)], isDistance = TRUE)
  
  system(paste0("mkdir /Volumes/big/buffer/Kang/singleCell/tsne_result/", keyword))
  distance.matrix <- 1 - cor(log10(dataNorm[rownames(dataNorm) %in% geneSet, ] + 1), method = "spearman")
  run_tsne_with_counts_distance(distance.matrix, keyword = paste0("/Volumes/big/buffer/Kang/singleCell/tsne_result/", keyword, "/pipe_distance"),
                                colors = cellType.color[colnames(dataNorm)], isDistance = TRUE)
  
  # run_tsne_with_counts_distance(dataNorm[rownames(dataNorm) %in% geneSet, ], keyword = paste0("/Volumes/big/buffer/Kang/singleCell/tsne_result/", keyword, "/pipe_counts"),
  #                               colors = cellType.color[colnames(dataNorm)], isDistance = FALSE)
  
}

testGeneSets_preDC <- function (dataNorm, geneSet, color, keyword= "") {
  
  sd.data <- apply(dataNorm,1, sd) 
  dataNorm <- subset(dataNorm, sd.data > 0) ## keep those genes that have variance > 0
  
  pdf(paste0(keyword, "_checkpreDCPattern.pdf"))
  par(mar = c(5, 5, 4, 2) + 0.1)
  sc.MDS(dataNorm[rownames(dataNorm) %in% geneSet, ], colors = color[colnames(dataNorm)])
  sc.PlotGeneHeatmap(dataNorm, hk.genes = geneSet,
                     col.colors = color[colnames(dataNorm)], log = T, do.zscore = F)
  dev.off()
  
  
  #   distance.matrix <- 1 -cor(log10(dataNorm[preDC.bimodal_genes, ] + 1), method = "spearman")
  #   run_tsne_with_counts_distance(distance.matrix, keyword = "/Volumes/big/buffer/Kang/singleCell/tsne_result/corIRF8TFs",
  #                                 colors = batch.colors[colnames(dataNorm)], isDistance = TRUE)
  #   
  #   run_tsne_with_counts_distance(distance.matrix, keyword = "/Volumes/big/buffer/Kang/singleCell/tsne_result/corIRF8TFs",
  #                                 colors = cellType.color[colnames(dataNorm)], isDistance = TRUE)
  
  system(paste0("mkdir /Volumes/big/buffer/Kang/singleCell/tsne_result/", keyword))
  distance.matrix <- 1 - cor(log10(dataNorm[rownames(dataNorm) %in% geneSet, ] + 1), method = "spearman")
  run_tsne_with_counts_distance(distance.matrix, keyword = paste0("/Volumes/big/buffer/Kang/singleCell/tsne_result/", keyword, "/pipe_distance"),
                                colors = color[colnames(dataNorm)], isDistance = TRUE, legend = F)
  
  # run_tsne_with_counts_distance(dataNorm[rownames(dataNorm) %in% geneSet, ], keyword = paste0("/Volumes/big/buffer/Kang/singleCell/tsne_result/", keyword, "/pipe_counts"),
  #                               colors = cellType.color[colnames(dataNorm)], isDistance = FALSE)
  
}


## calculate the differentially expressed genes between two groups of single cells
scdediff <- function(model, counts, group1, group2, prefix, n.cores=4){
  scdediff.file <- paste(prefix, ".scdediff.Rdata", sep="")
  if (!file.exists(scdediff.file)){
    groups <- c(group1, group2)
    names(groups) <- groups
    groups[group1] <- "group1"
    groups[group2] <- "group2"
    groups <- as.factor(groups)
    model <- model[names(groups), ]
    counts <- counts[, names(groups)]
    counts <- counts[rowSums(counts)>0, ]
    spikein.rows <- grep("ERCC-", row.names(counts))
    if (length(c(spikein.rows))>0){
      counts <- counts[-c(spikein.rows), ]
    }
    library(scde)
    o.prior <- scde.expression.prior(models=model, counts=counts, 
                                     length.out=400, show.plot=F)
    ediff <- scde.expression.difference(model, counts, o.prior, 
                                        groups=groups[colnames(counts)], n.randomizations=100, 
                                        n.cores=n.cores, verbose=1)
    save(ediff, file=scdediff.file)
    system("say I am a monkey")
  }
  load(scdediff.file)
  return(ediff[order(abs(ediff$Z), decreasing=T), ])
}






filterTF <- function(DEsigfile, DEfile, TFFile, tissue = "sp", prefix){
  ## filter ChEA TFs. CD69hi VS lo in Naive
  DE.re.sig <- read.csv(DEsigfile, row.names = 1)
  TFs_result  <- read.table(TFFile, header = T, sep = "\t")
  TFs_result <- subset(TFs_result, TFs_result$Adjusted.P.value < 0.05 & TFs_result$Combined.Score > 0) ## zscore should be at least negative, to the left side of the distribution
  
  ## compute fold change 
  inputSize <- as.numeric(nrow(DE.re.sig))
  overlapSta <- as.character(TFs_result$Overlap) 
  overlappingGeneNum <- as.numeric(matrix(unlist(strsplit(overlapSta, '/')), ncol = 2, byrow = T)[,1])
  TFTargetNum<- as.numeric((matrix(unlist(strsplit(overlapSta, '/')), ncol = 2, byrow = T)[,2]))
  
#   TFTargetNum[is.na(TFTargetNum)] <- 1
  log2FC <- log2((overlappingGeneNum/inputSize)/(TFTargetNum / 20000))
  out <- data.frame(TFs_result[, 1:2], log2FC, TFs_result[, 3:7])
  
  TopTFs <- matrix(unlist(strsplit(as.character(TFs_result[, 1]), '_')), ncol = 5, byrow = T)[,1]
  
  DE.re <- read.csv(DEfile, row.names = 1)
  DE.re <- DE.re[DE.re$pvalue < 0.05 & abs(DE.re$log2FoldChange) > 1, ]
  
  TFs.index <- which(TopTFs %in% rownames(DE.re))
  cellTypeInfo <- TFs_result[TFs.index, 1]
  
  usefulTF <- rownames(DE.re[TopTFs[TopTFs %in% rownames(DE.re)], ])
  
#   pdf(paste0(prefix,"_",tissue, "_TFHeatmap.pdf"))
#   sc.PlotGeneHeatmap( data.norm[, grepl(prefix, colnames(data.norm)) & grepl(tissue, colnames(data.norm))], hk.genes = unique(usefulTF), margins = c(12, 12),  log = T, do.zscore = T)
#   dev.off()
  
  filteredTF <- data.frame(DE.re[TopTFs[TopTFs %in% rownames(DE.re)], ], out[TFs.index, ])
  filteredTF <- filteredTF[order(rownames(filteredTF)), ]
  write.csv(filteredTF, paste0(prefix,"_",tissue, "_filteredTFResult.csv"))
  
  return(data.frame(TopTFs, out)) ## TF result with FC computed
}


CompueteChEAFoldChange <- function(chea.file, inputSize){
  TFs_result  <- read.table(chea.file, header = T, sep = "\t")
  #   TFs_result <- subset(TFs_result, TFs_result$Adjusted.P.value < 0.05)
  
  ## compute fold change 
  overlapSta <- as.character(TFs_result$Overlap) 
  overlappingGeneNum <- as.numeric(matrix(unlist(strsplit(overlapSta, '/')), ncol = 2, byrow = T)[,1])
  TFTargetNum<- as.numeric((matrix(unlist(strsplit(overlapSta, '/')), ncol = 2, byrow = T)[,2]))
  
  #   TFTargetNum[is.na(TFTargetNum)] <- 1
  
  log2FC <- log2((overlappingGeneNum/inputSize)/(TFTargetNum / 20000))
  out <- data.frame(TFs_result[, 1:2], log2FC, TFs_result[, 3:7])
  
  TFNames <- matrix(unlist(strsplit(as.character(TFs_result[, 1]), '_')), ncol = 5, byrow = T)[,1]
  
  write.csv(out, paste0(chea.file, "FCComputed.csv" ))
  return(data.frame(TFNames, out))
}

CompueteChEAFoldChange_singleCell <- function(chea.file, inputSize, keyword){
  database_size <- 48230
  TFs_result  <- read.table(chea.file, header = T, sep = "\t")
  #   TFs_result <- subset(TFs_result, TFs_result$Adjusted.P.value < 0.05)
  
  ## compute fold change 
  overlapSta <- as.character(TFs_result$Overlap) 
  overlappingGeneNum <- as.numeric(matrix(unlist(strsplit(overlapSta, '/')), ncol = 2, byrow = T)[,1])
  TFTargetNum<- as.numeric((matrix(unlist(strsplit(overlapSta, '/')), ncol = 2, byrow = T)[,2]))
  
  #   TFTargetNum[is.na(TFTargetNum)] <- 1
  
  log2FC <- log2((overlappingGeneNum/inputSize)/(TFTargetNum / database_size))
  out <- data.frame(TFs_result[, 1:2], log2FC, TFs_result[, 3:7])
  out$Overlap <- as.character(out$Overlap)
  
  TFNames <- matrix(unlist(strsplit(as.character(TFs_result[, 1]), '_')), ncol = 5, byrow = T)[,1]
  
  # TFNames <- TFNames %in% select ## make sure the TF is expressed in at least 10 cells of the dataset
  
  select_exp <- rownames(solid.sc_data.norm)[apply(solid.sc_data.norm[, grepl(keyword, cell.info[colnames(solid.sc_data.norm), "TYPE"])] > 0, 1, sum) > 10]
  write.csv(out, paste0(chea.file, "_FCComputed.csv" ), row.names = F)
  out.sig <- subset(out, out$log2FC > 2 & out$Adjusted.P.value < 0.05 & TFNames %in% select_exp & out$Combined.Score > 0 & out$Z.score < 0)
  out.sig <- out.sig[order(out.sig$log2FC, decreasing = T),]
  cat("# of enriched TFs ", nrow(out.sig))
  write.csv(out.sig, paste0(chea.file, "_FCComputed_sig.csv" ), row.names = F)
  
  
  # return(nrow(out.sig))
}


filterTF_DC <- function(DEsigfile, DEfile, TFFile, prefix){
  DE.re.sig <- read.csv(DEsigfile, row.names = 1)
  TFs_result  <- read.table(TFFile, header = T, sep = "\t")
  ## first filter by p.adj
  TFs_result <- subset(TFs_result, TFs_result$Adjusted.P.value < 0.05 & TFs_result$Combined.Score > 0) ## zscore should be at least negative, to the left side of the distribution
  
  ## compute fold change 
  inputSize <- as.numeric(nrow(DE.re.sig))
  database_size <- 48230
  overlapSta <- as.character(TFs_result$Overlap) 
  overlappingGeneNum <- as.numeric(matrix(unlist(strsplit(overlapSta, '/')), ncol = 2, byrow = T)[,1])
  TFTargetNum<- as.numeric((matrix(unlist(strsplit(overlapSta, '/')), ncol = 2, byrow = T)[,2]))
  
  log2FC <- log2((overlappingGeneNum / inputSize) / (TFTargetNum / database_size))
  out <- data.frame(TFs_result[, 1:2], log2FC, TFs_result[, 3:7])
  
  ## filter by fold change
  out <- subset(out, out$log2FC > 0 & out$Adjusted.P.value < 0.05)
  out <- out[order(out$log2FC, decreasing = T),]
  
  
  TopTFs <- matrix(unlist(strsplit(as.character(out[, 1]), '_')), ncol = 5, byrow = T)[,1]
  
  DE.re <- read.csv(DEfile, row.names = 1)
  DE.re <- DE.re[DE.re$pvalue < 0.05 & abs(DE.re$log2FoldChange) > 1, ]
  
  TFs.index <- which(TopTFs %in% rownames(DE.re))
  cellTypeInfo <- out[TFs.index, 1]
  
  filteredTF <- data.frame(DE.re[TopTFs[TopTFs %in% rownames(DE.re)], ], out[TFs.index, ])
  filteredTF <- filteredTF[order(rownames(filteredTF)), ]
  write.csv(filteredTF, paste0(prefix, "_filteredTFResult.csv"))
  
  final_TFs <- unique(TopTFs[TopTFs %in% rownames(DE.re)])
  
  return(final_TFs) ## TF result with FC computed
}







drawPairwiseHeatmap <- function(prefix = "JSM", normed.data, condition){
#   prefix <- "JSM"
  if(prefix == "euclidean"){
    sampleDists <- dist(t(normed.data))
  } else if(prefix == "spearman"){ 
    sampleDists <- as.matrix(1-cor(normed.data, method="spearman"))
  } else if(prefix == "JSM"){ 
    sampleDists <- JSD_distance_compute(input = normed.data, cut = 30, method = "JSM")
  }
  rownames(sampleDists) <- colnames(sampleDists) <- with(condition, paste(rownames(condition), cellType, CD69, tissue, sep=" : "))
  pdf(paste0(prefix, "heatmap.pdf"))
  library(gplots)
  library(RColorBrewer)
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  heatmap.2(sampleDists, trace="none", col=rev(hmcol), margin=c(13,13), cexCol=0.8)
  dev.off()
}


sc.PlotGeneHeatmap <- function(norm.counts, hk.genes, my_palette = NULL,
                               Rowv=TRUE, Colv=TRUE, gene.labels=NULL, 
                               cexRow=0.75, cexCol=0.75, binary=FALSE, 
                               colOrder=NULL, rowOrder=NULL, exp.cutoff=0, col.colors=NULL, row.colors = NULL, key.xlab = "value",
                               log = FALSE, do.zscore=FALSE, cluster = F, col.clustersNum = 3,  row.clustersNum = 3,
                               method = "complete", distance = "euclidean", margins = c(7, 5), prefix = ""){
  
  ## compute z score for a vector
  z.score <- function(x){
    
    mean <- mean(x)
    sd <- sd(x)
    zscore <- (x-mean)/sd
    
    return(zscore)
  }
  
  
  hclustfunc <- function(x){
    
    library(cba)
    hc <-  hclust(x, method = method) ## method is hierachical method, for example, complete, average, ward.D
    co <- order.optimal(x, hc$merge) ## optimize the clutering result by shuffiling the leaves while keeping the relative structure.
    
    hc$merge <- co$merge
    hc$order <- co$order 
    
    return(hc)
  }
  
  if (distance == "pearson") {
    ## distfunc <- function(x) as.dist((1 - cor(t(x), method = "spearman"))/2) 
    distfunc <- function(x) as.dist(1 - (cor(t(x), method = "pearson"))) 
    ##The standard call to dist computes the distance between the rows of the matrix provided,
    #     cor computes the correlation between columns of the provided matrix, 
    #     so the above example to work, you need to transpose the matrix: 
  }  else if (distance == "spearman") {
    distfunc <- function(x) as.dist((1 - cor(t(x), method = "spearman"))/2)
    # distfunc <- function(x) as.dist(1 - (cor(t(x), method = "spearman"))) 
    ##The standard call to dist computes the distance between the rows of the matrix provided,
    #     cor computes the correlation between columns of the provided matrix, 
    #     so the above example to work, you need to transpose the matrix: 
  } else {
    distfunc <- function(x) dist(x, method = distance)
  }
  
  ## subset the data with user-defined gene set, and keep the order of the input genes
  cat("CandidateGeneHeatmap\n")
  overlap_genes <- intersect(rownames(norm.counts), hk.genes)
  subset_genes <- hk.genes[hk.genes %in% overlap_genes]
  
  if (length(overlap_genes) < length(hk.genes)) {
    cat(paste0("reduced genes, # overlap genes", length(overlap_genes), " # of input genes: ", length(hk.genes), "\n"))
  }
  
  gene.counts <- norm.counts[subset_genes, ] ## use this keep the original order of genes 
  # sd.data <- apply(gene.counts, 1, sd)
  # gene.counts <- subset(gene.counts, sd.data>0) ## keep those genes that have variance > 0
  valid_genes <- rownames(gene.counts)
  
  if(log){
    gene.counts <- log10(gene.counts + 1)
  }
  if(do.zscore){
    z <- apply(gene.counts, 1, z.score)
    gene.counts <- t(z)
  }
  if (binary){
    gene.counts <- gene.counts > exp.cutoff
    gene.counts <- gene.counts + 0
  }
  if (!is.null(gene.labels)){
    new_labRow <- gene.labels 
  } else {
    new_labRow <- rownames(gene.counts)
  }
  
  if(is.null(my_palette))
    my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)
  
  library(gplots)
  gene.counts <- as.matrix(gene.counts)
  
  if(!is.null(colOrder)){
    gene.counts <- gene.counts[, colOrder]
    col.colors <- col.colors[colOrder]
    print(paste0("new column order: ", colnames(gene.counts)))
  }
  
  if(!is.null(rowOrder)){
    gene.counts <- gene.counts[rowOrder,]
    print(paste0("new row order: ", rowOrder))
  }
  
  print(paste0("# rows of data matrix: ", nrow(gene.counts)))
  print(paste0("# length of labrow: ", length(new_labRow)))
  
  if (cluster) {
    # obtain the row cluster result, usually it is about genes
    fit_row <- hclustfunc(distfunc(gene.counts))
    row_clusters <- cutree(fit_row, row.clustersNum)
    original_cluster_num <- row_clusters
    colorCandi <- c("orchid", rgb(246,190,84,max = 255), rgb(117,199,198,max = 255), rgb(94,102,176,max = 255))
    # colorCandi <- c("purple", "blue", "red", "green")
    if (row.clustersNum > length(colorCandi)) {
      colorCandi[5:length(row.clustersNum)] <- row.clustersNum[5:length(row.clustersNum)]
    }
    for (i in 1 : row.clustersNum) {
      row_clusters[row_clusters == i] <- colorCandi[i]
    }
    out <- data.frame(rownames(gene.counts), as.character(original_cluster_num), as.character(row_clusters))
    colnames(out) <- c("genes", "clusterNum", "clusterColor")
    write.csv(out, paste0(prefix, "_row_cluster_result.csv"), row.names = F)
    
    
    ## obtain the column cluster result, usually it is about samples
    fit_column <- hclustfunc(distfunc(t(gene.counts)))
    col_clusters <- cutree(fit_column, col.clustersNum) 
    col_clusters <- as.character(col_clusters)
    col_clusters[col_clusters == 2] <- "orange"
    # colorCandi <- c("gray10", "gray84", "gray59", "gray43")
    # if (clustersNum > length(colorCandi)){
    #   colorCandi[5:length(clustersNum)] <- clustersNum[5:length(clustersNum)]
    # }
    # for(i in 1 : clustersNum){
    #   clusters[clusters == i] <- colorCandi[i]
    # }
    out <- data.frame(colnames(gene.counts), as.character(col_clusters))
    colnames(out) <- c("samples", "clusterNum")
    write.csv(out, paste0(prefix, "_column_cluster_result.csv"), row.names = F)
    
    
    
    if (is.null(col.colors)) { 
      ## no specified color for the columns
      if (!is.null(row.colors)) { 
        ## specified color for the rows
        hm <- heatmap.2(gene.counts, scale="none", labRow = new_labRow, colRow = as.character(row_clusters),
                        hclust = hclustfunc, distfun=distfunc,
                        density.info = 'none', key.title=NA,
                        trace="none",
                        cexRow=cexRow, cexCol = cexCol,
                        Rowv=Rowv, Colv=Colv,
                        #lhei = c(2, 8),
                        col=my_palette, ColSideColors=as.character(col_clusters), 
                        RowSideColors = row.colors[valid_genes],
                        margins = margins)
        return(hm)
      } else {
        hm <- heatmap.2(gene.counts, scale="none", labRow = new_labRow, colRow = as.character(row_clusters),
                        hclust=hclustfunc, distfun=distfunc,
                        density.info = 'none', key.title=NA,
                        trace="none",
                        cexRow=cexRow, cexCol = cexCol,
                        Rowv=Rowv, Colv=Colv,
                        #lhei = c(2, 8),
                        col=my_palette, ColSideColors=as.character(col_clusters), 
                        RowSideColors = as.character(row_clusters), 
                        margins = margins)
        return(hm)
      }
    } else { ## specified column colors
      hm <- heatmap.2(gene.counts, scale="none", labRow = new_labRow, colRow = as.character(row_clusters),
                      hclust=hclustfunc, distfun=distfunc,
                      density.info = 'none', key.title=NA,
                      RowSideColors = row_clusters,
                      trace="none",
                      cexRow=cexRow, cexCol=cexCol,
                      Rowv=Rowv, Colv=Colv, 
                      ##lhei = c(2, 8),
                      col=my_palette, ColSideColors = col.colors, 
                      # RowSideColors = row.colors,
                      # [valid_genes], 
                      margins = margins)
    }
    return(hm)
  } else { ## cluster result end 
    ## no clustering required 
    if (!is.null(col.colors)) {
      if(!is.null(row.colors)){
        hm <- heatmap.2(gene.counts, scale="none", col=my_palette, 
                        labRow = new_labRow, 
                        RowSideColors = row.colors[valid_genes],
                        ColSideColors = col.colors,
                        keysize = 1,
                        key.title = NA,
                        density.info = 'none',
                        key.xlab = key.xlab,
                        dendrogram = 'both',
                        #lhei = c(2, 8),
                        trace="none", cexRow=cexRow, cexCol=cexCol, hclust=hclustfunc,
                        distfun = distfunc, Rowv=Rowv, Colv=Colv, margins = margins)
        return(hm)
      } else {
        hm <- heatmap.2(gene.counts, scale="none", col=my_palette,
                        labRow = new_labRow,
                        # RowSideColors = row.colors[valid_genes],
                        keysize = 1,
                        key.title = NA,
                        density.info = 'none',
                        key.xlab = key.xlab,
                        dendrogram = 'both',
                        #lhei = c(2, 8),
                        trace="none", cexRow=cexRow, cexCol=cexCol, hclust=hclustfunc,
                        distfun=distfunc, Rowv=Rowv, Colv=Colv,
                        ColSideColors=col.colors, margins = margins)
        return(hm)
      }
    } else {
      if(!is.null(row.colors)){
        hm <- heatmap.2(gene.counts, scale="none", labRow = new_labRow, ##colRow = as.character(row_clusters),
                        hclust=hclustfunc, distfun=distfunc,
                        density.info = 'none', key.title=NA,
                        trace="none",
                        cexRow=cexRow, cexCol = cexCol,
                        Rowv=Rowv, Colv=Colv,
                        #lhei = c(2, 8),
                        col=my_palette, 
                        # ColSideColors = col.colors, 
                        RowSideColors = row.colors[valid_genes],
                        margins = margins)
        return(hm)
      } else {
        hm <- heatmap.2(gene.counts, scale="none", labRow = new_labRow, ##colRow = as.character(row_clusters),
                        hclust=hclustfunc, distfun=distfunc,
                        density.info = 'none', key.title=NA,
                        trace="none",
                        cexRow=cexRow, cexCol = cexCol,
                        Rowv=Rowv, Colv=Colv,
                        #lhei = c(2, 8),
                        col=my_palette, 
                        # ColSideColors = col.colors, 
                        # RowSideColors = row.colors[valid_genes],
                        margins = margins)
        return(hm)
      }
      
    }
    
    return(hm)
  }
}




Get.cell.info <- function(counts){
  cell.info <- data.frame(NAME=colnames(counts), 
                          EXPERIMENT=unlist(strsplit(colnames(counts), "-"))[c(T,F)], 
                          OUTLIER=FALSE, TYPE="UNKNOWN")
  spikein.rows <- grep("ERCC-", row.names(counts))
  if (length(spikein.rows)>0){
    counts <- counts[-c(spikein.rows), ]
  }
  cell.info$NUMREADS <- colSums(counts)
  cell.info$NUMGENES <- apply(counts, 2, function(x){sum(x>0)})
  rownames(cell.info) <- cell.info$NAME
  return(cell.info)
  
}


do.DESeq2 <- function (counts, pheno, keyword = "") {
  
  library(DESeq2)
  sd.data <- apply(counts,1, sd) 
  counts <- subset(counts, sd.data>0) ## keep those genes that have variance > 0
  conditions <- data.frame(colnames(counts)) 
  conditions$pheno <- pheno
  conditions$pheno <- factor(conditions$pheno, levels=unique(pheno))
  
  dds <- DESeqDataSetFromMatrix(countData = counts, colData=conditions, design = ~ pheno +1 )
  dds <- DESeq(dds)
  
  result <- results(dds)
  result <- subset(result, !is.na(result$padj))
  result.sig <- subset(result, abs(log2FoldChange) > 1 & padj < 0.05)
  cat(nrow(result.sig), " genes are significant \n")
  write.csv(result.sig, paste0(keyword, "DESeq2.sig.result.csv"))
  write.csv(result, paste0(keyword, "DESeq2.result.csv"))
  
  pdf(paste0(keyword, "resultMAPlot.pdf"), 6, 6)
  DESeq2::plotMA(result, main="MA Plot", cex.lab = 1.3, cex.axis = 1.3, cex.main = 2)
  dev.off()
  
  
  tophat.rld <- rlog(dds)
  distsRL <- dist(t(assay(tophat.rld)))
  mat <- as.matrix(distsRL)
  rownames(mat) <- colnames(mat) <- colnames(counts)
  library(gplots)
  library(RColorBrewer)
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
#   rownames(mat) <- colnames(assay(tophat.rld))
  
  pdf(paste0(keyword,"heatmap.pdf"))
  heatmap.2(mat, trace="none", col=rev(hmcol), margin=c(13,13))
  dev.off()
  
  return(result)
  
}


testCorelationofGenesWithMDSDims <- function (data, fit.matrix, keyword = "") {
  
  data <- log10(data + 1)
  result <-  t(apply(data, 1, function(x){
    cor.result <-  cor.test(fit.matrix[ , 1], x)
    return(c(cor.result$p.value, cor.result$estimate))
  } ))
  
  result <- data.frame(result)
  colnames(result) <- c("pvalue", "spearman_rho")
  result$TF <- rownames(data) %in% TFs
  write.csv(result, paste0(keyword, "CorWithMDSDim1.csv"))
  
  selectedTFs <- subset(result, result[,1] < 0.05 & abs(result[,2]) > 0.2 & result$TF)
  write.csv(selectedTFs, paste0(keyword, "CorWithMDSDim1sigTF.csv"))
  
  selectedGenes <- subset(result, result[,1] < 0.05 & abs(result[,2]) > 0.2)
  write.csv(selectedGenes, paste0(keyword, "CorWithMDSDim1sigGenes.csv"))
  
  
  result <-  t(apply(data, 1, function(x){
    cor.result <-  cor.test(fit.matrix[ , 2], x)
    return(c(cor.result$p.value, cor.result$estimate))
  } ))
  
  result <- data.frame(result)
  colnames(result) <- c("pvalue", "spearman_rho")
  result$TF <- rownames(data) %in% TFs
  write.csv(result, paste0(keyword, "CorWithMDSDim2.csv"))
  
  selectedTFs <- subset(result, result[,1] < 0.05 & abs(result[,2]) > 0.2 & result$TF)
  write.csv(selectedTFs, paste0(keyword, "CorWithMDSDim2sigTF.csv"))
  
  selectedGenes <- subset(result, result[,1] < 0.05 & abs(result[,2]) > 0.2)
  write.csv(selectedGenes, paste0(keyword, "CorWithMDSDim2sigGenes.csv"))
  
}



checkGeneExpressionforFit <- function (data, fit.matrix) {
  genes <- c("E2F1","CDK1", "CDK2", "CD2", "ZFAND2A", "GFI1", "CDK4", "IRF8", "THBD", "CD1C", "BATF3", "CLEC9A", "KIT", "ITGAX","FLT3", "SPI1", "GAPDH","ERCC-00130", "ID2", "NFIL3", "CCR2", "NOTCH4", "IRF4", "SIGLEC6", "LY6H", "TCF3", "SOX4", "ZEB2", "CEBPD" ) ## "LY6C", "CD24"
  genes <- genes[genes %in% rownames(data)]
  
  for(i in 1: length(genes)){
    
    #   fit.matrix <- fit.MDS.preDC$points
    gene <- genes[i]
    cells <- as.character(colnames(data))
    rownames(fit.matrix) <- cells
    
    left <- rownames(fit.matrix)[fit.matrix[, 1] < 0]
    right <- rownames(fit.matrix)[fit.matrix[, 1] > 0]
    
    col.vector <- as.numeric(log10(data[gene, cells]+1))
    
    pdf(paste0("preDCMDS", "_", genes[i], "onHetGenes.pdf"))
    bp <- gradient_plot_withAVectorXYShapeMatrixInput(col.vector, fit.matrix, color.num = 3,
                                                      title=gene, shapes.progenitor = F, key.title = gene)
    print(bp)
    dev.off()
    
  }
}



getVariableGenesOwn <- function (star.counts, keyword="") {
  library( DESeq )
  library( genefilter )
  library( EBImage )
  library( statmod )
  # library(scLVM)
  # library( topGO )
  # library( org.At.tair.db )
  
#   star.counts <- star.counts[,grepl("PREDC", colnames(star.counts))]  
  
  
  options( max.print=300, width=100 )
  ERCCs <- grepl("ERCC-", rownames(star.counts))
  countsERCC <- star.counts[ERCCs, ]
  
  countsMmus <- star.counts[!ERCCs, ]
  
  
  
  countsMmus <- subset(countsMmus, rowMeans( countsMmus ) > 0 )
  lengthsMmus <- nrow(countsMmus)
  
  sfMmus <- estimateSizeFactorsForMatrix( countsMmus )
  sfERCC <- estimateSizeFactorsForMatrix( countsERCC )
  
  nCountsMmus <- t( t(countsMmus) / sfERCC )
  nCountsERCC <- t( t(countsERCC) / sfERCC )
  #normalise read counts
  #   nCountsERCC <- t( t(countsERCC) / sfERCC )
  
  #normalized counts (brennecke)
  meansMmus <- rowMeans( nCountsMmus )
  varsMmus <- rowVars( nCountsMmus )
  cv2Mmus <- varsMmus / meansMmus^2
  
  
  minMeanForFitA <- unname( quantile( meansMmus[ which( meansMmus > .3 ) ], .95 ) ) ## 0.95 in Bren paper
  useForFitA <- meansMmus >= minMeanForFitA
  fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansMmus[useForFitA] ),
                     cv2Mmus[useForFitA] )
#   ## We use the
#   glmgam.fit function from the statmod package to perform the regression as a GLM fit of the gamma family with log link.
#   The 'cbind' construct serves to produce a model matrix with an intercept.
  xi <- mean( 1 / sfMmus )
  
  a0 <- unname( fit$coefficients["a0"] )
  a1 <- unname( fit$coefficients["a1tilde"] - xi )
 cat("coefficient: ",  c( a0, a1 ))
  
  pdf(paste0(keyword, "CV2meanPlot.pdf"))
  # Prepare the plot (scales, grid, labels, etc.)
  plot( NULL, xaxt="n", yaxt="n",
        log="xy", xlim = c( 1e-1, 3e5 ), ylim = c( .005, 500 ),
        xlab = "average normalized read count", ylab = "squared coefficient of variation (CV^2)" )
  
  axis( 1, 10^(-1:5), c( "0.1", "1", "10", "100", "1000",
                         expression(10^4), expression(10^5) ) )
  axis(2, 10^(-2:1), c("0.01", "0.1", "1", "10"), las=2)
  abline( h=10^(-2:1), v=10^(-1:5), col="#D0D0D0", lwd=2 )
  # Add the data points
  points( meansMmus, cv2Mmus, pch=20, cex=.2, col="blue" )
  # Plot the fitted curve
  xg <- 10^seq( -2, 6, length.out=1000 )
  lines( xg, (xi+a1)/xg + a0, col="#FF000080", lwd=3 )
  # Plot quantile lines around the fit
  df <- ncol(nCountsMmus) - 1
  lines( xg, ( (xi+a1)/xg + a0 ) * qchisq( .975, df ) / df,
         col="#FF000080", lwd=2, lty="dashed" )
  lines( xg, ( (xi+a1)/xg + a0 ) * qchisq( .025, df ) / df,
         col="#FF000080", lwd=2, lty="dashed" )
  
  is_het <- ( (xi+a1)/meansMmus + a0 ) * qchisq( .975, df ) / df < cv2Mmus
  table(is_het) 
  points( meansMmus[is_het], cv2Mmus[is_het], pch=20, cex=.2, col="red" )
  dev.off()
 
 
points( meansMmus[sig], cv2Mmus[sig], pch=20, cex=.2, col="red" )
return(rownames(countsMmus)[is_het])


# meansERCC <- rowMeans( nCountsERCC )
# varsERCC <- rowVars( nCountsERCC )
# cv2ERCC <- varsERCC / meansERCC^2
# minMeanForFitA <- unname( quantile( meansERCC[ which( cv2ERCC > .3 ) ], .8 ) )
# useForFitA <- meansERCC >= minMeanForFitA
# minMeanForFitA
# table( useForFitA )
# 
# fitA <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansERCC[useForFitA] ),
#                     cv2ERCC[useForFitA] )
# residualA <- var( log( fitted.values(fitA) ) - log( cv2ERCC[useForFitA] ) )
# totalA <- var( log( cv2ERCC[useForFitA] ) )
# 1 - residualA / totalA
# 
# pdf("ERCCCVVSMean.pdf")
# plot( meansERCC, cv2ERCC, log="xy", col=1+useForFitA, main="ERCC" )
# xg <- 10^seq( -3, 5, length.out=100 )
# lines( xg, coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/xg )
# segments( meansERCC[useForFitA], cv2ERCC[useForFitA],
#           meansERCC[useForFitA], fitA$fitted.values, col="gray" )
# dev.off()
# 
# minBiolDisp <- .5^2
# xi <- mean( 1 / sfERCC )
# m <- ncol(countsMmus)
# psia1thetaA <- mean( 1 / sfERCC ) +
#   ( coefficients(fitA)["a1tilde"] - xi ) * mean( sfERCC / sfMmus )
# cv2thA <- coefficients(fitA)["a0"] + minBiolDisp + coefficients(fitA)["a0"] * minBiolDisp
# 
# testDenomA <- ( meansMmus * psia1thetaA + meansMmus^2 * cv2thA ) / ( 1 + cv2thA/m )
# pA <- 1 - pchisq( varsMmus * (m-1) / testDenomA, m-1 )
# padjA <- p.adjust( pA, "BH" )
# table( padjA < .1 )
# 
# 
# plot( NULL, xaxt="n", yaxt="n",
#       log="xy", xlim = c( 1e-1, 3e5 ), ylim = c( .005, 100 ),
#       xlab = "average normalized read count", ylab = "squared coefficient of variation (CV^2)" )
# axis( 1, 10^(-1:5), c( "0.1", "1", "10", "100", "1000",
#                        expression(10^4), expression(10^5) ) )
# axis( 2, 10^(-2:2), c( "0.01", "0.1", "1", "10" ,"100"), las=2 )
# abline( h=10^(-2:1), v=10^(-1:5), col="#D0D0D0", lwd=2 )
# # Plot the plant genes, use a different color if they are highly variable
# points( meansMmus, cv2Mmus, pch=20, cex=.2,
#         col = ifelse( padjA < .1, "#C0007090", "#70500040" ) )
# # Add the technical noise fit, as before
# xg <- 10^seq( -2, 6, length.out=1000 )
# lines( xg, coefficients(fitA)["a1tilde"] / xg + a0, col="#FF000080", lwd=3 )
# # Add a curve showing the expectation for the chosen biological CV^2 thershold
# lines( xg, psia1thetaA/xg + coefficients(fitA)["a0"] + minBiolDisp,
#        lty="dashed", col="#C0007090", lwd=3 )
# # Add the normalised ERCC points
# points( meansERCC, cv2ERCC, pch=20, cex=1, col="#0060B8A0" )

}




## create gradient transparant coloring using alpha, n is the number of color spectrum
colorRampAlpha <- function(..., n, alpha) {
  ## example: my_palette <- colorRampAlpha(c("blue", "white", "red"), n = 6,  alpha = 0.8)
  colors <- colorRampPalette(...)(n)
  paste(colors, sprintf("%x", ceiling(255*alpha)), sep="")
}


colCode <- function(values, endCol="blue", b=10) {
  rbPal <- colorRampPalette(c('white', endCol))
  x = rbPal(b)[as.numeric(cut(log(as.numeric(values)+0.1), breaks=b))]
  return(x)
}


transparantGray <- colorRampAlpha("gray", n =1 , alpha = 0.7)
transparantWhite <- colorRampAlpha("white", n =1 , alpha = 0.7)

cell.colors <- c("brown", "orange", "blue", "red", "green", "purple")

alpha <- 1
rgb.colors <- c(rgb(0, 0, 0, alpha), # CDP
                rgb(0, 0, 1, alpha), # CMP
                rgb(0, 1, 1, alpha), # GMDP
                rgb(0, 0.5, 0, alpha), # HSC
                rgb(1, 0, 0, alpha),  # LMPP
                rgb(1, 1, 0, alpha), # MDP
                rgb(1, 0, 1, alpha),  # MLP
                rgb(0, 1, 0, alpha), # MPP
                rgb(1, 0.5, 0, alpha) # BNKP
)

## fit MDS with distance as input
Plot.MDS.with.distance <- function(distances, cells=NULL, colors=NULL, shapes=NULL){
  # cells are the row names
  if (!is.null(cells)){
    fit <- cmdscale(distances[cells, cells], eig=TRUE, k=2)
  } else fit <- cmdscale(distances, eig=TRUE, k=2)
  x <- fit$points[, 1]
  y <- fit$points[, 2]
  if (is.null(colors)){
    colors <- rep("black", length(cells))
    names(colors) <- cells
  }
  if (is.null(shapes)){
    shapes <- rep(16, nrow(fit$points))
    names(shapes) <- rownames(fit$points)
  }
  #par(mar=c(1,1,1,18), xpd=T)
  #par(oma=c(1,1,1,1))
  plot(x,y, cex.lab=1.25, cex.axis=1.25, pch=shapes[rownames(fit$points)],
       col=colors[rownames(fit$points)],
       xlab="MDS Dim 1", ylab="MDS Dim 2", cex=1.5, main="")
  return(fit)
}


## compute iso map distance with k and epsilon
computeDistanceWithKandEpsilon <- function(k = 5, epsilon = 0.999){
  library(vegan)
  diffusion <- computeDiffusionDistMatrix(data.counts, log = T)
  #   diffusion <- as.matrix(dist(mds.fit$points))
  neighborhood.dist <- diffusion
  
  for(i in 1: nrow(diffusion)){
    line <- diffusion[i, ]
    cutoff <- sort(line)[k]
    cat("cutoff of k: ", cutoff, "\n")
    new.cutoff <- min(epsilon, k) ## give away points that are farther than epsilon
    line[line > new.cutoff] <- NA
    neighborhood.dist[i, ] <- line
  }
  iso.distance <- as.matrix(stepacross(neighborhood.dist, toolong = 1))
  return(iso.distance)
}


cleanDatabasebyExpressedGenes <- function (data.file, filter.genes, prefix="") {
  # data.file is the input database see example: /Volumes/big/buffer/TFdatabase/ENCODE_TF_ChIP-seq_2015.txt
  # return ChEA background data
  # return GSEA gmt gene sets file
  
  # data.file <- "/Volumes/big/buffer/TFdatabase/WikiPathways_2015.txt"
  
  encode_data <- readLines(data.file, n = -1, encoding="UTF-8")
    for(i in 1: length(encode_data)){
      line <- as.character(unlist(read.table(text = encode_data[i], sep="\t", stringsAsFactors = F)[1,]))
      members <- line[2:length(line)]
      members.exp <- members[members %in% filter.genes]
      pathwayname <- line[1]
      content <- as.character(c(line[1], "no", members.exp))
      write.table(t(content),paste0(prefix, "_gsea.filtered.gmt"), quote=F, sep="\t", row.names=F, col.names=F, append=TRUE)
    }
    
  ## and then revise the background file for ChEA
  out <- NULL
  for(i in 1: length(encode_data)){
    line <- as.character(unlist(read.table(text = encode_data[i], sep="\t", stringsAsFactors = F)[1,]))
    pathwayname <- line[1]
    members <- line[2:length(line)]
    members.exp <- members[members %in% filter.genes]
    
    if (length(members.exp) == 0)
      next
    
    pubmedID <- pathwayname
    temp <- cbind(pathwayname, pathwayname, members.exp, pubmedID)
    out <- rbind(out, temp)
  }
  
  write.table(out, paste0(prefix, "_chea.filtered.background.csv"), quote=F, sep=",", row.names=T, col.names=F, append=F) 
  #   write.table(out, "test.csv", quote=F, sep=",", row.names=T, col.names=F, append=F) 
  
}
cleanCPDBDatabasebyExpressedGenes <- function (data.file, filter.genes, prefix="") {
  # data.file is the input database see example: /Volumes/big/buffer/TFdatabase/ENCODE_TF_ChIP-seq_2015.txt
  # return ChEA background data, to do
  # return GSEA gmt gene sets file
  
  # data.file <- "/Volumes/big/buffer/TFdatabase/WikiPathways_2015.txt"
  options( warn = -1 )
  encode_data <- readLines(data.file, n = -1)
  for(i in 1: length(encode_data)){
    line <- as.character(unlist(read.table(text = encode_data[i], sep="\t", stringsAsFactors = F)[1,]))
    members <- as.character(unlist(strsplit(line[4],  ",")))
    members.exp <- members[members %in% filter.genes]
    if(length(members.exp) == 0)
      next
    pathwayname <- paste0(i, "_", line[1])
    sourcePath <- line[3]
    content <- as.character(c(pathwayname, sourcePath, members.exp))
    write.table(t(content),paste0(prefix, "_gsea.filtered.gmt"), quote=F, sep="\t", row.names=F, col.names=F, append=TRUE)
  }
  
#   ## and then revise the background file for ChEA
#   out <- NULL
#   for(i in 1: length(encode_data)){
#     line <- as.character(unlist(read.table(text = encode_data[i], sep="\t", stringsAsFactors = F)[1,]))
#     pathwayname <- line[1]
#     members <- line[2:length(line)]
#     members.exp <- members[members %in% filter.genes]
#     
#     if (length(members.exp) == 0)
#       next
#     
#     pubmedID <- pathwayname
#     temp <- cbind(pathwayname, pathwayname, members.exp, pubmedID)
#     out <- rbind(out, temp)
#   }
#   
#   write.table(out, paste0(prefix, "_chea.filtered.background.csv"), quote=F, sep=",", row.names=T, col.names=F, append=F) 
#   #   write.table(out, "test.csv", quote=F, sep=",", row.names=T, col.names=F, append=F) 
  
}


drawHeatmapWithDistance <- function (mat) {
  library(gplots)
  library(RColorBrewer)
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  heatmap.2(mat, trace="none", col=rev(hmcol), margin=c(13,13))
}



z.score <- function(x){
  
  mean <- mean(x)
  sd <- sd(x)
  zscore <- (x-mean)/sd
  
  return(zscore)
}



gradient_plot_withAVectorXYShape <- function (col.vector, mds.result, title="", gradientColor = "blue", shapes.progenitor=T, key.title = "Yield(log10)") {
  
  library(ggplot2)
  x <- mds.result$points[,1]
  y <- mds.result$points[,2]
  dat <- data.frame(x=x, y=y)
  bp <- qplot(x, y, data=dat, colour=col.vector, size = I(0.1), asp=1, main=title, xlab="MDS PC1", ylab="MDS PC2")  ## make it super small first to revise the shapes later 
  bp
  bp <- bp + scale_colour_gradient(name=key.title, low="blue", high=gradientColor) 

  if(shapes.progenitor){
  bp <- bp + geom_point(data = dat, size = I(3), aes(shape = factor(data.counts[,1]))) + labs(shape="Progenitor type") + scale_shape_manual(values=progenitor.shapes)
  } else {
  bp <- bp + geom_point(data = dat, shape = 16, solid=F, size = I(3))
  }
  bp
}

gradient_plot_withAVectorXYShapeMatrixInput <- function (col.vector, matrix, color.num = 2,  title="", dot_size = 4,  axis_textSize = 15, 
                                                         gradientColor = "blue", baseColor = "white", midColor = "white",
                                                         shapes.progenitor=T, shapes.vec = 16,  key.title = "Yield(log10)",
                                                         low.cutoff = 10^-10, xlab="", ylab = "", legend=T, legend.key.size = unit(1, "cm")) {
  
  library(ggplot2)
  x <- matrix[,1]
  y <- matrix[,2]
  dat <- data.frame(x=x, y=y)
  bp <- qplot(x, y, data=dat, colour=col.vector, asp=1, main="", xlab=xlab, ylab=ylab)  ## make it super small first to revise the shapes later 
  bp <- bp + geom_point(data = dat, shape = shapes.vec, size = I(dot_size))
  if(color.num == 2){
    ## this one is useful for the NI paper to plot the lineage specific pca or tnse.
    bp <- bp + scale_colour_gradient(name = key.title, low = baseColor, high = gradientColor, limits = c(low.cutoff, max(col.vector)), na.value = "gray")
  } else {
    my_palette <- colorRampAlpha(c( "gray","blue", "red"), n = 1000,  alpha = 0.9)
    # my_palette <- colorRampAlpha(c( "#EF00F2","black", "#F9FE0D"), n = 1000,  alpha = 0.7)
    bp <- bp + scale_colour_gradientn(name=key.title, colours = my_palette)
  }
  
  bp <- bp +  theme( plot.margin=unit(c(0.05,0.05,0.05,0.05), "cm"), axis.text.x = element_text(size = axis_textSize), axis.text.y = element_text(size = axis_textSize),             
                     axis.title.x = element_text(size = 22), axis.title.y = element_text(size = 22, angle = 90))
  
  if(legend){
    bp <- bp +  theme(legend.text = element_text(size = axis_textSize), legend.title = element_text(size = axis_textSize, colour = "black"))
    bp <- bp + theme(legend.background = element_rect(), legend.margin = unit(0.1, "cm"))
    bp <- bp + theme(legend.key.size = legend.key.size)
  } else {
    bp <- bp + theme(legend.position = "none") ## remove the legend
  }
  #   bp <- bp + theme_bw() ## make it white
  #   bp <- bp + theme( panel.grid.major = element_blank(),
  #                     panel.grid.minor = element_blank()
  #                     )
  #   bp <- bp + theme(axis.line = element_line(colour = "black"))
  #   bp
  return(bp)
}


computeDiffusionDimsMDS <- function (data, quantile.cutoff = 0.4, k = 2, legend.position = "bottomright", pdfPrefix="", log = TRUE) {
  
  n <- ncol(data)
  
  totalCellOutput <- rowSums(data[, 2:n])
  totalCellOutput[ totalCellOutput < 200] <- 0.5
  totalCellOutput[ totalCellOutput >= 200 & totalCellOutput < 800] <- 0.9
  totalCellOutput[ totalCellOutput >= 800 & totalCellOutput < 10000] <- 1.2
  totalCellOutput[ totalCellOutput >= 10000 &  totalCellOutput < 100000] <- 1.4
  totalCellOutput[ totalCellOutput >= 100000] <- 2
  
  alpha <- 0.65
  rgb.colors <- c(rgb(0, 0, 0, alpha), # CDP
                  rgb(0, 0, 1, alpha), # CMP
                  rgb(0, 1, 1, alpha), # GMDP
                  rgb(0, 0.5, 0, alpha), # HSC
                  rgb(1, 0, 0, alpha),  # LMPP
                  rgb(1, 1, 0, alpha), # MDP
                  rgb(1, 0, 1, alpha),  # MLP
                  rgb(0, 1, 0, alpha), # MPP
                  rgb(1, 0.5, 0, alpha) # BNKP
  )
  
  
  sampleDists <-  computeDiffusionDistMatrix(data = data, quantile.cutoff = quantile.cutoff, log = log)
  
  
  colors <- generateTransparentColors(data, alpha = 0.65)
  circle.colors <- generateTransparentColors(data, alpha = 1)
  
  if(k == 3){
    # MDS calling function
    fit <- cmdscale(sampleDists, eig=TRUE, k=k)
    mds.points <- fit$points[, 1:3]
    library(rgl)
    plot3d(mds.points, xlab = "MDS Dim1", ylab = "MDS Dim2", zlab = "MDS Dim3", col = colors, type='s', size=totalCellOutput)
    legend3d("topright", legend = c("CDP", "CMP", "GMDP", "HSC", "LMPP", "MDP", "MLP", "MPP", "BNKP"), pch = 16, col = c(seq(1:8), "saddlebrown"), cex=1.2, inset=c(0.02))
    x <- fit$points[, 1]
    y <- fit$points[, 2]
    z <- fit$points[, 3]
    
    pdf(paste0(pdfPrefix,"_PC12diffusionMap.pdf"))
    plot(x,y, cex.lab=1.25, cex.axis=1.25,
         col=circle.colors[rownames(fit$points)], pch=1,
         xlab="MDS Dim 1", ylab="MDS Dim 2", cex = totalCellOutput)
    
    points(x,y, cex.lab=1.25, cex.axis=1.25,
           col=colors[rownames(fit$points)], pch=16,
           xlab="MDS Dim 1", ylab="MDS Dim 2", cex = totalCellOutput)
    
    legend(legend.position, legend = c("CDP", "CMP", "GMDP", "HSC", "LMPP", "MDP", "MLP", "MPP", "BNKP"), pch = 16, col = rgb.colors, cex=0.7)
    dev.off()
    
    pdf(paste0(pdfPrefix,"_PC13diffusionMap.pdf"))
    plot(x,z, cex.lab=1.25, cex.axis=1.25,
         col=circle.colors[rownames(fit$points)], pch=1,
         xlab="MDS Dim 1", ylab="MDS Dim 3", cex = totalCellOutput)
    
    points(x,z, cex.lab=1.25, cex.axis=1.25,
           col=colors[rownames(fit$points)], pch=16, cex = totalCellOutput)
    
    legend(legend.position, legend = c("CDP", "CMP", "GMDP", "HSC", "LMPP", "MDP", "MLP", "MPP", "BNKP"), pch = 16, col = rgb.colors, cex=0.7)
    dev.off()
    
    pdf(paste0(pdfPrefix,"_PC23diffusionMap.pdf"))
    plot(y,z, cex.lab=1.25, cex.axis=1.25,
         col=circle.colors[rownames(fit$points)], pch=1,
         xlab="MDS Dim 2", ylab="MDS Dim 3", cex = totalCellOutput)
    
    points(y,z, cex.lab=1.25, cex.axis=1.25,
           col=colors[rownames(fit$points)], pch=16,
           cex = totalCellOutput)
    
    legend(legend.position, legend = c("CDP", "CMP", "GMDP", "HSC", "LMPP", "MDP", "MLP", "MPP", "BNKP"), pch = 16, col = rgb.colors, cex=0.7)
    dev.off()
  }
  
  
  if(k == 2) {
    fit <- cmdscale(sampleDists, eig=TRUE, k=k)
    x <- fit$points[, 1]
    y <- fit$points[, 2]
    
    pdf(paste0(pdfPrefix,"_2DdiffusionMap.pdf"))
    plot(x,y, cex.lab=1.25, cex.axis=1.25,
         col=circle.colors[rownames(fit$points)], pch=1,
         xlab="MDS Dim 1", ylab="MDS Dim 2", cex = totalCellOutput)
    
    points(x,y, cex.lab=1.25, cex.axis=1.25,
           col=colors[rownames(fit$points)], pch=16,
           xlab="MDS Dim 1", ylab="MDS Dim 2", cex = totalCellOutput)
    
    legend(legend.position, legend = c("CDP", "CMP", "GMDP", "HSC", "LMPP", "MDP", "MLP", "MPP", "BNKP"), pch = 16, col = rgb.colors, cex=0.6)
    dev.off()
  }
  
  
  return(fit)
}



plotPersonalizedDiffusionMap <- function (mds.fit, shapes = 16, colors = "blue", bordered.colors = "blue", dot.size = 1) {
  ## after obtaining a mds fit, plot it with different colors, size and shapes
#   plotPersonalizedDiffusionMap(mds.fit)
  x <- mds.fit$points[, 1]
  y <- mds.fit$points[, 2]
  
  plot(x,y, cex.lab=1.25, cex.axis=1.25, col=colors[rownames(mds.fit$points)], pch=shapes, xlab="MDS Dim 1", ylab="MDS Dim 2", cex = dot.size)
  points(x, y, cex.lab=1.25, cex.axis=1.25,
         col=bordered.colors[rownames(mds.fit$points)], pch=shapes-15,
         xlab="MDS Dim 1", ylab="MDS Dim 2", cex = dot.size)
}

plotPersonalizedt_SNEMap <- function (tsne.map, shapes = 16, colors = "blue", bordered.colors = "blue", dot.size = 1) {
  ## after obtaining a mds fit, plot it with different colors, size and shapes
  #   plotPersonalizedDiffusionMap(mds.fit)
  x <- tsne.map[, 1]
  y <- tsne.map[, 2]
  
  plot(x,y, cex.lab=1.25, cex.axis=1.25, col=colors, pch=shapes, xlab="t-SNE Dim 1", ylab="t-SNE Dim 2", cex = dot.size)
  points(x, y, cex.lab=1.25, cex.axis=1.25,
         col=bordered.colors, pch=shapes-15,
          cex = dot.size)
}



normalize.with.batch <- function(data, batch, method){
  ## batch is a vector that specify the batches for each sample.
#   sd.data <- apply(data,1, sd) 
#   mean.data <- apply(data,1,mean)
#   data <- subset(data, sd.data>0 & mean.data>2 ) ## keep those genes that have variance > 0 and mean count >2 
#   
  if(method=="limma"){
    library(limma)
    batch.norm.data <- removeBatchEffect(data, batch=batch)
  }
  else {
    library(sva)
    batch.norm.data <- ComBat(data, batch=batch, mod=NULL, numCovs = NULL, par.prior = TRUE,prior.plots = FALSE)
  }
  
  return(batch.norm.data)
}

generateTransparentColors <- function (data, alpha = 0.65) {
  
  colors <- rep("black", nrow(data))
  names(colors) <- rownames(data)
  
  rgb.colors <- c(rgb(0, 0, 0, alpha), # CDP
                  rgb(0, 0, 1, alpha), # CMP
                  rgb(0, 1, 1, alpha), # GMDP
                  rgb(0, 0.5, 0, alpha), # HSC
                  rgb(1, 0, 0, alpha),  # LMPP
                  rgb(1, 1, 0, alpha), # MDP
                  rgb(1, 0, 1, alpha),  # MLP
                  rgb(0, 1, 0, alpha), # MPP
                  rgb(1, 0.5, 0, alpha) # BNKP
  )
  
  colors[grepl("-CDP", rownames(data))] <- rgb.colors[1]
  colors[grepl("-CMP", rownames(data))] <- rgb.colors[2]
  colors[grepl("-GMDP", rownames(data))] <-rgb.colors[3]
  colors[grepl("-HSC", rownames(data))] <- rgb.colors[4] 
  colors[grepl("-LMPP", rownames(data))] <- rgb.colors[5]
  colors[grepl("-MDP", rownames(data))] <- rgb.colors[6]
  colors[grepl("-MLP", rownames(data))] <- rgb.colors[7]
  colors[grepl("-MPP", rownames(data))] <- rgb.colors[8]
  colors[grepl("-BNKP", rownames(data))] <- rgb.colors[9]
  
  colors
}


generateProgenitorShapes <- function (data) {
  
  shapes <- rep(16, nrow(data))
  names(shapes) <- rownames(data)
  
  progenitor.shapes <- c(15, # CDP
                         16, # CMP
                         17, # GMDP
                         18, # HSC
                         7,  # LMPP
                         8, # MDP
                         9,  # MLP
                         10, # MPP
                         11 # BNKP
  )
  
  shapes[which(data[, 1]=="CDP")] <- progenitor.shapes[1]
  shapes[which(data[, 1]=="CMP")] <- progenitor.shapes[2]
  shapes[which(data[, 1]=="GMDP")] <-progenitor.shapes[3] 
  shapes[which(data[, 1]=="HSC")] <- progenitor.shapes[4] 
  shapes[which(data[, 1]=="LMPP")] <- progenitor.shapes[5]
  shapes[which(data[, 1]=="MDP")] <- progenitor.shapes[6] ## this can include GMDP
  shapes[which(data[, 1]=="MLP")] <- progenitor.shapes[7]
  shapes[which(data[, 1]=="MPP")] <- progenitor.shapes[8]
  shapes[which(data[, 1]=="BNKP")] <- progenitor.shapes[9]
  
  shapes
}



generateTransparentColorsList <- function (cell.matrix, alpha = 0.65) {
  
  colors <- rep("black", nrow(cell.matrix))
  names(colors) <- rownames(cell.matrix)
  
  rgb.colors <- c(# rgb(0, 0, 0, alpha), # CDP
                  transparentColor("brown", alpha = alpha), # CDP
                  rgb(0, 0, 1, alpha), # CMP
                  rgb(0, 1, 1, alpha), # GMDP
                  rgb(0, 0, 0, alpha), # HSC
                  rgb(1, 0, 0, alpha),  # LMPP
                  rgb(1, 1, 0, alpha), # MDP
                  rgb(1, 0, 1, alpha),  # MLP
                  # rgb(0, 1, 0, alpha), # MPP
                  transparentColor("gray", alpha = alpha), # MPP
                  rgb(1, 0.5, 0, alpha) # BNKP
  )
  
  colors[which(cell.matrix[, 1]=="CDP")] <- rgb.colors[1]
  colors[which(cell.matrix[, 1]=="CMP")] <- rgb.colors[2]
  colors[which(cell.matrix[, 1]=="GMDP")] <-rgb.colors[3]
  colors[which(cell.matrix[, 1]=="HSC")] <- rgb.colors[4] 
  colors[which(cell.matrix[, 1]=="LMPP")] <- rgb.colors[5]
  colors[which(cell.matrix[, 1]=="MDP")] <- rgb.colors[6]
  colors[which(cell.matrix[, 1]=="MLP")] <- rgb.colors[7]
  colors[which(cell.matrix[, 1]=="MPP")] <- rgb.colors[8]
  colors[which(cell.matrix[, 1]=="BNKP")] <- rgb.colors[9]
  
  return(list(rgb.colors, colors))
}

sc.tsne <- function (heatmap.data, k=2, perplexity=30, col.colors="blue") {
  
  library(tsne)
  sampleDists <- as.dist(1-cor(log2(heatmap.data+1),method="spearman"))
  
#   cor.spearman <- cor(heatmap.data,method="spearman")
#   for(i in 1: ncol(heatmap.data))
#     cor.spearman[i,i] <- 0
  
#   sampleDists <- as.dist(cor.spearman)

  x <- tsne(sampleDists,k = k, perplexity = perplexity)
  plot(x, cex.lab=1.25, cex.axis=1.25, col=col.colors,pch=16,
       xlab="tsne Dim 1", ylab="tsne Dim 2", cex=1.5)  
  return(x)
}


CVVSMeanPlotERCC <- function (star.counts) {
  spikein.rows <- grep("ERCC-", row.names(star.counts))
  mapping.rows <- grep("^__", row.names(star.counts))
  genes.counts <- star.counts[-c(spikein.rows, mapping.rows), ]
  spikein.data <- star.counts[spikein.rows,]
  
  #2. calculate normalisation for counts
  countsMmus <- genes.counts
  countsERCC <- spikein.data
  lengthsMmus <- nrow(countsMmus)
  lengthsERCC <- nrow(countsERCC)
  
  
  sfERCC <- estimateSizeFactorsForMatrix( countsERCC )
  sfMmus <- sfERCC #also use ERCC size factor for endogenous genes
  
  
  #normalise read counts
  nCountsERCC <- t( t(countsERCC) / sfERCC )
  nCountsMmus <- t( t(countsMmus) / sfMmus )
  
  
  #normalized counts (brennecke)
  meansMmus <- rowMeans( nCountsMmus )
  varsMmus <- rowVars( nCountsMmus )
  cv2Mmus <- varsMmus / meansMmus^2
  
  meansERCC <- rowMeans( nCountsERCC )
  varsERCC <- rowVars( nCountsERCC )
  cv2ERCC <- varsERCC / meansERCC^2
  
  #Do fitting of technical noise
  
  #normalised counts (with size factor)
  minMeanForFitA <- unname( quantile( meansERCC[ which( cv2ERCC > .3 ) ], .1 ) )
  useForFitA <- meansERCC >= minMeanForFitA
  fitA <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansERCC[useForFitA] ),
                      cv2ERCC[useForFitA] )
  
  pdf("ERCCFitPlot.pdf")
  #plot fit
  plot( meansERCC, cv2ERCC, log="xy", col=1+useForFitA)
  
  
  xg <- 10^seq( -3, 5, length.out=100 )
  lines( xg, coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/xg )
  segments( meansERCC[useForFitA], cv2ERCC[useForFitA],
            meansERCC[useForFitA], fitA$fitted.values, col="gray" )
  dev.off()
}


get.pca.loadings <- function (pca, topNum = 100, prefix = "") {
  # data <- log2(data+1)
  # sd.data <- apply(data,1, sd) 
  # data <- subset(data, sd.data>0) ## keep those genes that have variance > 0
  # pca <- prcomp(t(data),retx=TRUE,center=TRUE,scale=TRUE,cor=TRUE)
  # 
  # write.csv(pca$rotation[,1:2], "pca.rotations.csv")
  PC1 <- abs(pca$rotation[,1])
  PC2 <- abs(pca$rotation[,2])
  PC3 <- abs(pca$rotation[,3])
  
  t.PC1 <- sort(PC1, decreasing = T)[topNum]
  t.PC2 <- sort(PC2, decreasing = T)[topNum]
  t.PC3 <- sort(PC3, decreasing = T)[topNum]
  
  PC1.dim <- PC1[PC1>=t.PC1]
  PC2.dim <- PC2[PC2>=t.PC2]
  PC3.dim <- PC3[PC3>=t.PC2]
  
  write.csv(pca$rotation[,1][names(PC1) %in% names(PC1.dim)], paste0(prefix,"_PC1.loadings.csv"))
  write.csv(pca$rotation[,2][names(PC2) %in% names(PC2.dim)], paste0(prefix,"_PC2.loadings.csv"))
  write.csv(pca$rotation[,3][names(PC3) %in% names(PC3.dim)], paste0(prefix,"_PC3.loadings.csv"))
  
  # write.table(names(PC1)[PC1>t.PC1],"PC1.genes.txt", row.names = F, col.names = F, quote=F )
  # write.table(names(PC2)[PC2>t.PC2],"PC2.genes.txt", row.names = F, col.names = F, quote=F )
  # write.table(names(PC3)[PC3>t.PC3],"PC3.genes.txt", row.names = F, col.names = F, quote=F )
  # 
  # write.csv(pca$x[,1], "PC1.csv")
  # write.csv(pca$x[,2], "PC2.csv")
  # write.csv(pca$x[,3], "PC3.csv")
  
}


sc.getCodingGenes <- function(){
  coding.genes <- read.table("/Volumes/big/buffer/coding.genes.txt")[,1]
  return(coding.genes)
}

sc.topFive <- function(counts){
  qc.rows <- grep("^__", row.names(counts))
  spikein.rows <- grep("ERCC-", row.names(counts))
  if (length(c(qc.rows, spikein.rows))>0){
    counts <- counts[-c(qc.rows, spikein.rows), ]
  }
  fraction.counts <- sc.NormalizeData(counts, method="sum")
  topfivenames <- apply(fraction.counts, 2, function(x){head(names(sort(x, decreasing=T)), 5)})
  topfivefractions <- apply(fraction.counts, 2, function(x){head(sort(x, decreasing=T), 5)})
  topfivelist <- vector('list', ncol(counts))
  for (i in 1:ncol(counts)){
    sample.name <- names(counts)[i]
    topfivelist[[i]] <- data.frame(SAMPLE=rep(sample.name, 5), 
                                   FRACTION=topfivefractions[, sample.name], 
                                   GENE=topfivenames[, sample.name], 
                                   RANK=1:5)
  }
  topfive <- do.call('rbind', topfivelist)
}


sc.plotTopFive <- function(topfive){
  par(mar=c(10,4,2,2))
  boxplot(FRACTION ~ GENE, subset(topfive, RANK==1), varwidth=T, las=2, ylim=c(0,1), 
          ylab="Fraction of uniquely mapped reads that map to top gene")
}



sc.QC <- function(counts, hk.genes=NULL, other.genes=NULL){
  library(plyr)
  qc.rows <- grep("^__", row.names(counts))
  if (length(qc.rows) > 0){
    counts <- counts[-qc.rows, ]
  }
  spikein.rows <- grep("ERCC-", row.names(counts))
  htseq.summary <- as.data.frame(t(apply(counts, 2, function(x){
    c(R_UNIQ_MAP_star=sum(x), 
      R_UNIQ_MAP_ERCC_star=sum(x[spikein.rows]),
      R_UNIQ_MAP_GENES_star=sum(x[-spikein.rows]),
      G0_GENES_star=sum(x[-spikein.rows]>0), 
      G1_GENES_star=sum(x[-spikein.rows]>1),
      G5_GENES_star=sum(x[-spikein.rows]>5),
      G0_ERCC_star=sum(x[spikein.rows]>0), 
      G1_ERCC_star=sum(x[spikein.rows]>1), 
      G5_ERCC_star=sum(x[spikein.rows]>5))
  })))
  htseq.summary$P_UNIQ_MAP_ERCC_star <- htseq.summary$R_UNIQ_MAP_ERCC_star / htseq.summary$R_UNIQ_MAP_star
  if (!is.null(hk.genes)){
    hk.rows <- which(hk.genes %in% rownames(counts))
    hk.info <- as.data.frame(t(apply(counts, 2, function(x){
      c(G0_HK_GENES_star=sum(x[hk.rows]>0),
      G5_HK_GENES_star=sum(x[hk.rows]>5))
    })))
    htseq.summary <- cbind(htseq.summary, hk.info)
  }
  if (!is.null(other.genes)){
    other.rows <- which(other.genes %in% rownames(counts))
    other.info <- as.data.frame(t(apply(counts, 2, function(x){
      c(G0_OTHER_GENES_star=sum(x[other.rows]>0),
        G5_OTHER_GENES_star=sum(x[other.rows]>5))
    })))
    htseq.summary <- cbind(htseq.summary, other.info)
  }
  htseq.summary$SAMPLE_NAME <- colnames(counts)
  htseq.summary$Expressed_Genes_Per <- htseq.summary$G1_GENES_star/as.numeric(nrow(counts))
  return(htseq.summary)
}

sc.plotRawReads <- function(htseq.summary){
  hist(log10(htseq.summary$R_TOT), main="",
       cex.lab=1.5, cex.axis=1.5, cex.main=1, pch=16, breaks=20,
       cex.sub=1.5, xlab="log10(# raw reads per cell)", ylab="Count")
}

sc.plotExpressedGenes <- function(htseq.summary){
  hist(htseq.summary$G1_GENES_star, main="",
       cex.lab=1.5, cex.axis=1.5, cex.main=1, pch=16, breaks=20,
       cex.sub=1.5, xlab="# genes that have more than 1 reads mapped", ylab="Count")
}


sc.plotGeneReads <- function(htseq.summary){
  hist(log10(htseq.summary$R_UNIQ_MAP_GENES_star), main="",
       cex.lab=1.5, cex.axis=1.5, cex.main=1, pch=16, breaks=20,
       cex.sub=1.5, xlab="log10(# reads mapped to genes per cell)", ylab="Count")
}
sc.plotMappingPercentage <- function(htseq.summary){
  hist(htseq.summary$M_P_star, main="",
       cex.lab=1.5, cex.axis=1.5, cex.main=1, pch=16, breaks=20, 
       cex.sub=1.5, xlim=c(0, 1), xlab="Mapping % per cell", ylab="Count")
}
sc.plotSpikeinPercentage <- function(htseq.summary){
  hist(htseq.summary$P_UNIQ_MAP_ERCC_star, main="",
       cex.lab=1.5, cex.axis=1.5, cex.main=1, pch=16, breaks=20,
       cex.sub=1.5, xlim=c(0, 1), xlab="Mapping % to spike-in per cell", ylab="Count")
}
neReadsvsSpikeinPercentage <- function(htseq.summary){
  plot(htseq.summary$P_UNIQ_MAP_ERCC_star, log10(htseq.summary$R_UNIQ_MAP_GENES_star),
       cex.lab=1.5, cex.axis=1.5, cex.main=1, pch=16, 
       cex.sub=1.5, xlim=c(0, 1), xlab="% uniquely mapped reads mapping to spike-in per cell",
       ylab="log10 (# reads mapping to genes per cell)")
}
neReadsvsSpikeinPercentage <- function(htseq.summary){
  plot(htseq.summary$P_UNIQ_MAP_ERCC_star, log10(htseq.summary$R_UNIQ_MAP_GENES_star),
       cex.lab=1.5, cex.axis=1.5, cex.main=1, pch=16, 
       cex.sub=1.5, xlim=c(0, 1), xlab="% uniquely mapped reads mapping to spike-in per cell",
       ylab="log10 (# reads mapping to genes per cell)")
}

sc.plotNumberGenesMoreThan5vsGeneReads <- function(htseq.summary, colors=NULL){
  if (is.null(colors)){
    colors <- rep("black", nrow(htseq.summary))
    names(colors) <- rownames(htseq.summary)
  }
  plot(log10(htseq.summary$R_UNIQ_MAP_GENES_star), htseq.summary$G5_GENES_star,
       cex.lab=1.5, cex.axis=1.5, cex.main=1, col=colors[rownames(htseq.summary)],
       cex.sub=1.5, xlab="log10(# reads uniquely mapping to genes per cell)",
       ylab="# genes detected per cell (>5 counts)")
}

sc.plotNumberGenesMoreThan1vsGeneReads <- function(htseq.summary, colors=NULL){
  if (is.null(colors)){
    colors <- rep("black", nrow(htseq.summary))
    names(colors) <- rownames(htseq.summary)
  }
  plot(log10(htseq.summary$R_UNIQ_MAP_GENES_star), pch=16,htseq.summary$G1_GENES_star,
       cex.lab=1.5, cex.axis=1.5, cex.main=1, col=colors[rownames(htseq.summary)],
       cex.sub=1.5, xlab="log10(# reads mapped to genes per cell)",
       ylab="# genes detected per cell (>1 counts)")
}

sc.SpikeInRegression <- function(counts, colors=NULL, do.plot=T){
  library(plyr)
  spikein.rows <- grep("ERCC-", row.names(counts))
  spikein.counts <- counts[spikein.rows, ]
  ercc.controls <- read.delim("ERCC_Controls_Analysis.txt", 
                              header=T, check.names=F)
  rownames(ercc.controls) <- ercc.controls[, "ERCC ID"]
  concentration <- ercc.controls[row.names(spikein.counts), 
                                 "concentration in Mix 2 (attomoles/ul)"]
  if (do.plot) {
    plot(c(),c(),xlim=range(log10(concentration)), ylim=range(log10(spikein.counts+1)), 
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5,  
       xlab="log10(concentration in spike-in mix)", 
       ylab="log10(# reads of spike-in transcript)")
  }
  if (is.null(colors)){
    colors <- rep("black", ncol(counts))
    names(colors) <- colnames(counts)
  }
  slopes <- rep(0, ncol(counts))
  names(slopes) <- colnames(counts)
  for (cell.name in colnames(counts)){
    fit <- lm(log10(spikein.counts[, cell.name]+1) ~ log10(concentration))
    if (do.plot){
      abline(fit, col=colors[cell.name])
      points(log10(concentration), log10(spikein.counts[, cell.name]+1), col=colors[cell.name])
    }
    slopes[cell.name] <- fit$coefficients[1]
  }
  return(slopes)
}

sc.GetDropoutFitConc <- function(counts){
  library(plyr)
  spikein.rows <- grep("ERCC-", row.names(counts))
  spikein.counts <- counts[spikein.rows, ]
  ercc.controls <- read.delim("ERCC_Controls_Analysis.txt", 
                              header=T, check.names=F)
  rownames(ercc.controls) <- ercc.controls[, "ERCC ID"]
  concentration <- ercc.controls[row.names(spikein.counts), 
                                 "concentration in Mix 2 (attomoles/ul)"]
  dropout.threshold <- 0
  dropout <- logical(0)
  conc <- numeric(0)
  for (i in 1:ncol(counts)){
    s.spikein.counts <- as.numeric(spikein.counts[, i])  
    this.dropout <- s.spikein.counts  <= dropout.threshold
    dropout <- c(dropout, this.dropout)
    conc <- c(conc, concentration)
  }
  dropout.df <- as.data.frame(cbind(conc, dropout))
  glm.fit <- glm(dropout ~ conc, data=dropout.df, family="binomial")
  return(glm.fit)
}

sc.GetDropoutFit <- function(counts){
  library(plyr)
  spikein.rows <- grep("ERCC-", row.names(counts))
  spikein.counts <- counts[spikein.rows, ]
  ercc.controls <- read.delim("ERCC_Controls_Analysis.txt", 
                              header=T, check.names=F)
  rownames(ercc.controls) <- ercc.controls[, "ERCC ID"]
  concentration <- ercc.controls[row.names(spikein.counts), 
                                 "concentration in Mix 2 (attomoles/ul)"]
  dropout.threshold <- 0
  dropout <- logical(0)
  pred.count <- numeric(0)
  for (i in 1:ncol(counts)){
    s.spikein.counts <- as.numeric(spikein.counts[, i])
    lm.fit <- lm(s.spikein.counts ~ concentration - 1)    
    this.dropout <- s.spikein.counts  <= dropout.threshold
    this.pred.count <- round(predict(lm.fit, data.frame(concentration=concentration)), 2)
    dropout <- c(dropout, this.dropout)
    pred.count <- c(pred.count, this.pred.count)
    print(names(counts)[i])
    print(max(this.pred.count[this.dropout]))
    browser()
  }
  dropout.df <- as.data.frame(cbind(pred.count, dropout))
  glm.fit <- glm(dropout ~ pred.count, data=dropout.df, family="binomial")
  return(glm.fit)
}

sc.GetDropoutPrediction <- function(counts, predict){
  library(plyr)
  spikein.rows <- grep("ERCC-", row.names(counts))
  spikein.counts <- counts[spikein.rows, , drop=F]
  ercc.controls <- read.delim("ERCC_Controls_Analysis.txt", 
                              header=T, check.names=F)
  rownames(ercc.controls) <- ercc.controls[, "ERCC ID"]
  concentration <- ercc.controls[row.names(spikein.counts), 
                                 "concentration in Mix 2 (attomoles/ul)"]
  dropout.threshold <- 2
  lm.fit <- lm(spikein.counts[, 1] ~ concentration - 1)    
  plot(concentration, spikein.counts[, 1], main=names(counts)[1])
  abline(lm.fit)
  pred.count <- round(predict(lm.fit, data.frame(concentration=concentration)), 2)
  s.spikein.counts <- spikein.counts[, 1]
  dropout <- as.numeric(s.spikein.counts <= dropout.threshold)
  dropout.df <- arrange(as.data.frame(cbind(pred.count, dropout)), pred.count)
  glm.fit <- glm(dropout ~ pred.count, data=dropout.df, family="binomial")
  predict.df <- data.frame(pred.count=predict)
  predictions <- predict(glm.fit, predict.df, type="response")
  return(predictions)
}

sc.plotDropoutProbability <- function(counts, colors=NULL, do.plot=T, normalize=F){
  library(plyr)
#   counts <- spikein.data
#   normalize <- F
#   do.plot=T
#   colors=NULL
  
  sample.names <- colnames(counts)
  spikein.rows <- grep("ERCC-", row.names(counts))
  spikein.counts <- counts[spikein.rows, ]
  if (normalize){
    spikein.counts <- round(apply(spikein.counts, 2, function(x){x/sum(x)*1000000}))
  }
  ercc.controls <- read.delim("ERCC_Controls_Analysis.txt", 
                              header=T, check.names=F)
  rownames(ercc.controls) <- ercc.controls[, "ERCC ID"]
  concentration <- ercc.controls[row.names(spikein.counts), 
                                 "concentration in Mix 2 (attomoles/ul)"]
  dropout.threshold <- 0
  predictions <- matrix(nrow=10000, ncol=length(sample.names))
  counts.seq <- seq(1, 10000)
  colnames(predictions) <- sample.names
  if (is.null(colors)){
    colors <- rep("black", ncol(counts))
    names(colors) <- colnames(counts)
  }
  for (i in 1:length(sample.names)){
    s.name <- sample.names[i]
    s.spikein.counts <- spikein.counts[, s.name]
    lm.fit <- lm(s.spikein.counts ~ concentration - 1)    
    pred.count <- round(predict(lm.fit, data.frame(concentration=concentration)), 2)
    dropout <- as.numeric(s.spikein.counts <= dropout.threshold)
    dropout.df <- arrange(as.data.frame(cbind(pred.count, dropout)), pred.count)
    glm.fit <- glm(dropout ~ pred.count, data=dropout.df, family="binomial")
    predict.df <- data.frame(pred.count=counts.seq)
    predictions[, i] <- predict(glm.fit, predict.df, type="response")
  }
  plot(c(),c(),xlim=range(counts.seq), ylim=c(0,1), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, 
       xlab="counts", ylab="drop-out probability")
  for (cell.name in colnames(counts)){
    if (do.plot){
      lines(counts.seq, predictions[, cell.name], col=colors[cell.name])
    }
  }
  predictions500 <- round(predictions[500, ], 3)
  return(predictions500)
}



sc.PlotPCA <- function(counts, colors=NULL, shapes=NULL){
  if (is.null(colors)) {
    colors <- rep("black", ncol(counts))
    names(colors) <- colnames(counts)
  }
  if (is.null(shapes)) {
    shapes <- rep(16, ncol(counts))
    names(shapes) <- colnames(counts)
  }
  cat("PlotPCA\n")
  PCA.results <- PCA(t(counts), scale.unit=T, ncp=4, graph=F)
  par(mfrow=c(2, 2))
  for (i in 1:2){
    for (j in (i+1):3){
      PCi <- PCA.results$ind$coord[, i]
      PCj <- PCA.results$ind$coord[, j]
      plot(PCi, PCj, xlab=i, ylab=j, col=colors[row.names(PCA.results$ind$coord)], 
                                     pch=shapes[row.names(PCA.results$ind$coord)])
    }
  }
  return(PCA.results)
}

sc.GetHKGenes <- function(){
  
  prefix <- "/Users/wm2313/Dropbox (CGC)/SingleCell_T/PilotExperiment/Bioinformatics/"
  
  if(!file.exists(paste0(prefix,"gene.scores.txt"))){
    if (!file.exists(paste0(prefix, "HK_exons.csv"))){
      stop("sc.GetHKGenes error: HK_exons.csv not found")
    }
    cat("Finding high expression housekeeping genes..\t")
    hk.exons <- read.csv(paste0(prefix,"HK_exons.csv"))[, 1:8]
    colnames(hk.exons) <- c("Gene_name", "Refseq", "Refseq_exon_number", "chromosome", 
                            "exon_start", "exon_end", "geometric_mean_rpkm", "std_log2_RPKM")
    library(plyr)
    gene.score <- ddply(hk.exons, "Gene_name", 
                        function(x){c("Score"=
                                        sum(x$geometric_mean_rpkm*(x$exon_end - x$exon_start + 1)) / 
                                        sum(x$exon_end - x$exon_start+1))})
    
    write.table(gene.score, "gene.scores.txt", row.names=F, quote=F)
  } else {
    cat("Reading table of high expression housekeeping genes..\t")
    gene.score <- read.table(paste0(prefix,"gene.scores.txt"), header=T)
  }
  cat("done\n")
  library(plyr)
  gene.score <- arrange(gene.score, Score, decreasing=T)
  return(gene.score)
}


sc.NormalizeData <- function(counts, method=NULL, genes=NULL, colData = NULL, design = NULL){
  cat("Normalizing", method, "..\t")
  if (method=="median"){
    norm.counts <- apply(counts, 2, function(x){x/median(x[x>5])})
  }
  if (method=="edgeR"){
    library(edgeR)
    edgeR.factors <- calcNormFactors(counts)
    norm.counts <- t(t(counts) * edgeR.factors)
  }
  if (method=="DESeq"){
    library(DESeq2)
    sf <- estimateSizeFactorsForMatrix(counts)
    norm.counts <- t(t(counts) / sf)
  } 
  if(method == "vsd"){
    #### construct DESeqDataSet object
    print("please provide phenotype dataframe")
    # dds <- DESeqDataSetFromMatrix(as.matrix(ILN_data), colData=ILN_pheno, design = ~1 + factor(subject) + factor(tissue) + factor(cellType))
    dds <- DESeqDataSetFromMatrix(as.matrix(counts), colData=colData, design = design)
    dds <- estimateSizeFactors(dds)
    dds <- estimateDispersions(dds)
    # vsd <- rlog(dds, blind = FALSE, fast = FALSE)
    vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
    
    # data.sub_subject_norm <- assay(rld)
    data.sub_vsdNorm <- assay(vsd)
    colnames(data.sub_vsdNorm) <- colnames(ILN_data)
  }
  if (method=="sum"){
    norm.counts <- apply(counts, 2, function(x){x/sum(x)})
  }
  if (method=="hkgenes"){
    hk.counts <- counts[genes, ]
    hk.sums <- colSums(hk.counts)
    hk.median.sum <- median(hk.sums)
    # if there are no counts for any of the housekeeping genes, we don't normalize
    norm.counts <- apply(counts, 2, function(x){x*hk.median.sum/
                    ifelse(!sum(x[genes])==0, sum(x[genes]), hk.median.sum)})
  }
  cat("done\n")
  return(norm.counts)
}

JSD_distance_compute <- function(input, cut=40, method = "JSM"){
  ## The input matrix should have min value > 0 
  ## or else Na will be returned
  
  jensen_shannon <- function(p, q){
    ## JSD = H(0.5 *(p +q)) - 0.5H(p) - 0.5H(q)
    # H(X) = \sum x_i * log2(x_i)
    
    p = p / sum(p)
    q = q / sum(q)
    Hj = shannon.entropy(0.5 *(p+q)) 
    Hp = shannon.entropy(p) 
    Hq = shannon.entropy(q)
    
    jsd = Hj - 0.5*(Hp+Hq)
    return(jsd)
  }
  
  
  shannon.entropy <- function(p)
  {
    if (min(p) < 0 || sum(p) <= 0)
      return(NA)
    p.norm <- p[p>0]/sum(p)
    -sum(log2(p.norm)*p.norm)
  }
  
  
  nsample = ncol(input) 
  x =  matrix(nrow=nsample, ncol=nsample)
  z = input[rowMeans(input) > cut, ]
  
  # JSD
  for (i in 1:nsample) {
    for (j in i:nsample) {
      
      if(method == "pearson"){
        
        x[i,j] = cor(z[, i], z[, j], method = "pearson")
        
      }else if(method == "JSM"){
        x[i,j] = (jensen_shannon(z[,i ], z[,j ]))^0.5
      }
      x[j, i] = x[i, j]
      
    }
  }
  
  
  return(x)
  
}


sc.MDS_distanceInput <- function(distances, cells = NULL, txt = F,
                                 dot_size = 1.5, colors = NULL, shapes = 16,
                                 legend = F, bound = 0.1, alpha = 0.5){
  
  if (!is.null(cells)) {
    fit <- cmdscale(distances[cells, cells], eig=TRUE, k=2)
  } else fit <- cmdscale(distances, eig=TRUE, k=2)
  
  x <- fit$points[, 1]
  y <- fit$points[, 2]
  
  v1 = min(x) - bound
  v2 = max(x) + bound
  v3 = min(y) - bound
  v4 = max(y) + bound
  
  if (txt) {
    if (is.null(colors)) {
      colors <- rep("blue", ncol(counts))
    }
    
    plot(x,y, xlim=c(v1,v2), ylim=c(v3,v4), cex.lab=1.5, cex.axis=1.25, pch=shapes,
         col=colors,
         # xlab="MDS Dim 1", ylab="MDS Dim 2",
         cex=dot_size)
    #     abline(h = 0, lty = "dotted")
    #     abline(v = 0, lty = "dotted")
    
    text(x, y, labels=rownames(fit$points), col="black", cex=txt_cex, pos=textpos)
    
  } else {
    if (is.null(colors)) {
      colors <- rep("blue", ncol(counts))
      names(colors) <- colnames(counts)
    }
    
    trans_colors <- colors[rownames(fit$points)]
    for(i in 1: length(trans_colors)){
      trans_colors[i] <- colorRampAlpha(trans_colors[i], n = 1, alpha = alpha)
    }
    
    plot(x,y, xlim=c(v1,v2), ylim=c(v3,v4),cex.lab=2, cex.axis=1.5, pch=shapes -15 ,
         col=colors[rownames(fit$points)],
         xlab="", ylab="",
         cex=dot_size)
    # xlab="", ylab="", cex=dot_size)
    points(x,y, cex.lab=2, cex.axis=1.5, pch=shapes,
           col=trans_colors[rownames(fit$points)],
           cex=dot_size)
    #     abline(h = 0, lty = "dotted")
    #     abline(v = 0, lty = "dotted")
  }
  
  return(fit)
}



sc.MDS <- function (counts, k = 2, colors = NULL, shapes=16, txt=FALSE, dot_size = 2,  alpha = 0.7,  log=T, bound = 0.2, distance = "spearman", textpos = 3, txt_cex = 1 ){
  # coloring scheme:
  # if the labels need to be printed, then the dots are always black, but the colors of txt can be changed
  # if no txt, the dots will be colored by the color vector.
  
  par(mar=c(5, 6, 4, 1)+.1)
  sd.data <- apply(counts, 1, sd)
  counts <- subset(counts, sd.data > 0) # & rowMeans(counts) > 1) ## keep those genes that have variance > 0

  
  if (log) {
    counts <- log(counts+1)
  } 
  
  if (distance == "spearman") {
    sampleDists <- as.dist(1-cor(counts, method="spearman"))
    fit <- cmdscale(sampleDists, eig=TRUE, k=k)
    rownames(fit$points) <- colnames(counts)
  } else if (distance == "euclidean") {
    sampleDists <- as.dist(dist(counts))
    fit <- cmdscale(sampleDists, eig=TRUE, k=k)
    rownames(fit$points) <- colnames(counts)
  } else if (distance == "JSD") {
    sampleDists <- JSD_distance_compute(counts, cut = 10 )
    fit <- cmdscale(sampleDists, eig=TRUE, k=k)
    rownames(fit$points) <- colnames(counts)
  }
  
  
  x <- fit$points[, 1]
  y <- fit$points[, 2]
  
  v1 = min(x)-bound
  v2 = max(x)+bound
  v3 = min(y)-bound
  v4 = max(y)+bound
  
  if (txt) {
    if (is.null(colors)) {
      colors <- rep("blue", ncol(counts))
    }
    
    plot(x,y, xlim=c(v1,v2), ylim=c(v3,v4), cex.lab=1.5, cex.axis=1.25, pch=shapes,
         col=colors,
         xlab="MDS Dim 1", ylab="MDS Dim 2", cex=dot_size)
    #     abline(h = 0, lty = "dotted")
    #     abline(v = 0, lty = "dotted")
    
    text(x, y, labels=rownames(fit$points), col="black", cex=txt_cex, pos=textpos)
    
  } else {
    if (is.null(colors)) {
      colors <- rep("blue", ncol(counts))
      names(colors) <- colnames(counts)
    }
    trans_colors <- colors[rownames(fit$points)]
    for (i in 1: length(trans_colors)) {
      trans_colors[i] <- colorRampAlpha(trans_colors[i], n = 1, alpha = alpha)
    }
    
    plot(x, y, xlim=c(v1, v2), ylim=c(v3, v4), cex.lab = 2, cex.axis = 1.5, pch = shapes - 15 ,
         col = colors[rownames(fit$points)],
         xlab="MDS Dim 1", ylab="MDS Dim 2", cex = dot_size)
    # xlab="", ylab="", cex=dot_size)
    points(x, y, cex.lab = 2, cex.axis = 1.5, pch = shapes,
           col = trans_colors[rownames(fit$points)],
           cex = dot_size)
    #     abline(h = 0, lty = "dotted")
    #     abline(v = 0, lty = "dotted")
  }
  
  return(fit)
}


sc.MDS.lineage <- function(counts, colors=NULL, shapes=NULL, txt=FALSE){
  sampleDists <- as.dist(1-cor(counts+1, 
                               method="spearman"))
  fit <- cmdscale(sampleDists, eig=TRUE, k=3)
  x <- fit$points[, 1]
  y <- fit$points[, 2]
  z <- fit$points[, 3]
  x_bound1 <- min(x)-0.2
  x_bound2 <- max(x)+0.1
  
  if (is.null(colors)) {
    colors <- rep("black", ncol(counts))
    names(colors) <- colnames(counts)
  }
  if (is.null(shapes)) {
    shapes <- rep(16, ncol(counts))
    names(shapes) <- colnames(counts)
  }
  plot(x,y, xlim=c(x_bound1, x_bound2), cex.lab=1.25, cex.axis=1.25, pch=shapes[rownames(fit$points)],
       col="blue",
       xlab="MDS Dim 1", ylab="MDS Dim 2", cex=1.2)
  if(txt) text(x,y,labels=rownames(fit$points), cex=1.2)
  abline(h=0, lty = "dotted")
  abline(v=0, lty = "dotted")
  
  pdf("MDS.pc13.pdf",5,5)
  plot(x,z, xlim=c(x_bound1, x_bound2), cex.lab=1.25, cex.axis=1.25, pch=shapes[rownames(fit$points)],
       col="blue",
       xlab="MDS Dim 1", ylab="MDS Dim 3", cex=1.25)
  if(txt) text(x,z,labels=rownames(fit$points))
  abline(h=0, lty = "dotted")
  abline(v=0, lty = "dotted")
  dev.off()
  
  pdf("MDS.pc23.pdf",5,5)
  plot(y,z, xlim=c(x_bound1, x_bound2), cex.lab=1.25, cex.axis=1.25, pch=shapes[rownames(fit$points)],
       col="blue",
       xlab="MDS Dim 2", ylab="MDS Dim 3", cex=1.25)
  if(txt) text(y,z,labels=rownames(fit$points))
  abline(h=0, lty = "dotted")
  abline(v=0, lty = "dotted")
  dev.off()
  
#   pdf("3d.plot.pdf",4,4)
#  library(scatterplot3d)
#  scatterplot3d(fit$points[,1:3], xlab="PC1", ylab="PC2", zlab="PC3", main="Principle component analysis", pch=16, color=colors)
#   text(x,y,z,labels=rownames(fit$points))
#   dev.off()
  
  #browser()
  return(fit)
}

pca <- function (data, log = TRUE, colors = "blue", shape = 16, dot_size = 2, bound = 5, 
                 txt = F, text_size = 1, alpha = 0.7, title = "", prefix = "") {
  if(log){
    data <- log10(data+1)
  }
  
  sd.data <- apply(data,1, sd) 
  data <- subset(data, sd.data > 0) ## keep those genes that have variance > 0
  pca <- prcomp(t(data), retx = TRUE, center = TRUE, scale = TRUE) #, cor=TRUE)
  
  var.vector <- pca$sdev^2
  prop.vector <- signif(var.vector / sum(var.vector), 3) * 100
  
  
  xdata <- pca$x
  v1 = min(pca$x[,1])-bound
  v2 = max(pca$x[,1])+bound
  v3 = min(pca$x[,2])-bound
  v4 = max(pca$x[,2])+bound
  
  v5 = min(pca$x[,3])-bound
  v6 = max(pca$x[,3])+bound
  
  x <- xdata[,1]
  y <- xdata[,2]
  
  trans_colors <- colors
  for(i in 1: length(trans_colors)){
    trans_colors[i] <- colorRampAlpha(trans_colors[i], n = 1, alpha = alpha)
  }
  
  pdf(paste0(prefix, "_PC12PCA.pdf"))
  par(mar=c(5,6,4,1)+.1)
  plot(xdata[,1], xdata[,2], cex = dot_size, xlim=c(v1,v2), ylim=c(v3,v4), cex.lab=2, cex.axis=2,
       xlab=paste0("PC1 ","(", prop.vector[1], "%)"), ylab=paste0("PC2 ", "(",prop.vector[2], "%)"),
       main=title, col = colors[colnames(data)], pch=shape - 15, asp = 1)
  
  points(x, y, cex.lab = 2, cex.axis = 1.5, pch = shape,
         col = trans_colors,
         cex = dot_size)
  # abline(h = 0, lty = "dotted")
  # abline(v = 0, lty = "dotted")
  if (txt) 
    text(xdata[,1], xdata[,2], labels=colnames(data), col="black", asp = T, cex=text_size, pos=2)
  dev.off()
  
  pdf(paste0(prefix, "_PC13PCA.pdf"))
  plot(xdata[,1], xdata[,3],xlim=c(v1,v2), ylim=c(v5,v6), cex.lab=1.25, cex.axis=1.25, xlab=paste0("PC1 ","(", prop.vector[1], "%)"), ylab=paste0("PC3 ", "(",prop.vector[3], "%)"),
       type="p", main=title, col=colors, pch = shape -15 ,  cex = dot_size)
  points(xdata[,1], xdata[,3],xlim=c(v1,v2), ylim=c(v5,v6), cex.lab=1.25, cex.axis=1.25, 
         type="p",  col=trans_colors, pch = shape,  cex = dot_size)
  abline(h = 0, lty = "dotted")
  abline(v = 0, lty = "dotted")
  if (txt)
    text(xdata[,1], xdata[,3], labels=colnames(data), col=colors, asp = T, cex=0.7,pos=3)
  dev.off()
  
  pdf(paste0(prefix, "_PC23PCA.pdf"))
  plot(xdata[,2], xdata[,3],xlim=c(v3,v4), ylim=c(v5,v6), cex = dot_size, cex.lab=1.25, cex.axis=1.25, xlab=paste0("PC2 ","(", prop.vector[2], "%)"), ylab=paste0("PC3 ", "(",prop.vector[3], "%)"),
       type="p", main=title, col=colors, pch=shape -15 )
  points(xdata[, 2], xdata[,3],xlim=c(v3,v4), ylim=c(v5,v6), cex.lab=1.25, cex.axis=1.25, 
         type="p",  col=trans_colors, pch = shape,  cex = dot_size)
  abline(h = 0, lty = "dotted")
  abline(v = 0, lty = "dotted")
  if (txt)
    text(xdata[,2], xdata[,3], labels=colnames(data), col=colors, asp = T, cex=0.7, pos=3)
  dev.off()
  
  return(pca)
}



hcluster <- function (data, method = "complete", cex = 0.6, colors = NULL, optimize = FALSE) {
  
  distances <- dist(t(data) ) 
  hier.clust <- hclust( distances, method=method )
  
  if(optimize){
    hclustfunc <- function(x){
      library(cba)
      hc <-  hclust(x, method = method)
      co <- order.optimal(x, hc$merge)
      
      ho <- hc
      ho$merge <- co$merge
      ho$order <- co$order 
      
      return(ho)
    }
    distfunc <- function(x) dist(x)
    
    # obtain the clusters
    fit <- hclustfunc(distfunc(t(data)))
    
    hier.clust <- fit
  } else {
    
    if (!is.null(colors)) {
      library(dendextend)
      
      labelCol <- function(x) {
        if (is.leaf(x)) {
          ## fetch label
          label <- attr(x, "label")
          code <- substr(label, 1, 1)
          ## use the following line to reset the label to one letter code
          # attr(x, "label") <- code
          attr(x, "nodePar") <- list(lab.col=colors)
        }
        return(x)
        
      }
      
      d <- dendrapply(as.dendrogram(hier.clust), labelCol)
      plot(d)
      
    } else {
      plot(hier.clust, cex=cex, col=colors)
    }
  }
  return(hier.clust)
}

gcluster <- function(data, cutoff=40, method="pearson") {
  
  jensen_shannon <- function(p, q){
    p = p / sum(p)
    q = q / sum(q)
    Hj = shannon.entropy(0.5 *(p+q)) 
    Hp = shannon.entropy(p) 
    Hq = shannon.entropy(q)
    
    jsd = Hj - 0.5*(Hp+Hq)
    return(jsd)
  }
  
  
  shannon.entropy <- function(p)
  {
    if (min(p) < 0 || sum(p) <= 0)
      return(NA)
    p.norm <- p[p>0]/sum(p)
    -sum(log2(p.norm)*p.norm)
  }
  
  library(gplots)
  nsample = ncol(data) 
  x =  matrix(nrow=nsample, ncol=nsample, dimnames=list(colnames(data), colnames(data)))
  
  # only take rows with mean count > cutoff
  z = data[rowMeans(data) > cutoff, ]
  
  # JSD
  for (i in 1:nsample) {
    for (j in i:nsample) {
      
      if(method == "pearson"){
        
        x[i,j] = cor(z[i], z[j], method = "pearson")
        
      }else if(method == "JSM"){
        
        x[i,j] = (jensen_shannon(z[i], z[j]))^0.5
      }
      x[j, i] = x[i, j]
      
    }
  }
  
  #   pdf(out, 6, 6)
  #  par(family="Menlo")
  par(mar = c(5, 4, 4, 2) + 0.1 + c(15, 0, 0, 15))
  
  #   RowSideColors <- c(rep("black", 28), rep("blue",10)) ## it encodes the color according to row positions
  #   print(length(RowSideColors))
  print(nrow(x))
  my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
  heatmap.2(x,sepwidth=c(0.05,0.05),col=my_palette,margins=c(8,8),scale="none", key=T, keysize=1,trace="none", dendrogram="both", density.info = "none",cexRow=0.9,cexCol=0.9, key.xlab=method,srtCol = NULL,labCol=NULL)
  #  heatmap.2(x, scale="none", key=T, trace="none", dendrogram="row", density.info = "none" , keysize= 5)
  ## heatmap.2(t(log(as.matrix(data))), margins=c(5,10),cexRow=0.8, cexCol=1,key=T,trace="none",sepcolor='white')
  #   dev.off()
  return(x)
}

MDS <- function (star.counts,htseq.summary) {
  sampleDists <- as.dist(1-cor(log2(star.counts+1), 
                               method="spearman"))
  fit <- cmdscale(sampleDists, eig=TRUE, k=2)
  x <- fit$points[, 1]
  y <- fit$points[, 2]
  shapes <- rep(16, nrow(fit$points))
  names(shapes) <- rownames(fit$points)
  low <- subset(htseq.summary, R_UNIQ_MAP_GENES_star < 100000)$SAMPLE_NAME
  mid <- subset(htseq.summary, R_UNIQ_MAP_GENES_star < 1000000)$SAMPLE_NAME
  mid <- setdiff(mid, low)
  high <- subset(htseq.summary, R_UNIQ_MAP_GENES_star >= 1000000)$SAMPLE_NAME
  
  shapes[low] <- 17
  shapes[mid] <- 15
  shapes[high] <- 16
  
  CD69 <- names(star.counts)[star.counts["CD69", ]>1]
  CCR7 <- names(star.counts)[star.counts["CCR7", ]>1]
  # CD45RA <- names(star.counts)[star.counts["SELL", ]>1]
  
  
  # can we differentiate naive and CM from EM?
  colors <- rep("grey", nrow(fit$points))
  names(colors) <- rownames(fit$points)
  colors[CCR7] <- "red"
  colors[CD69] <- "blue"
  colors[intersect(CCR7, CD69)] <- "purple"
  # pdf("MDS.spearman.pdf", width=9, height=5)
  plot(x,y, cex.lab=1.25, cex.axis=1.25, col=colors[rownames(fit$points)],
       pch=shapes[rownames(fit$points)],
       xlab="MDS Dim 1", ylab="MDS Dim 2", cex=1.5)
  legend("bottomright", 
         legend=c("neither", "CCR7", "CD69", "both"), pch=16, 
         col=c("grey", "red", "blue", "purple"))
}


iso.map <- function (star.counts, cutoff=1, colors="blue", k = 2) {
  sampleDists <- as.dist(1-cor(log10(star.counts+1), 
                               method="spearman"))
  library(vegan)
  iso <- isomap(sampleDists, k=k)
  par(mar=c(0,0,2,0))
  par(oma=c(0,0,2,0))
  
  #   CD69 <- names(star.counts)[star.counts["CD69", ]>cutoff]
  #   CCR7 <- names(star.counts)[star.counts["CCR7", ]>cutoff]
  #   # CD45RA <- names(star.counts)[star.counts["SELL", ]>1]
  #   # can we differentiate naive and CM from EM?
  #   colors <- rep("grey", ncol(star.counts))
  #   names(colors) <- colnames(star.counts)
  #   colors[CCR7] <- "red"
  #   colors[CD69] <- "blue"
  #   colors[intersect(CCR7, CD69)] <- "purple"
  #   
  cell.names=colnames(star.counts)
  #   pdf("isomap.pdf",15,8)
  par(mfrow=c(1,2))
  plot(iso, pch=16, cex=1.25, cex.lab=1.25, cex.axis=1.25, cex.main=1.25, 
       col=colors[cell.names])
#   legend("bottomright", 
#          legend=c("neither", "CCR7", "CD69", "both"), pch=16, 
#          col=c("grey", "red", "blue", "purple"))
  
  plot(iso, type="n")
  ordilabel(iso, labels=cell.names)
  #   mtext(paste("Isomap plot of cell expression data", plot.label), 
  #         outer = TRUE, cex = 1.5)
  #   dev.off()
  
  par(mfrow=c(1,1))
  #   dev.off()
}


computeDiffusionDistMatrix <- function (data, quantile.cutoff = 0.4, log = TRUE) {
  ## compute diffusion matrix with Gaussion smoothing
  n <- ncol(data)
  if(log){
    dist.matrix <- as.matrix(dist(log10(data[, 2:n] +1)))
  } else {
    dist.matrix <- as.matrix(dist(data[, 2:n]))
  }
  ## Gaussian Kernel
  epilon <- var(dist(log10(data[, 2:n] + 1 )))
  dist.diffusion <- exp(-dist.matrix / (2 * epilon)) 
  
  ## only take care of the neighbors
  for(i in 1: nrow(dist.diffusion)){
    dist.diffusion[i,i] <- 0
    cutoff <- quantile(dist.diffusion[i, ], quantile.cutoff) # seems 0.4 is the best, only keep the top 60% similar cells
    dist.diffusion[i, dist.diffusion[i, ] < cutoff ] <- 0
  }
  
  ## make sure cell i to j transition summation probability is 1
  row.sums <- rowSums(dist.diffusion)
  dist.diffusion <- dist.diffusion / row.sums
  sampleDists <-  1 - dist.diffusion
  
  return(sampleDists)
}

computeISOMap.distance <- function (data, quantile.cutoff = 0.4, k = 2) {
  
  ## compute the geodesic distance among samples
  
  sampleDists <-  computeDiffusionDistMatrix(data, quantile.cutoff = quantile.cutoff )
  
  library(vegan)
  iso.distance <- as.matrix(isomapdist(sampleDists, k = k, fragmentedOK = F))
  
  #iso.distance <- as.matrix(isomapdist(sampleDists.zscore, epsilon = -2, fragmentedOK = T))
  return(iso.distance)
}



iso.map.singlePlot <- function (data, quantile.cutoff = 0.4, colors="blue", k = 2) {
  
  ## plot an iso map
  
  sampleDists <-  computeDiffusionDistMatrix(data, quantile.cutoff = quantile.cutoff )

  library(vegan)
  iso <- isomap(sampleDists, k=k)
  par(mar=c(0,0,2,0))
  par(oma=c(0,0,2,0))
  
  cell.names=rownames(data)
  plot(iso, pch=16, cex=1.25, cex.lab=1.25, cex.axis=1.25, cex.main=1.25, 
       col=colors[cell.names])
}



rank.and.normalize.vector <- function (x) {
  x <- (rank(x)-.5)/length(x)
  x <- qnorm(x)
}
sc.rank.and.normalize <- function (x) {
  if( is.vector(x) )
    return( rank.and.normalize.vector(x) )
  if( is.data.frame(x) ) {
    d <- NULL
    for (v in x) {
      if( is.null(d) )
        d <- data.frame( rank.and.normalize(v) )
      else 
        d <- data.frame(d, rank.and.normalize(v))
    }
    names(d) <- names(x)
    return(d)
  }
  stop("Data type not handled")
}

sc.PrintPCAGenes <- function(pca.results){
  dimension.PCA.allgenes <- dimdesc(pca.results, axes=c(1,2,3))
  for (i in 1:3){
    cat(i, "\n")
    dim.i <- as.data.frame(dimension.PCA.allgenes[[i]])
    print(head(dim.i, 15))
  }
}

sc.scdediff <- function(o.ifm, counts, groups, prefix, cells.exclude=NULL, n.cores=4, batch=NULL){
  scdediff.file <- paste(prefix, ".scdediff.txt", sep="")
  # calculate models
  if (!file.exists(scdediff.file)){
    library(scde)
    counts <- counts[rowSums(counts)>0, ]
    counts <- counts[, colSums(counts)>1e4]
    valid.cells <- rownames(o.ifm[o.ifm$corr.a>0, ])
    keep.cells <- setdiff(intersect(valid.cells, colnames(counts)), cells.exclude)
    o.ifm <- o.ifm[keep.cells, ]
    counts <- counts[, keep.cells]
    library(scde)
    o.prior <- scde.expression.prior(models=o.ifm, counts=counts, 
                                     length.out=400, show.plot=F)
    ediff <- scde.expression.difference(o.ifm, counts, o.prior, 
                                        groups=groups[row.names(o.ifm)], n.randomizations=100, 
                                        n.cores=n.cores, batch=batch, verbose=1)
    # scde.browse.diffexp(ediff, o.ifm, counts, o.prior ,groups=groups[row.names(o.ifm)], batch=batch[row.names(o.ifm)])
    #scde.test.gene.expression.difference("CLEC9A", models=o.ifm, counts=counts, prior=o.prior, batch=batch, groups=groups)
    if (!is.null(batch)){
      save(ediff, file=paste0(scdediff.file, ".Rdata"))
      write.table(ediff$batch.adjusted[order(abs(ediff$batch.adjusted$Z), decreasing=T), ], 
                  file=scdediff.file, 
                  row.names=T, col.names=T, sep="\t", quote=F)
    } else {
        write.table(ediff[order(abs(ediff$Z), decreasing=T), ], file=scdediff.file, 
                            row.names=T, col.names=T, sep="\t", quote=F)
    }
  } 
  if (!is.null(batch)){
  
    load(paste0(scdediff.file, ".Rdata"))
    b <- ediff$results[order(abs(ediff$results$Z), decreasing=T), ]
    write.table(b, file=scdediff.file, row.names=T, col.names=T, sep="\t", quote=F)
    browser()
    
    return(ediff$results[order(abs(ediff$results$Z), decreasing=T), ])
  } else {
    ediff <- read.table(scdediff.file, header=T, sep="\t")
    return(ediff)
  }
}

sc.scdeMDSheatmap <- function(scde.model, counts, gene, fit, shapeshe){
  
  library(scde)
  scde.model <- subset(scde.model, corr.a > 0)
  a.o.fpm <- scde.expression.magnitude(scde.model, 
                                     counts=counts[gene, rownames(scde.model)])
  min.a.o.fpm <- min(a.o.fpm[a.o.fpm > -Inf])
  a.o.fpm[a.o.fpm==-Inf] = min.a.o.fpm - 1
  my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
  
 # for (gene in c("SORL1", "ALOX5AP", "PTGS2", "C1orf54", "CLEC9A", "ZEB2", "PRDM1", 
#                 "SLAMF7", "ABCG1", "ITGAL", "IL2RG", "CCR7", "ITGAE", 
#                 "LILRB1", "TNFSF13", "FCER1A", "CD19", "CD244", "SELL", "TCF4", "JUN", "SMAD4", 
#                 "TGIF2", "SOX4", "FOXO3", "MTF1", "TCF4", "PLAUR", "HDAC2")){
    gene.data <- a.o.fpm
    norm.gene.data <- (gene.data - min(gene.data)) / (max(gene.data)-min(gene.data))
    norm.gene.data <- round(norm.gene.data*298)+1
    
  
    x <- fit$points[, 1]
    y <- fit$points[, 2]
    x <- x[names(x)%in% colnames(norm.gene.data)]
    y <- y[names(y) %in% colnames(norm.gene.data)]
    plot(x,y,pch=shapes[names(x)], col=my_palette[norm.gene.data[, names(x)]], main=gene)
  
}

sc.scdeheatmap <- function(scde.model, counts, genes, ColSideColors=NULL, dendrogram='both', Rowv=T, Colv=T){
  library(scde)
  scde.model <- subset(scde.model, corr.a > 0)
  o.fpm <- scde.expression.magnitude(scde.model, 
                                     counts=counts[genes, rownames(scde.model)])
  min.o.fpm <- min(o.fpm[o.fpm > -Inf])
  o.fpm.mod <- apply(o.fpm, 2, function(x){x[x==-Inf] = min.o.fpm - 1; return(x)})
  library(gplots)
  my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
  if (!is.null(ColSideColors)){
    ColSideColors <- ColSideColors[colnames(o.fpm.mod)]
    hm <- heatmap.2(o.fpm.mod, scale="none", col=my_palette, 
                    trace="none", cexRow=0.8, cexCol=0.75, 
                    ColSideColors=ColSideColors, dendrogram=dendrogram, Rowv=Rowv, Colv=Colv) 
  } else {
    hm <- heatmap.2(o.fpm.mod, scale="none", col=my_palette, 
                    trace="none", cexRow=0.8, cexCol=0.75, dendrogram=dendrogram, Rowv=Rowv, Colv=Colv) 
  }
    
}

sc.scdeCorrectedCounts <- function(scde.model, counts){
  
  library(scde)
  scde.model <- subset(scde.model, corr.a > 0)
  cells <- intersect(cells, rownames(scde.model))
  o.fpm <- scde.expression.magnitude(scde.model[cells, ], 
                                     counts=counts[genes, cells])
  min.o.fpm <- min(o.fpm[o.fpm > -Inf])
  o.fpm.mod <- apply(o.fpm, 2, function(x){x[x==-Inf] = min.o.fpm - 1; return(x)})
}


Get.scde.normalized.expression <- function(counts, o.ifm, cells=NULL){
  ## this one is returning log scale read counts
  library(scde)
  if (is.null(cells)){
    cells <- intersect(colnames(counts), row.names(o.ifm))
  }
  o.fpm <- scde.expression.magnitude(o.ifm[cells, ], 
                                     counts=counts[, cells])
}

sc.GetscdeModel <- function(counts, prefix, groups=NULL, n.cores=4){
  o.ifm.file <- paste(prefix, ".scdemodel.Rdata", sep="")
  # calculate models
  if (!file.exists(o.ifm.file)){
    counts <- counts[rowSums(counts)>0, ]
    counts <- counts[, colSums(counts)>1e4]
    
  # should I be excluding spike in transcripts?  
    spikein.rows <- grep("ERCC-", row.names(counts))
    mapping.rows <- grep("^__", row.names(counts))
    if (length(c(spikein.rows, mapping.rows))>0){
      counts <- counts[-c(spikein.rows, mapping.rows), ]
    }
    require(scde)
    o.ifm <- scde.error.models(counts=counts,groups=groups, n.cores=n.cores,
                               threshold.segmentation=T,
                               save.crossfit.plots=F,save.model.plots=F,verbose=1);
    save(o.ifm, file=o.ifm.file)
  } else load(o.ifm.file)
  o.ifm <- o.ifm[row.names(o.ifm) %in% names(counts), ]
  return(o.ifm)
}

sc.scdeMRWDist <- function(scde.model, counts, prefix, n.cores=4){
  spikein.rows <- grep("ERCC-", row.names(counts))
  mapping.rows <- grep("^__", row.names(counts))
  if (length(c(spikein.rows, mapping.rows))>0){
    counts <- counts[-c(spikein.rows, mapping.rows), ]
  }
  distance.file <- paste(prefix, ".scdeMRWdist.txt", sep="")
  if (!file.exists(distance.file)){
    library(scde)
    o.ifm <- subset(scde.model, corr.a > 0)
    cd <- counts[, rownames(o.ifm)]
    cd <- cd[rowSums(cd)>0, ]
    cd <- cd[, colSums(cd)>10]
    o.ifm <- o.ifm[colnames(cd), ]
    cell.names <- colnames(cd)
    names(cell.names) <- cell.names
    o.prior <- scde.expression.prior(models=o.ifm, counts=cd, length.out=400, show.plot=F)
    p.self.fail <- scde.failure.probability(models=o.ifm, counts=cd)
  # reclculate posteriors with the individual posterior modes 
    jp <- scde.posteriors(models=o.ifm,cd,o.prior,return.individual.posterior.modes=T,n.cores=n.cores)
  # find joint posterior modes for each gene - a measure of MLE of group-average expression
    jp$jp.modes <- log(as.numeric(colnames(jp$jp)))[max.col(jp$jp)]
    p.mode.fail <- scde.failure.probability(models=o.ifm,magnitudes=jp$jp.modes)
  # weight matrix
    matw <- 1-sqrt(p.self.fail*sqrt(p.self.fail*p.mode.fail))
  # magnitude matrix (using individual posterior modes here)
    mat <- log10(exp(jp$modes)+1);
  # weighted distance
    library(boot)
    mode.fail.dist <- as.dist(1-do.call(rbind,mclapply(cell.names,function(nam1) {
    unlist(lapply(cell.names,function(nam2) {
      corr(cbind(mat[,nam1],mat[,nam2]),w=sqrt(sqrt(matw[,nam1]*matw[,nam2])))
    }))
    },mc.cores=n.cores)));
    write.table(as.matrix(mode.fail.dist), distance.file, sep="\t", col.names=T, quote=F)
  }
  mode.fail.dist <- read.table(distance.file, header=T, check.names=F)
  return(mode.fail.dist)
}


Get.gene.id.to.name.dict <- function(){
  gene.name.dict <- read.table("/Users/wm2313/Dropbox (CGC)/SingleCell_DC/DCAnalysis/Data/humanEnsemblandERCC92_geneid_genename_dictionary.txt")
  gene.names <- gene.name.dict$V2
  names(gene.names) <- gene.name.dict$V1
  return(gene.names)
}

Get.gene.name.to.id.dict <- function(){
  gene.name.dict <- read.table("/Users/wm2313/Dropbox (CGC)/SingleCell_DC/DCAnalysis/Data/humanEnsemblandERCC92_geneid_genename_dictionary.txt")
  gene.names <- gene.name.dict$V1
  names(gene.names) <- gene.name.dict$V2
  return(gene.names)
}




