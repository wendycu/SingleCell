options(stringsAsFactors=F)

# variant call in MHC gene
# do heatmap using marker genes and housekeeping genes
# use multiple housekeeping genes to find dead cells
# make sure not clustering on properties of libraries
# correlate between different clusterings
# correlate between Isomap results from this and last time
# double chn
# check number of genes
# compare STAR results to Bowtie results
# sweave report

rank.and.normalize.vector <- function (x) {
  x <- (rank(x)-.5)/length(x)
  x <- qnorm(x)
}
rank.and.normalize <- function (x) {
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



DCAnalysis <- function(){
  for (i in 3:1){
    if (i==1){counts.file <- "dc.counts.txt"
              prefix <- "single_cell_dc_plots_round1/"
              plot.data.label <- "round 1"}
    else if (i==2){counts.file <- "dc.counts.round2.txt"
          prefix <- "single_cell_dc_plots_round2/"
          plot.data.label <- "round 2"}
    else if (i==3){counts.file <- "dc.counts.combined.txt"
                   prefix <- "single_cell_dc_plots_combined/"
                   plot.data.label <- "combined"
    }
  raw.counts <- read.table(counts.file, header=T)
  row.names(raw.counts) <- raw.counts$GENE
  label.counts <- raw.counts[, -1]
  labels <- colnames(label.counts)
  names(labels) <- labels
  colors <- labels
  THBD.expr <- label.counts["THBD", ]
  CD1C.expr <- label.counts["CD1C", ]
  labels[THBD.expr>0] <- "BDCA3"
  labels[CD1C.expr>0] <- "BDCA1"
  labels[CD1C.expr>0 & THBD.expr>0] <- "Both"
  labels[CD1C.expr==0 & THBD.expr==0] <- "None"
  colors[labels%in%"BDCA3"] <- "red"
  colors[labels%in%"BDCA1"] <- "blue"
  colors[labels%in%"Both"] <- "purple"
  colors[labels%in%"None"] <- "black"
  
  dc.genes <- c("CD34", "KIT", "THBD", "BTLA", "CLEC9A", "IRF8", 
                "BATF3", "IDO2", "IDO1", "CD1C", "IRF4", "ETS2", 
                "ID2", "FLT3", "ZBTB46", 
                "SIRPA", "CX3CR1", "CD86", "ITGAX")
  gene.labels <- c("Both", "Both", "(BDCA3) BDCA3", "BDCA3", "BDCA3", "BDCA3", 
                   "BDCA3", "BDCA3", "BDCA3", "(BDCA1) BDCA1", "BDCA1", "BDCA1", 
                   "BDCA1", "(CD135) Both", "Both", 
                   "(CD172a) Both", "Both", "Both", "(CD11c) Both")
  housekeeping.genes <- read.table("housekeeping_genes.txt", header=T)$GENE
  groups=NULL
  cells.exclude=NULL
  if (i==3){
    mds.file <- paste(prefix, "DESeq.spearmanmds_coords.txt", sep="")
    mds <- read.table(mds.file, header=T)
    groups <- row.names(mds)
    names(groups) <- groups
    groups[mds$V1>0] <- "group1"
    groups[mds$V1<=0] <- "group2"
    groups <- as.factor(groups)
    cells.exclude <- c("A3", "B5")
  }
  Analysis(raw.counts, prefix, labels, colors, dc.genes, 
           gene.labels, housekeeping.genes, groups, cells.exclude, 
           plot.data.label=plot.data.label)
}
}

scde.diff <- function(counts, groups, prefix, cells.exclude){
  library(scde)
  counts <- counts[, names(groups)]
  n.cores <- 2
  # should I be excluding spike in transcripts?  
  counts <- counts[-grep("ERCC-", row.names(counts)), ]
  o.ifm.file <- paste(prefix, "o.ifm.groups.Rdata", sep="")
  if (!file.exists(o.ifm.file)){
  o.ifm <- scde.error.models(counts=counts, groups=groups, n.cores=n.cores,
                             threshold.segmentation=T, save.crossfit.plots=F, 
                             save.model.plots=F, verbose=1)
  save(o.ifm, file=o.ifm.file)
  } else load(o.ifm.file)
  counts <- counts[, which(!colnames(counts) %in% cells.exclude)]
  o.ifm <- o.ifm[which(!row.names(o.ifm) %in% cells.exclude), ]
  o.prior <- scde.expression.prior(models=o.ifm, counts=counts, 
                                   length.out=400, show.plot=F)
  ediff <- scde.expression.difference(o.ifm, counts, o.prior, 
                groups=groups[row.names(o.ifm)], n.randomizations=100, n.cores=n.cores, verbose=1)
  head(ediff[order(ediff$Z, decreasing=T), ])
  tail(ediff[order(ediff$Z, decreasing=T), ])
  write.table(round(head(ediff[order(ediff$Z, decreasing=T), ], 20), 4), 
              paste(prefix, "overexpressed.txt", sep=""), sep="\t", row.names=T, quote=F)
  write.table(round(tail(ediff[order(ediff$Z, decreasing=T), ], 20), 4), 
              paste(prefix, "underexpressed.txt", sep=""), sep="\t", row.names=T, quote=F)
  browser()
}


#bdca3.genes <- c("THBD", "BTLA", "CLEC9A", "IRF8", "BATF3", "IDO2", "IDO1")
#bdca1.genes - c("CD1C", "IRF4", "ETS2", "ID2")
# CD34   expressed by all myeloid and lymphoid progenitors
# KIT    (CD117) expressed by all myeloid and lymphoid progenitors
# XCR1
# TLR3
# THBD   (CD141) marker for BDCA3
# BTLA   more expressed in BDCA3
# CLEC9A marker for BDCA3
# IRF8   expressed in BDCA3 (interferon regulatory factor 8)
# BATF3  transcription factor expressed in BDCA3
# in BDC
# IDO1   more expressed in BDCA3


# CD1C   marker for BDCA1
# IRF4   expressed in BDCA1
# ETS2   transcription factor expressed in BDCA1
# ID2    expressed in BDCA1
# FLT3   (CD135)  expressed in both
# ZBTB46 expressed in both
# SIRPA  (CD172a) expressed in both Signal-regulatory protein alpha
# CX3CR1 expressed in both CX3X chemokine receptor 1
# CD86   expressed in both   signals T cells
# ITGAX  (CD11c) integrin, alpha X expressed in both

# CD14  zero expression 


n.cores <- 4

CandidateGeneHeatmap <- function(log.counts, 
                     candidate.genes, candidate.gene.labels, housekeeping.genes, prefix){
  cat("CandidateGeneHeatmap\n")
  housekeeping.gene.subset <- housekeeping.genes[housekeeping.genes %in% row.names(log.counts)]
  housekeeping.gene.counts <- log.counts[c(housekeeping.gene.subset), ]
  housekeeping.gene.mean.ex <- rowMeans(housekeeping.gene.counts)
  housekeeping.genes.hiex <- housekeeping.gene.subset[housekeeping.gene.mean.ex >= 
                                    quantile(housekeeping.gene.mean.ex, 0.5)]
  hiex.housekeeping.genes.counts <- log.counts[housekeeping.genes.hiex, ]
  #
  #heatmap.counts <- log.counts[c(candidate.genes, housekeeping.gene.subset), ]
  #row.names(heatmap.counts) <- c(paste(candidate.genes, candidate.gene.labels), housekeeping.gene.subset)
  my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
  library(gplots)
  # plot housekeeping genes
  pdf(paste(prefix, "housekeeping.genes.pdf", sep=""))
  heatmap.2(hiex.housekeeping.genes.counts, scale="none", col=my_palette, 
            trace="none", cexRow=0.75, cexCol=0.75)
  dev.off()
  # plot candidate genes
  candidate.gene.counts <- log.counts[candidate.genes, ]
  row.names(candidate.gene.counts) <- paste(candidate.genes, candidate.gene.labels, sep="-")
  pdf(paste(prefix, "candidate.genes.pdf", sep=""))
  par(oma=c(1,2,2,6))
  heatmap.2(candidate.gene.counts, scale="none", col=my_palette, 
            trace="none", cexRow=0.75, cexCol=0.75)
  dev.off()
  # no both heatmap
  pdf(paste(prefix, "candidate.genes.no.both.pdf", sep=""))
  par(oma=c(1,2,2,6))
  heatmap.2(candidate.gene.counts[!grepl("Both", row.names(candidate.gene.counts)), ], scale="none", col=my_palette, 
            trace="none", cexRow=0.75, cexCol=0.75)
  dev.off()
  # candidate and housekeeping genes
  gene.counts <- log.counts[c(candidate.genes, housekeeping.genes.hiex), ]
  row.names(gene.counts) <- c(paste(candidate.genes, candidate.gene.labels, sep="-"), housekeeping.genes.hiex)
  pdf(paste(prefix, "candidate.and.housekeeping.genes.pdf", sep=""), height=10, width=10)
  par(oma=c(1,2,2,6))
  heatmap.2(gene.counts, scale="none", col=my_palette, 
            trace="none", cexRow=0.4, cexCol=0.4)
  dev.off()
}

Normalize <- function(counts){
  cat("Normalize\n")
  library(DESeq2)
  sf <- estimateSizeFactorsForMatrix(counts)
  norm.counts <- t(t(counts) / sf)
  return(norm.counts)
}

GetscdeDistances <- function(counts, scde.model, prefix){
  cat("GetscdeDistances\n")
  scde.distances.file <- paste(prefix, "scde.distances.Rdata", sep="")
  # calculate models
  if (!file.exists(scde.distances.file)){
    o.fpm <- scde.expression.magnitude(scde.model, counts=counts);
    require(boot)
    k <- 0.95;
    cell.names <- colnames(counts)
    names(cell.names) <- cell.names
    reciprocal.dist <- as.dist(1 - do.call(rbind, mclapply(cell.names, function(nam1) {
      unlist(lapply(cell.names,function(nam2) {
        # reciprocal probabilities
        f1 <- scde.failure.probability(models=scde.model[nam1,,drop=F],magnitudes=o.fpm[,nam2])
        f2 <- scde.failure.probability(models=scde.model[nam2,,drop=F],magnitudes=o.fpm[,nam1])
        # weight factor
        pnf <- sqrt((1-f1)*(1-f2))*k +(1-k); 
        boot::corr(log10(cbind(counts[,nam1],counts[,nam2])+1),w=pnf)
      }))
    },mc.cores=n.cores)),upper=F)
    save(reciprocal.dist, file=scde.distances.file)
  } else load(scde.distances.file)
  return(reciprocal.dist)
}

GetscdeModel <- function(counts, prefix){
  cat("GetscdeModel\n")
  require(scde)
  o.ifm.file <- paste(prefix, "o.ifm.Rdata", sep="")
  # calculate models
  
  if (!file.exists(o.ifm.file)){
    o.ifm <- scde.error.models(counts=counts,n.cores=n.cores,
                               threshold.segmentation=T,
                               save.crossfit.plots=F,save.model.plots=F,verbose=1);
    save(o.ifm, file=o.ifm.file)
  } else load(o.ifm.file)
  o.ifm <- o.ifm[row.names(o.ifm) %in% names(counts), ]
  o.prior <- scde.expression.prior(models=o.ifm, counts=counts, length.out=400, show.plot=F)
  o.fail.curves <- scde.failure.probability(o.ifm, magnitudes=log((10^o.prior$x)-1))
  pdf(paste(prefix, "Dropoutcurves.pdf", sep=""))
  par(mfrow=c(1,1),mar = c(3.5,3.5,5,0.5), mgp = c(2.0,0.65,0));
  plot(c(),c(), xlim=range(o.prior$x),ylim=c(0,1), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, 
       xlab="expression magnitude (log10)",ylab="drop-out probability", main="Drop-out curves for each cell")
  invisible(apply(o.fail.curves, 2, 
                  function(y) lines(x=o.prior$x, y=y, col="orange")))
  dev.off()
  return(o.ifm)
  #p.self.fail <- scde.failure.probability(models=o.ifm, counts=counts)
}

SampleDistancePlots <- function(sampleDists, prefix, plot.label=NULL, 
                                labels=NULL, colors=NULL, plot.data.label=NULL){
  cat("SampleDistancePlots\n")
  sampleDistMatrix <- as.matrix( sampleDists )
  cell.names <- row.names(sampleDistMatrix)
  if (is.null(colors)){
    colors=rep("black", nrow(sampleDistMatrix))
    names(colors) <- row.names(sampleDistMatrix)
  }
  if (is.null(labels)){
    labels <- row.names(sampleDistMatrix)
    names(labels) <- labels
  } else {
    row.names(sampleDistMatrix) <- paste(row.names(sampleDistMatrix), labels[row.names(sampleDistMatrix)], sep="-")
    colnames(sampleDistMatrix) <- paste(colnames(sampleDistMatrix), labels[row.names(sampleDistMatrix)], sep="-")
  }
  #colnames(sampleDistMatrix) <- NULL
  library( "gplots" )
  library( "RColorBrewer" )
  my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
  pdf(paste(prefix, "sample_heatmap.pdf"))
  heatmap.2( sampleDistMatrix, trace="none", col=my_palette)
  dev.off()
  fit <- cmdscale(sampleDists, eig=TRUE, k=2)
  x <- fit$points[, 1]
  y <- fit$points[, 2]
  write.table(fit$points, file=paste(prefix, "mds_coords.txt", sep=""), 
              row.names=T, quote=F)
  pdf(paste(prefix, "sample_mds.pdf", sep=""))
  plot(x,y, type='n', cex.lab=1.25, cex.axis=1.25, 
       main=paste("MDS plot of cell expression data,\ndistance method ", plot.label, ", ", 
                  plot.data.label, " data", sep=""))
  text(x, y, paste(rownames(fit$points), labels[rownames(fit$points)]), cex=0.5, 
       col=colors[rownames(fit$points)])
  dev.off()
  methods <- c("ward.D", "ward.D2", "complete")
  # add colors of leaves 
  # http://stackoverflow.com/questions/10571266/colouring-branches-in-a-dendrogram-in-r
#   for (method in methods){
#     CLUST <- hclust(sampleDists, method=method)
#     pdf(paste(prefix, ".", method, ".", "HierarchicalClustering.pdf", sep=""))
#     par(cex=0.5)
#     plot(CLUST, xlab="", labels=paste(cell.names, labels[cell.names]), sub="", 
#        axes=F, ylab="", 
#        main=paste("Hierarchical Clustering of cell \nexpression data, method", 
#                   method, plot.label))
#     dev.off()
#   }
  for (method in methods){
    CLUST <- hclust(sampleDists, method=method)
    dhc <- as.dendrogram(CLUST)
    dL <- dendrapply(dhc, function(n){
      if(is.leaf(n)){
        #labelCol <- colors[cell.names];
        #labelCol <- paste("#",substring(digest(attr(n,"label")),1,6), sep="");
        labelCol <- colors[attr(n, "label")]
        #attr(n, "edgePar") <- list(col = labelCol);
        attr(n, "label") <- paste(attr(n, "label"), labels[attr(n, "label")])
        attr(n, "nodePar") <- list(pch = NA, lab.col = labelCol, lab.cex = 0.75);
      }
      n;
    });
    pdf(paste(prefix, ".", method, ".", "HierarchicalClustering.pdf", sep=""))
    par(cex=0.5)
    plot(dL, xlab="", sub="", axes=F, ylab="", cex.main=1.7,
         main=paste("Hierarchical Clustering of counts data,\ndistance method ", 
                    plot.label, ", ", plot.data.label, " data", sep=""))
    dev.off()
    
  }
  library(vegan)
  iso <- isomap(sampleDists, k=3)
  pdf(paste(prefix, "isomap.pdf", sep=""), width=10, height=10)
  par(mfrow=c(1,2))
  #plot(iso, type="n", main="Isomap plot of cell expression data")
  #ordilabel(iso, fill=colors[cell.names], labels=cell.names)
  par(mar=c(0,0,2,0))
  par(oma=c(0,0,2,0))
  plot(iso, pch=16, cex=1.25, cex.lab=1.25, cex.axis=1.25, cex.main=1.25, 
       col=colors[cell.names])
  plot(iso, type="n")
  ordilabel(iso, labels=cell.names)
  mtext(paste("Isomap plot of cell expression data", plot.label), 
        outer = TRUE, cex = 1.5)
  dev.off()
}

SelectVariableGenes <- function(norm.counts){
  cat("SelectVariableGenes\n")
  library(FactoMineR)
  library(plyr)
  library(genefilter)
  num.genes.hiex <- rowSums(norm.counts>10)
  norm.counts <- norm.counts[num.genes.hiex>=3 & num.genes.hiex < 80, ]
  genes.var <- rowVars(norm.counts)
  norm.counts <- norm.counts[genes.var > 0.2, ]
  PCA.results <- PCA(t(norm.counts), scale.unit=T, ncp=4, graph=F)
  dimension.PCA.allgenes <- dimdesc(PCA.results, axes=c(1,2,3,4))
  var.genes <- c()
  for (i in 1:4){
    dim.i <- as.data.frame(dimension.PCA.allgenes[[i]])
    dim.i$gene <- row.names(dim.i)
    pos.corr <- arrange(subset(dim.i, quanti.correlation>0), quanti.p.value)
    neg.corr <- arrange(subset(dim.i, quanti.correlation<0), quanti.p.value)
    var.genes <- c(var.genes, head(pos.corr$gene, 18), head(neg.corr$gene, 18))
  }
  return(unique(var.genes))
}

QuakeAnalysis <- function(){
  counts.file <- "day18.quake.counts.txt"
  raw.counts <- read.table(counts.file, header=T)
  prefix <- "quake_plots/"
  day18.info <- read.delim("day18numreadsinfo.txt", stringsAsFactors=F)
  labels <- day18.info$putative_cell_type
  names(labels) <- day18.info$cell_name
  colors <- labels
  colors[labels%in%"AT2"] <- "red"
  colors[labels%in%"AT1"] <- "green"
  colors[labels%in%"Clara"] <- "blue"
  colors[labels%in%"BP"] <- "purple"
  colors[labels%in%"Ciliated"] <- "orange"
  colors <- as.character(colors)
  Analysis(raw.counts, prefix, labels, colors)
}

PlotPCA <- function(norm.counts, prefix, colors=NULL, threed=F){
  if (is.null(colors)) colors="black"
  cat("PlotPCA\n")
  PCA.results <- PCA(t(norm.counts), scale.unit=T, ncp=4, graph=F)
  pdf(paste(prefix, "PCAplot.pdf", sep=""))
  par(mfrow=c(2, 2))
  for (i in 1:2){
    for (j in (i+1):3){
      PCi <- PCA.results$ind$coord[, i]
      PCj <- PCA.results$ind$coord[, j]
      plot(PCi, PCj, xlab=i, ylab=j, col=colors[row.names(PCA.results$ind$coord)])
    }
  }
  dev.off()
  if (threed){
    library(rgl)
    plot3d(PCA.results$ind$coord[, 1:3])
  browser()
  }
}

FailureProbabilityPlots <- function(scde.model, gene.counts, prefix){
  cat("FailureProbabilityPlots\n")
  o.prior <- scde.expression.prior(models=scde.model,counts=gene.counts,length.out=400,show.plot=F)
  o.fail.curves <- scde.failure.probability(scde.model, magnitudes=log((10^o.prior$x)-1))
  total.gene.counts <- colSums(gene.counts)
  pdf(paste(prefix, "FailureProbabilityPlot.pdf", sep=""))
  par(mfrow=c(2,2))
  plot(total.gene.counts, o.fail.curves[50, ], main=paste("expression", signif(log((10^o.prior$x[50])), 2)))
  plot(total.gene.counts, o.fail.curves[100, ], main=paste("expression", signif(log((10^o.prior$x[100])), 2)))
  plot(total.gene.counts, o.fail.curves[150, ], main=paste("expression", signif(log((10^o.prior$x[150])), 2)))
  plot(total.gene.counts, o.fail.curves[200, ], main=paste("expression", signif(log((10^o.prior$x[200])), 2)))
  dev.off()
}

AnalyzeDropout <- function(counts, exclude.cells=c("E11", "E8")){
  ercc.controls <- read.delim("ERCC_Controls_Analysis.txt", 
                              header=T, check.names=F)
  cell.names <- setdiff(colnames(counts), exclude.cells)
  spike.in.rows <- grep("ERCC-", row.names(counts))
  ercc <- counts[spike.in.rows, ]
  ercc.concs <- ercc.controls[, c(2, 5)]
  ercc[, "ERCC ID"] <- rownames(ercc)
  ercc <- merge(ercc, ercc.concs)
  numcolwise(sum)(ercc)
  ercc <- ercc[, !colnames(ercc) %in% exclude.cells]
  numnonzero <- rowSums(ercc[, cell.names]>0)
  plot(log10(ercc[, "concentration in Mix 2 (attomoles/ul)"]), numnonzero)
}

Analysis <- function(raw.counts, prefix, labels=NULL, colors=NULL, 
                     candidate.genes=NULL, candidate.gene.labels=NULL, 
                     housekeeping.genes=NULL, groups=NULL, cells.exclude=NULL, 
                     plot.data.label=NULL){
  row.names(raw.counts) <- raw.counts$GENE
  controls <- grep("NTC|PTC", colnames(raw.counts))
  counts <- raw.counts[, -c(1, controls)]
  #AnalyzeDropout(counts)
  if (!is.null(groups)){
    scde.diff(counts, groups, prefix, cells.exclude)
  }
  browser()
  #QCPlot(counts, prefix)
  spike.in.rows <- grep("ERCC-", row.names(counts))
  spike.in.counts <- counts[spike.in.rows, ]
  raw.gene.counts <- counts[-spike.in.rows, ]
  total.gene.counts <- colSums(raw.gene.counts)
  total.spike.in.counts <- colSums(spike.in.counts)
  gene.counts <- raw.gene.counts[rowSums(raw.gene.counts) > 0, ]
  gene.counts <- gene.counts[, total.gene.counts > 10000]
  nonzero.genes <- apply(gene.counts, 2, function(x){sum(x>0)})
  gene.counts <- gene.counts[, nonzero.genes>1000]
  
  norm.gene.counts <- as.matrix(Normalize(gene.counts)) 
  
  #rank.pca.counts <- norm.gene.counts
  #threshold <- 10
  #num.above.threshold <- apply(rank.pca.counts, 1, function(x){sum(x>threshold)})
  
  if (!is.null(candidate.genes)){
    CandidateGeneHeatmap(log2(norm.gene.counts+1), 
                         candidate.genes, candidate.gene.labels, 
                         housekeeping.genes, paste(prefix, "log", sep=""))
    #CandidateGeneHeatmap(1 * (norm.gene.counts>0), candidate.genes, candidate.gene.labels, 
    #                     housekeeping.genes, paste(prefix, "binary", sep=""))
  }
  PCA.variable.genes <- SelectVariableGenes(log2(norm.gene.counts+1))
  PlotPCA(log2(norm.gene.counts[PCA.variable.genes, ]+1), prefix, colors)
  GeneVariationPlot(norm.gene.counts, prefix, PCA.variable.genes)
  cor.methods <- c("pearson", "spearman")
  for (cor.method in cor.methods){
    DESeq.sampleDists <- as.dist(1-cor(log2(norm.gene.counts+1), method=cor.method))
    SampleDistancePlots(DESeq.sampleDists, paste(prefix, "DESeq.", cor.method, sep=""), 
                        cor.method, labels, colors, plot.data.label=plot.data.label)
  }
  scde.model <- GetscdeModel(gene.counts, prefix)
  scde.sampleDists <- GetscdeDistances(gene.counts, scde.model, prefix)
  SampleDistancePlots(scde.sampleDists, paste(prefix, "scde", sep=""), "scde", labels, colors)
  FailureProbabilityPlots(scde.model, gene.counts, prefix)
}

GeneVariationPlot <- function(norm.gene.counts, prefix, topVarGenes=NULL){
  cat("GeneVariationPlot\n")
  library( "genefilter" )
  if (is.null(topVarGenes)){
    topVarGenes <- head( order( rowVars( norm.gene.counts ), decreasing=TRUE ), 200 )
  }
  pdf(paste(prefix, "gene.variation.plot.pdf", sep=""))
  library("RColorBrewer")
  heatmap.2( norm.gene.counts[ topVarGenes, ], scale="row",
             trace="none", dendrogram="column",
             col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
  dev.off()
}

QCPlot <- function(counts, prefix){
  cat("QCPlot\n")
  #summary.csv <- read.csv(summary.csv.file)
  ercc.controls <- read.delim("ERCC_Controls_Analysis.txt", 
                              header=T, check.names=F)
  library(plyr)
  spikein.rows <- grep("ERCC-", row.names(counts))
  spikein.counts <- counts[spikein.rows, ]
  raw.gene.counts <- counts[-spikein.rows, ]
  # CV plot
  total.gene.counts <- colSums(raw.gene.counts)
  total.spikein.counts <- colSums(spikein.counts)
  # exclude cells with less than 10,000 counts to genes
  passed.cells <- counts[, total.gene.counts > 10000]
  passed.cells.gene.counts <- passed.cells[-spikein.rows, ]
  passed.cells.spikein.counts <- passed.cells[spikein.rows, ]
  # CV plot
  library(genefilter)
  means <- rowMeans(passed.cells.gene.counts)
  sds <- rowSds(passed.cells.gene.counts)
  ercc.means <- rowMeans(passed.cells.spikein.counts)
  ercc.sds <- rowSds(passed.cells.spikein.counts)
  cvs <- sds/means
  ercc.cvs <- ercc.sds / ercc.means
  pdf(paste(prefix, "CVplot.pdf", sep=""), width=6, height=5)
  par(mar=c(5, 5, 5, 2))
  plot(log10(means), cvs, pch=".", xlim=c(0, 6), xlab="log10(mean)", ylab="log10(CV)", 
       cex.lab=1.5, cex.axis=1.5, cex.main=1.2, main="Coefficient of variation vs. mean \nfor each transcript (spike-in transcripts in blue)")
  points(log10(ercc.means), ercc.cvs, col="blue", pch=16)
  dev.off()
  htseq.summary <- t(apply(counts[, -1], 2, function(x){
    c(TotalCounts=sum(x), 
      TotalERCCCounts=sum(x[spikein.rows]), 
      TotalGeneCounts=sum(x[-spikein.rows]),
      NonzeroGenes=sum(x[-spikein.rows]>0), 
      NonzeroERCC=sum(x[spikein.rows]>0))
  }))
  htseq.summary <- as.data.frame(htseq.summary)
  sample.names <- rownames(htseq.summary)
  htseq.summary$GenePercentage <- htseq.summary$TotalGeneCounts / htseq.summary$TotalCounts 
  htseq.summary[, "Sample Name"] <- rownames(htseq.summary)
  #summary <- merge(htseq.summary, summary.csv)
  #   pdf("GeneMappingPerc.pdf", width=8, height=3.5)
  #   plot(htseq.summary$TotalCounts, htseq.summary$GenePercentage, 
  #        xlab="Total num mapped reads", ylab="Percent mapping to genes", 
  #        cex.lab=1.5, cex.axis=1.5, cex.main=1.5, 
  #        cex.sub=1.5, main="")
  #   outliers <- subset(htseq.summary, GenePercentage < 0.05 | TotalCounts < 500000)
  #   text(outliers$TotalCounts, outliers$GenePercentage, outliers[, "Sample Name"], 
  #        cex=0.5, pos=1)
  #   dev.off()
  spikein.counts$gene <- row.names(spikein.counts)
  ercc.counts <- merge(spikein.counts, ercc.controls, 
                          by.x = "gene", by.y = "ERCC ID")
  #   pdf("ERCCslopes.pdf")
  #   plot(0, 0, type='n', xlim=c(0, 11), ylim=c(0, 13), xlab="log(ERCC mix concentration+1)", 
  #        ylab="log(counts+1)")
  #   for (sa in sample.names){
  #     ercc <- ercc.counts[, sa]
  #     ercc.mix2 <- ercc.counts[, "concentration in Mix 2 (attomoles/ul)"]
  #     fit <- lm(log(ercc+1) ~ log(ercc.mix2+1))
  #     if (coef(fit)[2] < 1.2) {
  #       text(10, coef(fit)[1] + 10*coef(fit)[2], sa)
  #     }
  #     abline(fit)
  #   }
  #   dev.off()
  pdf(paste(prefix, "ERCCscatterplot.pdf", sep=""), height=6, width=4)
  par(mfrow=c(3,2))
  par(mar=c(5, 5, 2, 2))
  par(oma=c(2,2,8,2))
  samples <- sample(colnames(counts), 12)
  samples1 <- samples[1:6]
  samples2 <- samples[7:12]
  for (i in 1:6){
    sample1 <- samples1[i]
    sample2 <- samples2[i]
    plot(log10(ercc.counts[, sample1]+1), log10(ercc.counts[, sample2]+1), 
         xlab=paste("log10 cts,", sample1), 
         ylab=paste("log10 cts,", sample2), cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
    
  }
  mtext("Counts in pairs of samples\n for each spike-in transcript", outer = TRUE, cex = 1)  
  dev.off()
  pdf(paste(prefix, "Nonzerogenes.pdf", sep=""), width=4, height=4)
  par(mar=c(5, 5, 2, 2))
  plot(htseq.summary$TotalGeneCounts, htseq.summary$NonzeroGenes, cex.lab=1.5, cex.axis=1.5, cex.main=1, 
       cex.sub=1.5, main="# genes detected vs. \n# pairs mapping to genes for 67 cells",
       xlab="Total # pairs mapping to genes", ylab="# genes detected")
  dev.off()
  #summary.csv$MappingPercentage <- summary.csv$Number.of.Mapped.Pairs/summary.csv$Number.of.Pairs
#   pdf(paste(prefix, "NumReadPairs.pdf", sep=""),  width=6, height=6)
#   par(mar=c(5, 5, 2, 2))
#   hist(summary.csv$Number.of.Pairs, breaks=20, 
#        xlab="Number of read pairs", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, 
#        cex.sub=1.5, main="Number of read pairs for 67 cells")
#   dev.off()
#   pdf(paste(prefix, "MappingPercentage.pdf", sep=""),  width=6, height=6)
#   par(mar=c(5, 5, 2, 2))
#   hist(summary.csv$MappingPercentage, breaks=20, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, 
#        cex.sub=1.5, 
#        xlab="Mapping %", main="Mapping percentages for all cells")
#   dev.off()
  pdf(paste(prefix, "MappingPercentageSpikeIn.pdf", sep=""), width=6, height=6)
  par(mar=c(5, 5, 2, 2))
  hist(htseq.summary$TotalERCCCounts/htseq.summary$TotalCounts, breaks=20, cex.lab=1.5, cex.axis=1.5, cex.main=1, 
       cex.sub=1.5, xlim=c(0, 1), 
       main="Percentage of mapped pairs mapping to spike-in for all cells", xlab="% mapped pairs mapping to spike-in")
  dev.off()
  pdf(paste(prefix, "NumberPairsMappedtoGenes.pdf", sep=""), width=6, height=6)
  par(mar=c(5, 5, 2, 2))
  hist(htseq.summary$TotalGeneCounts, breaks=20, cex.lab=1.5, cex.axis=1.5, cex.main=1, 
       cex.sub=1.5, 
       main="Number of pairs mapped to genes for all cells", xlab="Number of pairs mapped to genes")
  dev.off()
}
