## isomap analysis on bulk RNA-Seq data


library(vegan)

doIsomap <- function (dist, k=3, col, l)
{
	
	iso <- isomap(dist, k=k)
	pdf("isomap.pdf", width = 20, height = 10)
#	par(mfrow=c(1,1))
	
	par(mfrow=c(1,2))

 	 par(mar=c(0,0,2,0))
 #	 par(oma=c(0,0,2,0))
 	plot(iso, col=colors, pch=1, cex=1.25, cex.lab=1.25, cex.axis=1.25, cex.main=1.25)

  	plot(iso, type="n")
  	ordilabel(iso, labels=l)
  #	mtext(paste("Isomap plot of expression data", plot.label), outer = TRUE, cex = 1.5)

	
#	text(iso$points, l)
	dev.off()
	
}



weightedDist <- function(data, weight) {
	raw.cor <- cor(as.numeric(data))

	
}

plotMDS <- function(distance, colors, cell.labels){
  CLUST <- hclust(distance)
  pdf("mds.pdf")
  par(mfrow=c(2,1))
  par(cex=0.5)
  plot(CLUST, labels=cell.labels, main="", xlab="", sub="", axes=F, ylab="")
  par(cex=1)
  fit <- cmdscale(distance, eig=TRUE, k=2)
  x <- fit$points[, 1]
  y <- fit$points[, 2]
  plot(x,y, type='n')
  text(x,y,rownames(fit$points), cex=0.5, col=colors[rownames(fit$points)])
  dev.off()
}


library(getopt)

spec <- matrix(c(
        'input'     , 'i', 1, "character", "input csv or tsv file (required)",
        'type'     , 't', 2, "character", "csv or tsv (optional); default: tsv",
        "kn", "k", 2, "integer", "[default: 3] k neighbors in Isomap",
        'cellcyle','c',2, "character", "[optional] cell cycle genes", 
        'help'   , 'h', 0, "logical",   "this help"
),ncol=5,byrow=T)

opt = getopt(spec);

if (!is.null(opt$help) || is.null(opt$input)) {
    cat(paste(getopt(spec, usage=T),"\n"));
    q();
}

itype = "t"
if (!is.null(opt$type)) {
   itype = opt$type
   
}

kn = 3
if (!is.null(opt$kn)) {
	kn = opt$kn
}


cfile = opt$input

if( itype == "c")  {
    counts = read.table(cfile, header=T, sep=",")
} else if (itype == "t") {
  counts = read.table(cfile, header=T, sep="\t")                    
} else {
  counts = read.table(cfile, header=T)                      
}

dist = weightedDist(counts)

colors = rep("black", ncol(counts))
# colors[grep("control", colnames(counts))] <- "blue"
# colors[grep("noniso", colnames(counts))] <- "red"

ll= strsplit(colnames(counts), split='.', fixed=T)
l = colnames(counts)
for (i in 1:ncol(counts)) {
	l[i] = ll[[i]][1]
	
}


plotMDS(as.dist(1- dist), colors, colnames(counts))



doIsomap(1 - dist, col=colors, k = kn, l = l)