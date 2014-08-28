## isomap analysis on bulk RNA-Seq data


library(vegan)

doIsomap <- function (dist, k=3, col, l)
{
	
	iso <- isomap(dist, k=k)
	pdf("isomap.pdf", width = 30, height = 15)
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


shannon.entropy <- function(p)
{
	if (min(p) < 0 || sum(p) <= 0)
	   return(NA)
	   p.norm <- p[p>0]/sum(p)
	   -sum(log2(p.norm)*p.norm)
}


jsm <- function(p,q){
    # Jensen-Shannon metric
	       ## JSD = H(0.5 *(p +q)) - 0.5H(p) - 0.5H(q)
 	       # H(X) = \sum x_i * log2(x_i)
#	       p = p[p >0 & q >0]
#	       q = q[p>0 & q>0]
	       p = p / sum(p)
	       q = q / sum(q)
	       Hj = shannon.entropy(0.5 *(p+q)) 
	       Hp = shannon.entropy(p) 
	       Hq = shannon.entropy(q)
	       
	       jsd = Hj - 0.5*(Hp+Hq)
#	       cat(Hj, Hp, Hq, jsd, "\n")
	       return(jsd**0.5)
}

distJSM <- function(data, weight) {
	n = ncol(data)
	d <- matrix(nrow = n, ncol = n)
## 	print(n)
	for (i in 1:(n-1)) {
	    d[i, i] = 0
	    for (j in (i+1):n) {
	    	d[i,j] = jsm(as.numeric(data[,i]), as.numeric(data[,j]))
		d[j,i] = d[i,j]
	    }
	}	
	return(d)
}

weightedDist <- function(data, weight) {
	raw.cor <- 1 - cor(data, method = "spearman")
}

plotMDS <- function(distance, colors, cell.labels){
  CLUST <- hclust(distance)
  pdf("mds.pdf", width=10, height=20)
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

#dist = distJSM(counts)

colors = rep("black", ncol(counts))

# positive Sox9
colors[counts["Nkx2-2", ] >= 1] <- "red"
colors[counts["Relb", ] >= 1] <- "blue"
colors[counts["Relb", ] >= 1 & counts["Nkx2-2", ] >= 1 ] <- "purple"

# colors[grep("control", colnames(counts))] <- "blue"
# colors[grep("noniso", colnames(counts))] <- "red"

ll= strsplit(colnames(counts), split='.', fixed=T)
l = colnames(counts)
for (i in 1:ncol(counts)) {
	l[i] = ll[[i]][1]
	
}


plotMDS(as.dist(dist), colors, colnames(counts))



doIsomap(dist, col=colors, k = kn, l = l)
