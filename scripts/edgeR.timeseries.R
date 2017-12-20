#####################################
##
## Script to perform DE on 2 different conditions, based on edgeR Negative Binomial model
## and the count from htseq-count, in a time-series design with 2 time points: 0h-2h.
## 
## Takes 1 argument: which counts to use(1: raw, 2: censored, 3:rasta). The design matrix is hardcoded.
##
## Adapted for the RASTA analysis and Andrea Schaeffer's datasets. The analysis will be
## performed 3 times for raw counts, censored and rasta corrected.
##
######################################
library("edgeR")
library("RColorBrewer")
library("RColorBrewer")
library("gplots")
library("ggplot2")
library("xtable")

##
## get name patterns from command line
##
args   <- commandArgs(trailingOnly=TRUE)
prefix <- "MEF_"
suffix <- "_mm9.bam.csv"
#cpath  <- "/fsimb/groups/imb-bioinfocf/projects/cfb_internal/rasta/results/rasta.good/"
#col    <- 1	# first column: raw counts
cpath  <- args[1]
col    <- as.integer(args[2])	# column to use for gene counts
main   <- args[3]	# basename for output files

pdf(paste0(main,".pdf"))

##
## read count tables from htseq-count
##
files    <- sort(list.files(pattern=paste0(".*",suffix),path=cpath))
samples  <- gsub(paste0(suffix,"$"),"",files)
countTable <- lapply(paste0(cpath,"/",files),read.table,row.names=1,sep=",",head=T)
countTable <- lapply(countTable,function(x) x[,col,drop=F])	# select the counts only from this column
countTable <- Reduce(function(x,y) {
	x <- merge(x,y,by=0,all=T,sort=F)
	rownames(x) <- x[,1]
	x[,-1]
} ,countTable)
colnames(countTable) <- samples
countTable <- apply(countTable,2,function(x) ifelse(is.na(x),0,x))

##
## create the design and contrasts matrix
##
group  <- factor(c(rep("WT_0h",3),rep("Ing1_0h",3),rep("WT_2h",3),rep("Ing1_2h",3)))
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
conts  <- makeContrasts(Ing1vsWT=(Ing1_2h-Ing1_0h)-(WT_2h-WT_0h),levels=design)

##
## Create digital gene expression, estimate dispersions and calculate linear model for DE
##
# DGE and compute effective library size estimation by TMM normalization
y <- DGEList(counts=countTable,group=group)
y <- calcNormFactors(y)

# filter uniformative genes
m <- 1e6 * t(t(y$counts) / y$samples$lib.size)	# per library size and million reads normalized counts matrix
y <- y[rowSums(m > 1) >= 2,] 	# select genes with RPM > 1 in more than 1 sample

# dispersions
y <- estimateGLMCommonDisp(y,design)	# overall dispersion for the dataset
y <- estimateGLMTrendedDisp(y,design)	# gene-wise dispersion estimates
y <- estimateGLMTagwiseDisp(y,design)	# tagwise dipersion, strongly recommended in multifactor designs

# call DE (standard comparison between 2 conditions)
fit <- glmFit(y,design)	# fit model
lrt <- glmLRT(fit,contrast=conts[,"Ing1vsWT"])

##
## sanity checks
##
m <- log(m[rownames(y),] + 1)	# log transform the counts table. And add 1 tag to every gene to avoid log(0) counts
v <- apply(m,1,sd,na.rm=T)		# get top variant genes

# MDS and variance estimation
plotMDS.DGEList(y,cex=.3,col=brewer.pal(length(levels(group)),"Accent")[group])
plotBCV(y)	# variance along log gene count-per-milion

# Heatmap of top variant 50 genes of the counts-per-milion table
hmcol <- colorRampPalette(brewer.pal(length(levels(group)),"Oranges"))(100)
heatmap.2(m[rev(order(v))[1:50],],col=hmcol,trace="none",margin=c(10,6))

# Heatmap of sample to sample distances
dists <- dist(t(m))
mat <- as.matrix(dists)
heatmap.2(mat,trace="none",col=rev(hmcol),margin=c(13,13))

# MA plots
# get DE genes (p.adjust='BH', pval<.05)
de <- decideTestsDGE(lrt)
degenes <- rownames(y)[as.logical(de)]

# MA plot
plotSmear(lrt,de.tags=degenes,main=main)	# MA plot
abline(h=c(-1,1),col="blue")	# indicate 2-fold changes in the MA plot
abline(v=0,col="blue")			# indicate >1 counts-per-million

# save the DE results into a file (whole gene list in a csv file, 50 top in a tex file)
write.csv(lrt$table,file=paste0(main,".csv"))
f <- file(paste0(main,".tex"))
writeLines(print(xtable(as.data.frame(topTags(lrt,n=50)))),f)
close(f)

dev.off()
