#############################
##
## Venn diagram between the 4 DE analysis
##
#############################
library(VennDiagram)
library(RColorBrewer)

# read input files
BASEDIR = "/fsimb/groups/imb-bioinfocf/projects/cfb_internal/rasta/results/DEanalysis"
files    <- list.files(BASEDIR,pattern="*.csv")
captions <- gsub("\\.csv$","",files)
data <- lapply(files,function(f) read.csv(paste0(BASEDIR,"/",f),head=T,row.names=1))

# calculate DE genes based on pvals
data <- lapply(data,function(x) { x$fdr <- p.adjust(x$PValue,method="BH"); x })
de   <- lapply(data,function(x) rownames(x)[x$fdr < .01])
names(de) <- captions

# venn and save output table with all the genes reported
venn.diagram(de,fill=brewer.pal(4,"Set1"),alpha=.5,filename="venn.tiff")

out <- Reduce(function(x,y) {
	z <- merge(x,y,by=0,all=T,sort=F)
	rownames(z) <- z[,1]
	z[,-1] }, lapply(data,function(x) x[unique(unlist(de)),]))

colnames(out) <- gsub("\\.x$",  paste0(".",captions[1]),colnames(out))
colnames(out) <- gsub("\\.y$",  paste0(".",captions[2]),colnames(out))
colnames(out) <- gsub("\\.x.1$",paste0(".",captions[3]),colnames(out))
colnames(out) <- gsub("\\.y.1$",paste0(".",captions[4]),colnames(out))

write.csv(out[,order(colnames(out))],file="venn.csv")
