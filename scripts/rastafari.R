#############################
##
## RASTA implementation on dupradar
##
#############################
library(ShortRead)		# read bam files
library(rtracklayer)	# read gtf files
library(GenomicRanges)	# set operations with GRanges objects
library(Rcpp)			# call C++ code within R
library(VGAM)			# Zero Truncated Poisson distrib
library(parallel)		# parallel processing of list objects

##
## read input files
##
#BAM <- "/home/sergi/IMB/u-ssayols-dupradaR/dupRadaR/test/MEF_24_mm9.bad.bam"
#GTF <- "/home/sergi/IMB/u-ssayols-dupradaR/dupRadaR/test/genes.gtf"
args <- commandArgs(trailingOnly = TRUE)
BAM <- args[1]
OUT <- args[2]
GTF <- "/fsimb/groups/imb-bioinfocf/projects/cfb_internal/rasta/test/genes.gtf"
STRANDNESS <- FALSE	# FIXME: 'reverse' htseq-count like
CORES <- 4

aln <- readGAlignments(BAM)
gtf <- import.gff(GTF,format="gtf",asRangedData=F,feature.type="exon")
gtf <- unlist(reduce(split(gtf,elementMetadata(gtf)$gene_id)))	# flatten reference GTF file (per gene)

##
## construct table with reads per gene position:
## 
##  1. Read name (text) <- NULL
##  2. Sample 1 read counts (numeric)
##  3. Sample 2 read counts (numeric) <- NULL
##  4. Gene name (text)
##  5. Read start site (numeric)
##  6. Read end site (numeric)
##

# remove multimappers and add the gene name to the uniquely mapping reads
hits <- countOverlaps(aln,gtf,ignore.strand=!STRANDNESS)	# calculate reads mapping to multiple exons
aln  <- aln[hits==1]	# and get rid of them
mm   <- as.data.frame(findOverlaps(aln,gtf,ignore.strand=!STRANDNESS))	# match read <--> exon
mm$geneid  <- names(gtf)[mm$subjectHits]
aln  <- as.data.frame(aln)[,c("seqnames","strand","start","end")]
aln$geneid <- mm$geneid[match(as.integer(rownames(aln)),mm$queryHits)]

# go through all reads and check how many start at the same position (over simplification of duplicates?)
aln <- if(STRANDNESS) aln[order(aln$strand,aln$seqnames,aln$start),] else aln[order(aln$seqnames,aln$start),]

cppFunction("
	IntegerVector readCounts(IntegerVector start, int n) {
	//PRE: start MUST be sorted by CHR-START

		IntegerVector hits(n);
		hits[n-1] = 1;	// last read, we're not gonna check if it has duplicates, so at least it has 1 read

		int i = 0;
		while(i < n-1) {

			// check for user interrupt every 1M reads processed
			if(i % 1000000 == 0) Rcpp::checkUserInterrupt();

			// check how many reads share the same start positions
			int j = 1;
			while(i+j <= n && start[i] == start[i+j]) j++;

			// and note the number of reads sharing the same position for ALL these reads
			for(int k=i; k < i+j; k++) hits[k] = j;
			i += j;
		}

		return(hits);
	}")
aln$RC <- readCounts(aln$start,n=nrow(aln))

#readCountsR <- function(start,n) {
#
#	hits = integer(n)
#	i = 1
#
#	while(i < n) {
#		j = 1
#		while(i+j <= n && start[i] == start[i+j]) { j = j + 1 }
#		for(k in i:(i+j-1)) hits[k] = j-1;
#		i = i + j
#	}
#
#	hits
#}
#system.time(readCountsR(sort(aln$start),n=nrow(aln)))

##
## RASTA analysis
##
aln <- split(aln,aln$geneid)	# split the reads by geneid (it will take a while...)
#aln <- mclapply(aln,function(x) x[!duplicated(x$start),],mc.cores=CORES)	# keep only unique start positions
aln <- lapply(aln,function(x) x[!duplicated(x$start),])	# keep only unique start positions
#counts <- mclapply(1:length(aln),function(i) { 
counts <- lapply(1:length(aln),function(i) {
	
#	cat(i,fill=T)	# DEBUG
	x <- aln[[i]]
	read.counts <- x$RC
	uniq.counts <- unique(x$RC)
	read.counts.censored <- ifelse(read.counts > 1,1,read.counts)

	# estimate special cases: only 1 position covered
	if(length(read.counts) == 1) return(c(sum(read.counts),sum(read.counts.censored),sum(read.counts)))

	# estimate special cases: only 1 unique read count (usually no duplicates across the gene)
	if(length(uniq.counts) == 1) {
		r <- try(fit <- vglm(read.counts ~ 1,pospoisson))	# take all reads into the ZTP distrib without cluster distinction
	} else {
		# hierarchical clustering using complete linkage and Canberra distance
		hc <- hclust(dist(uniq.counts,method="canberra"),method="complete")
		ct <- cutree(hc,k=2)

		# zero trucanted poisson distrib calculated on the subcluster with less counts
		cmeans <- c(mean(read.counts[read.counts %in% uniq.counts[ct==1]]),
					mean(read.counts[read.counts %in% uniq.counts[ct==2]])) # **BUG
		min.group <- which(cmeans==min(cmeans))
		r <- try(fit <- vglm(read.counts[read.counts %in% uniq.counts[ct==min.group]] ~ 1, pospoisson))
	}

	pospois.rate <- if(class(r) != "try-error") Coef(fit) else 0
	read.counts.norm <- sapply(read.counts,function(x) { min(x,max(1,qpois(.95,pospois.rate))) })
	
	return(c(sum(read.counts),sum(read.counts.censored),sum(read.counts.norm)))

})
#},mc.cores=CORES)

##
## and plot to compare the read counts before and after correction
##
pdf(paste0("./",OUT,"/",basename(BAM),".pdf"))

x <- matrix(unlist(counts),ncol=3,byrow=T)
rownames(x) <- names(aln)
o <- order(x[,1])
write.csv(x[o,],file=paste0("./",OUT,"/",basename(BAM),".csv"))

plot (1:nrow(x),log2(x[o,1]),type='l',col='red',lwd=5,main=basename(BAM))	# uncorrected
lines(1:nrow(x),log2(x[o,2]),col='green',lwd=.1)	# censored
lines(1:nrow(x),log2(x[o,3]),col='blue',lwd=.1)		# rasta
legend("topleft",fill=c("red","green","blue"),legend=c("no correction","censored","rasta"))

dev.off()
