###############################
##
## call dupradar (R version) to make duprate plots
##
###############################

library(dupRadar)

args <- commandArgs(trailingOnly = TRUE)
BAM <- args[1]
OUT <- args[2]
GTF <- "/fsimb/groups/imb-bioinfocf/projects/cfb_internal/rasta/test/genes.gtf"
BAMUTIL  <- "/fsimb/groups/imb-bioinfocf/common-tools/dependencies/bamUtil/latest/bamUtil/bin"
STRANDED <- 2
PAIRED   <- FALSE
THREADS  <- 4

pdf(paste0("./",OUT,"/",basename(BAM),".dupradar.pdf"))
bamDuprm <- markDuplicates(dupremover="bamutil",BAM,path=BAMUTIL,rminput=FALSE)
dm <- analyzeDuprates(bamDuprm,GTF,STRANDED,PAIRED,THREADS)
duprateExpPlot(DupMat=dm)
dev.off()
