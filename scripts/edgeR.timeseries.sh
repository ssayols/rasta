#!/bin/bash
R=/fsimb/groups/imb-bioinfocf/common-tools/dependencies/R/latest/bin
BASEDIR=/fsimb/groups/imb-bioinfocf/projects/cfb_internal/rasta

echo "${R}/Rscript ${BASEDIR}/edgeR.timeseries.R ${BASEDIR}/results/rasta.bad 1 high_dups_unc" | bsub -n1 -cwd $(pwd) -J high_dups_unc
echo "${R}/Rscript ${BASEDIR}/edgeR.timeseries.R ${BASEDIR}/results/rasta.bad 3 high_dups_cor" | bsub -n1 -cwd $(pwd) -J high_dups_cor
echo "${R}/Rscript ${BASEDIR}/edgeR.timeseries.R ${BASEDIR}/results/rasta.good 1 low_dups_unc" | bsub -n1 -cwd $(pwd) -J low_dups_unc
echo "${R}/Rscript ${BASEDIR}/edgeR.timeseries.R ${BASEDIR}/results/rasta.good 3 low_dups_cor" | bsub -n1 -cwd $(pwd) -J low_dups_cor
