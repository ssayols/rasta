#!/bin/bash
R=/fsimb/groups/imb-bioinfocf/common-tools/dependencies/R/latest/bin
BASEDIR=/fsimb/groups/imb-bioinfocf/projects/cfb_internal/rasta
for f in ./mapped/good/*.bam; do echo "${R}/Rscript ${BASEDIR}/dupradar.R ${BASEDIR}/${f} ./test/good" | bsub -n4 -cwd $(pwd) -J $(basename $f); done
for f in ./mapped/bad/*.bam; do echo "${R}/Rscript ${BASEDIR}/dupradar.R ${BASEDIR}/${f} ./test/bad" | bsub -n4 -cwd $(pwd) -J $(basename $f); done
