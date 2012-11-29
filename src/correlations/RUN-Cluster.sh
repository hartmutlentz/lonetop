#!/bin/sh
# request "bin/sh" as shell for job
#$ -S /bin/sh

### cd to directory where the job was submitted:
cd $SGE_O_WORKDIR

LANG=C
LC_ALL=C
export LANG
export LC_ALL

echo "$JOB_NAME $JOB_ID started at `date`"

### put the things to call your program(s) below this line
/users/tsd/lentz/accessibility/lonetop/src/correlations/MatrixList.py
# Danach: qsub -o $PWD -e $PWD HU-Cluster.sh
### put the things to call your program(s) above this line

echo "done at `date`"
