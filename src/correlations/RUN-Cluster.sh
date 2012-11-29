#!/bin/sh
# request "bin/sh" as shell for job
#$ -S /bin/sh

### cd to directory where the job was submitted:
# on beta-cl
cd $SGE_O_WORKDIR
# on dong:
#cd $PBS_O_WORKDIR

LANG=C
LC_ALL=C
export LANG
export LC_ALL

echo "$JOB_NAME $JOB_ID started at `date`"

### put the things to call your program(s) below this line
#
#/users/tsd/lentz/accessibility/lonetop/src/correlations/MatrixList.py
./MatrixList.py
#
# Danach: qsub -o $PWD -e $PWD HU-Cluster.sh
# oder qsub -l mem=50Mb RUN-Cluster.sh
### put the things to call your program(s) above this line

echo "done at `date`"
