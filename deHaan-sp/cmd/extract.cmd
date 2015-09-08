#!/bin/bash

#BSUB -q normal
#BSUB -J sp-main
#BSUB -o sp.out
#BSUB -e sp.out
#BSUB -W 100:00              # wall-clock time (hrs:mins)
#BSUB -n 96                  # number of tasks in job
#BSUB -R "span[ptile=16]"     


source /usr/share/Modules/init/bash
 
module purge
module load /shared/modules/intel-comp-2013.1
module load /shared/modules/open-mpi-1.6.5
module load /shared/modules/intel-softs-io
module load /shared/modules/R

export R_HOME_DIR=/install/shared/intel-soft/math/R/BUILD/R-3.1.0/lib64/R
export R_PROFILE=${R_HOME_DIR}/library/Rmpi/Rprofile

export WORKDIR=/gpfs/users/rchailan/phd/R/SpaceTimeExtreme/deHaan-sp/inputs/ww3

cd ${WORKDIR}
ncks -X 3.05,3.40,42.40,43.10 -X 3.5,5.0,43.20,43.583 -X 3.13,3.60,43.10,43.30 megagol2015a-gol-cleaned.nc megagol2015a-gol-cleaned-coastband.nc
