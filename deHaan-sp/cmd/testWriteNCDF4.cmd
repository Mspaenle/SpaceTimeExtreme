#!/bin/bash

#BSUB -q normal
#BSUB -J wNCDF4
#BSUB -o WriteNCDF4.out
#BSUB -e WriteNCDF4.out
#BSUB -W 4:10              # wall-clock time (hrs:mins)
#BSUB -n 1                  # number of tasks in job
#BSUB -R "span[ptile=16]"     


source /usr/share/Modules/init/bash
 
module purge
module load /shared/modules/intel-comp-2013.1
module load /shared/modules/open-mpi-1.6.5
module load /shared/modules/open-mpi-ncdf4.3.1-HDF5par
module load /shared/modules/R

export R_HOME_DIR=/install/shared/intel-soft/math/R/BUILD/R-3.1.0/lib64/R
export R_PROFILE=${R_HOME_DIR}/library/Rmpi/Rprofile

export WORKDIR=/gpfs/users/rchailan/phd/R/SpaceTimeExtreme/deHaan-sp/R/scripts

cd ${WORKDIR}
if [ -f  WriteNCDF4.out ]; then
	rm WriteNCDF4.out
fi
mpiexec -np 1 R --no-save -q < testWriteNCDF.R
