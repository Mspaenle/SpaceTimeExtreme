#!/bin/bash

#BSUB -q normal
#BSUB -J reverse
#BSUB -o reverse.out
#BSUB -e reverse.out
#BSUB -W 100:00              # wall-clock time (hrs:mins)
#BSUB -n 16                  # number of tasks in job
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
ncap2 -O -v -s 'tp=1/fp;' megagol2015a-gol-cleaned2.nc megagol2015a-gol-cleaned2-tp.nc
ncatted -a globwave_name,tp,o,c,peak_period megagol2015a-gol-cleaned2-tp.nc
ncatted -a long_name,tp,o,c,"peak period" megagol2015a-gol-cleaned2-tp.nc
ncatted -a units,tp,o,c,"s" megagol2015a-gol-cleaned2-tp.nc
ncatted -a standard_name,tp,o,c,"sea_surface_wave_peak_period" megagol2015a-gol-cleaned2-tp.nc
ncks -A megagol2015a-gol-cleaned2.nc  megagol2015a-gol-cleaned3.nc
ncks -A megagol2015a-gol-cleaned2-tp.nc  megagol2015a-gol-cleaned3.nc
