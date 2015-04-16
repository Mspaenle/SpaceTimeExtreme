#!/bin/bash

# export INDIR=$PWD/inputs
# export ORIGINDIR=$INDIR/megagol-2015a
# export WORKDIR=work
export OUTDIR=/gpfs2/scratch/chailanr/extracted-data-from-megagol2015a
export ENVELOP="3.,5.,42.25.,43.60"
export GRIDDEDDIR=/gpfs2/scratch/chailanr/mirmidon-toolbox/scripts/megagol-autorun/outputs/gridded
export PREFIXOUTFILE=MEGAGOL2015a-OUNF-
export SUFFIXOUTFILE=.nc

# source utilities
source resources/sh/utility.sh

# starting notification
rightnow
log "notice" "STARTING... $d"

ncks --version
log $? "NCO version check"

mkdir $OUTDIR

for year in $(ls $GRIDDEDDIR); do 
	# for any file in a origin-year, extract the GOL AREA
	for originfile in $(ls $GRIDDEDDIR/$year); do

		if [[ $originfile =~ $PREFIXOUTFILE[0-9]*$SUFFIXOUTFILE ]]; then
			ncatted -O -a coordinates,hs,c,c,"longitude latitude"  $GRIDDEDDIR/$year/$originfile
			log $? "COORDINATES Attributes $originfile"
			ncatted -O -a coordinates,dir,c,c,"longitude latitude" $GRIDDEDDIR/$year/$originfile
			log $? "COORDINATES Attributes $originfile"
			ncatted -O -a coordinates,fp,c,c,"longitude latitude"  $GRIDDEDDIR/$year/$originfile
			log $? "COORDINATES Attributes $originfile"
			ncatted -O -a coordinates,t01,c,c,"longitude latitude" $GRIDDEDDIR/$year/$originfile
			log $? "COORDINATES Attributes $originfile"
			ncatted -O -a coordinates,t02,c,c,"longitude latitude" $GRIDDEDDIR/$year/$originfile
			log $? "COORDINATES Attributes $originfile"
			
			mkdir $OUTDIR/$year

			ncks -O -X $ENVELOP -v hs,dir,fp,t01 $GRIDDEDDIR/$year/$originfile $OUTDIR/$year/$originfile
			log $? "NCKS -X ENVELOP $originfile"
		fi
	done
done
