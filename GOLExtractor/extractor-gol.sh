#!/bin/bash
export INDIR=$PWD/inputs
export ORIGINDIR=$INDIR/megagol-2015a
export WORKDIR=work
export OUTDIR=outputs

export ENVELOP="3.,5.,42.25.,43.60"

# source utilities
source resources/sh/utility.sh

# starting notification
rightnow
log "notice" "STARTING... $d"

# for any file, extract the GOL AREA
for originfile in $(ls $ORIGINDIR); do
	ncatted -O -a coordinates,hs,c,c,"longitude latitude" $ORIGINDIR/$originfile
	ncatted -O -a coordinates,dir,c,c,"longitude latitude" $ORIGINDIR/$originfile
	log $? "COORDINATES Attributes $originfile"
	
	/Users/rchailan/Applications/nco/unstable/bin/ncks -O -X $ENVELOP -v hs,dir $ORIGINDIR/$originfile $OUTDIR/$originfile
	log $? "NCKS -X ENVELOP $originfile"
done



# ncatted -O -a coordinates,hs,c,c,"longitude latitude" out.nc 
# ncks -v hs,longitude,latitude out.nc out2.nc
# ncwa -a time out2.nc out3.nc
# ncks -C -O -x -v time out3.nc out-withouttime.nc
# ncks -X 9.,10.,43.,44. -v hs out-withouttime.nc out2.nc