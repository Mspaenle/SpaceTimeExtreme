#!/bin/bash

#source utilities
source ./resources/sh/utility.sh
rightnow
log "notice" "STARTING... $d"

## PARAMETERS ##
work="work"
# bathy="input/sirocco.europe.grd"
output="output"
outfile=${output}"/combined"
#infile="input/megagol2015a-gol.nc"
infile="input/megagol2015-a-gol.nc"

checkdir $work
log $? "$work folder clear"
checkdir $output
log $? "$outputs folder clear"
checkoutfile $outfile
log $? "$outfile clean"

palette="input/shallow-water.cpt"
rm $palette
if [ ! -f $palette ] ; then
	gmt makecpt -Z-500/0/50 -Csealand -G-5000/0 -M > $palette
	log $? "create palette"
fi

png_resolution=800
projection=-JM20c
gmt gmtset PS_MEDIA a3
gmt gmtset MAP_ANNOT_ORTHO ver_text
gmt gmtset MAP_FRAME_TYPE fancy+
gmt gmtset MAP_FRAME_WIDTH 0.08c
gmt gmtset FONT_TITLE Helvetica
gmt gmtset FONT_TITLE 14p
gmt gmtset MAP_TITLE_OFFSET 0.3c
gmt gmtset FONT_LABEL 8p
gmt gmtset MAP_LABEL_OFFSET 0.2c
gmt gmtset FORMAT_FLOAT_OUT %8.8f
gmt gmtset MAP_ANNOT_OBLIQUE 32
## END PARAMETERS ##

# Set envelope
envelope="-R3/5/42.25/43.60" #envelope considered for fitmaxstab v1

#Basemap param
rose=-T3.13/43.5/0.4i/2
echelle=-L3.13/43.45/34/10k:"Kilometers":

## PROCESS NODE EXTRACTION ##
ncks -O -d time,1,1 ${infile} ${work}/foo.nc
log $? "extract 1 time layer"

#lon 
ncks -v longitude ${work}/foo.nc |sed '1,11d'  > $work/longitude.ncks
awk '{print $1"="$2}' $work/longitude.ncks > $work/longitude.ncks2
awk 'BEGIN { FS = "=" } ;{print $1"	"$3}' $work/longitude.ncks2 > $work/longitude
log $? "long extraction"

#lat
ncks -v latitude ${work}/foo.nc |sed '1,11d'  > $work/latitude.ncks
awk '{print $1"="$2}' $work/latitude.ncks > $work/latitude.ncks2
awk 'BEGIN { FS = "=" } ;{print $1"	"$3}' $work/latitude.ncks2 > $work/latitude
log $? "lat extraction"

join $work/longitude $work/latitude > $work/lonlat.xy
awk  '{print $2"	"$3"   "$1}' $work/lonlat.xy > $work/nodes.xy
log $? "join node in a lonlat file"

awk 'NR%20==0' $work/nodes.xy > $work/labeltmp.xy
log $? "extract - only 1 over 20 - label points"
sed 's/node//g' $work/labeltmp.xy > $work/label.xy

## PLOT ##
# gmt pscoast ${projection} ${envelope} -Df -G#d9bb7a -C#d9bb7a -N1/0.2p,#0000FF,solid  -P -K  >> ${outfile}.ps
gmt pscoast ${projection} ${envelope} -Df -G200 -C200 -N1/0.2p,#000000,solid -P -K  >> ${outfile}.ps
log $? "pscoast"

gmt	psbasemap ${projection} ${envelope} ${echelle} ${rose} -Bf0.1a0.2:longitude:/f0.05a0.1:latitude:/:."Node names":WeSn -P -K -O >> ${outfile}.ps
log $? "psbasemap"

 for i in $(seq 1 3); do
 	hyperslabs="input/hyperslabs-"${i}".dat"
 	gmt	psxy ${hyperslabs} ${projection} ${envelope}  -W2,black -L -P -K -O >> ${outfile}.ps
 	log $? "psxy ${hyperslabs}"
 done

gmt psxy ${envelope} ${projection} input/sites.xyz -W1,blue -S+0.2 -P -K -O >> ${outfile}.ps
#gmt pstext ${envelope} ${projection} input/sites.xyz -D0.3/0.3  -F+a39+f6p,Helvetica,black -P -O >> ${outfile}.ps
 gmt psxy ${envelope} ${projection} $work/nodes.xy -S+0.05 -P -O >> ${outfile}.ps
log $? "psxy nodes location"

# gmt pstext ${envelope} ${projection} $work/label.xy  -F+a39+f6p,Helvetica,black -P -O >> ${outfile}.ps
# log $? "pstext labels"

gmt ps2raster -E$png_resolution -A -Tg -P ${outfile}.ps
log $? "psraster"

rm ${outfile}.ps
log $? "clean ps"



