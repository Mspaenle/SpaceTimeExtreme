#!/bin/bash

#source utilities
source ./resources/sh/utility.sh
rightnow
log "notice" "STARTING... $d"

## PARAMETERS ##
work="work"
output="output"
outfile=$output/sites-info.dat
sites=input/sites.xyz


################
## GMT PARAMS ##
################
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
gmt gmtset PROJ_ELLIPSOID Clarke-1880-IGN
gmt gmtset PROJ_LENGTH_UNIT c
###############

########################
## Projection  params ##
########################
lonREF=3.0
latREF=46.5
latsec1=44.0
latsec2=49.0
offsetX=700000
offsetY=6600000
projL=-Jl${lonREF}/${latREF}/${latsec1}/${latsec2}
paramJ1=${projL}/1:1
paramJ2=-Jl${lonREF}/${latREF}/${latsec1}/${latsec2}/1:850000
offset=-C${offsetX}/${offsetY}
projX=-Jx1:850000
########################

# Set envelope
envelope="-R3/5/42.25/43.60" #envelope considered for fitmaxstab v1

#Basemap param
rose=-T3.13/43.5/0.4i/2
echelle=-L3.13/43.45/34/10k:"Kilometers":

## PROCESS NODE GEOM ##

# map projection
awk -F ";" '{print $1,$2,$3}' ${sites} | sed 's/node\[//' | sed 's/\]//' > ${work}/sites.xyz
log $? "extract xyz file"

gmt mapproject ${work}/sites.xyz $paramJ1 ${envelope} -F -S $offset > ${work}/sites-lamb93.xyz
log $? "mapproject xyz coords (lamb93)"

join -1 3 -2 3 ${work}/sites-lamb93.xyz ${work}/sites.xyz > ${work}/sites-all.xyz
log $? "join lam93 and WGS84 coords"


# geom computations
cp ${work}/sites-all.xyz ${work}/sites-all-bis.xyz
nbsitestot=$(wc -l ${work}/sites-all.xyz | awk '{print $1}')
awk -v pi=3.14159265359 -v outfile=${outfile} -v nbsitestot=${nbsitestot} '
	BEGIN { 
		printf("S1\tS1.x\tS1.y\tS1.lon\tS1.lat\tS2\tS2.x\tS2.y\tS2.lon\tS2.lat\tdist.h\talpha\tlabel\t\n") > outfile ;
	}
	FNR==NR {S1[FNR]=$1; S1x[FNR]=$2; S1y[FNR]=$3; S1lon[FNR]=$4; S1lat[FNR]=$5;next}
	{
		s1=S1[FNR];s1x=S1x[FNR];s1y=S1y[FNR];s1lon=S1lon[FNR];s1lat=S1lat[FNR];
		for (i = FNR+1; i <= nbsitestot; i++) {
			s2=S1[i]; s2x=S1x[i]; s2y=S1y[i]; s2lon=S1lon[i]; s2lat=S1lat[i];
			dx = s1x-s2x;
            dy = s1y-s2y;
            dij = sqrt(dx*dx+dy*dy);
            cosa = (s2x-s1x)/dij;
            sina = (s2y-s1y)/dij;
            atanalpha = atan2(sina/cosa,1);
            aDeg = -(atanalpha * 180/pi - 90);
            

			if (45 < aDeg && aDeg <= 135) {lab="WE";}
			else {lab = "NS";}

			printf("%d\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%s\n",s1,s1x,s1y,s1lon,s1lat,s2,s2x,s2y,s2lon,s2lat,dij,aDeg,lab) > outfile
		}
	}
' ${work}/sites-all.xyz ${work}/sites-all-bis.xyz
# /*if (22.5 < aDeg && aDeg <= 67.5) {lab="NE";}
#             else if (67.5 < aDeg && aDeg <= 112.5) {lab="E";}
#             else if (112.5 < aDeg && aDeg <= 157.5) {lab="SE";}
#             else {lab="N";}*/
	
# if ( 0 < aDeg && aDeg <= 90) {lab = "NE";}
# else {lab = "SE";}


## PLOT ##

# gmt ps2raster -E$png_resolution -A -Tg -P ${outfile}.ps
# log $? "psraster"

# rm ${outfile}.ps
# log $? "clean ps"



