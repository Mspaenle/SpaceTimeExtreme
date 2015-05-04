#! /bin/bash

source "properties.sh"
source "resources/sh/utility.sh"
rightnow
log "notice" "STARTING... $d"

################
## PARAMETERS ##
################
work="work"
bathy="inputs/sirocco.europe.grd"
envelope="-R3/5/42.25/43.60"
envelopelon="-R3/5"
envelopelat="-R42/44"

equidist="10k"
outfile="outputs/impact"

#####################################
## Create Baseline from BATHY file ##
#####################################
if [ "$HASBASELINE" = false ] ; then
    log "notice" "No baseline found. Create a baseline with mode $BASELINEMODE"
    case $BASELINEMODE in 
    	1 )
			echo "TODO";;
		2 ) 
			log 0 "H0 value is $H0"
			gmt grd2xyz ${bathy} ${envelope} > ${work}/bathy-neg.xyz
			log $? "grd2xyz"
			awk '{ print $3 }' ${work}/bathy-neg.xyz > ${work}/bathy-neg.z
			gmt math ${work}/bathy-neg.z NEG = ${work}/bathy.z
			awk '{ print $1"	"$2 }' ${work}/bathy-neg.xyz > ${work}/bathy-neg.xy
			paste ${work}/bathy-neg.xy ${work}/bathy.z > ${work}/bathy.xyz
			log $? "inverse bathy point"
			rm ${work}/bathy-neg.xy ${work}/bathy-neg.z ${work}/bathy-neg.xyz ${work}/bathy.z
			log $? "clean tmp bathy files"

			isobathmin=$(awk "BEGIN {print $H0-0.1; exit}")
			isobathmax=$(awk "BEGIN {print $H0+0.1; exit}")
			gmt gmtselect ${envelope} -Z${isobathmin}/${isobathmax} ${work}/bathy.xyz > ${work}/isobath.xyz
			log $? "iso-value point"

			envelopelon=$(gmt minmax -I0.001 ${work}/isobath.xyz |awk 'BEGIN {FS="/"} {print $1"/"$2 }' )
			# echo $envelopelon
			gmt greenspline ${work}/isobath.xyz ${envelopelon} -I0.01 -Cn10 -St0.95 -G${work}/baseline_spline 
			# gmt greenspline ${work}/isobath.xyz ${envelopelon} -I0.01 -Cn10 -St0.95 -G${work}/baseline_spline -V
			log $? "greenspline"

			gmt grdtrack ${work}/isobath.xyz -Ar -C0.1/0.01/0.05 -G${bathy} -R${envelope}| awk '$0 ~ /^>/ {print $7}' | sed 's/\// /g' > ${work}/track.xyz
			log $? "grdtrack"

			log $? "baseline created";;
		3 )
			echo "TODO";;
	esac
fi

if [ "${GMT}" = true ]; then

	#GMT PARAMS 
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
	######

	palette="work/shallow-water.cpt"
	gmt makecpt -Z-500/0/50 -Csealand -G-5000/0 -M > $palette
	gmt	grdimage ${bathy} $projection $envelope -C$palette -P -K > ${outfile}.ps
	log $? "grdimage"
	gmt pscoast $projection $envelope -Df -G#d9bb7a -Cwhite -N1/0.2p,#0000FF,solid  -P -K -O >> ${outfile}.ps
	log $? "pscoast"
	gmt	psbasemap $envelope $projection $echelle $rose -Bf0.5a1:longitude:/f0.25a0.5:latitude:/:."bla":WeSn -P -O -K >> ${outfile}.ps
	log $? "psbasemap"
	gmt	grdcontour $bathy $envelope -S $projection -C$contourfile -A60+gwhite+f4 -Wcthinnest,black,solid -Wathinner,black,solid -P -O -K >> ${outfile}.ps
	log $? "grdcontour"

	gmt psxy ${work}/isobath.xyz $projection $envelope -S+0.4c -W1p,red -P -O -K >> ${outfile}.ps
	gmt psxy ${work}/baseline_spline $projection $envelope -S+0.4c -W1p,blue -P -O -K >> ${outfile}.ps
	gmt psxy ${work}/track.xyz $projection $envelope -S+0.4c -W1p,black -P -O -K >> ${outfile}.ps
	log $? "psxy"

	gmt	psscale -D21/9/17.5/0.3 -C$palette -B500:"":/:"Depth(m)": -E -O  >> ${outfile}.ps
	log $? "psscale"

	gmt ps2raster -E$png_resolution -A -Tg -P ${outfile}.ps
	log $? "psraster"
	rm ${outfile}.ps
fi



rightnow
log "notice" "That's all folks !  $d"