#! /bin/bash

source "properties.sh"
source "resources/sh/utility.sh"
rightnow
log "notice" "STARTING... $d"

################
## PARAMETERS ##
################
work="work"
output="outputs"
input="inputs"
bathy="inputs/sirocco.europe.grd"
envelope="-R3/5/42.25/43.60"
envelopelon="-R3/5"
envelopelat="-R42/44"
envelopelambert93="-R700000.0000/865629.1328/6128081.8176/6279561.7280"
outfile="outputs/impact"
################

################
## GMT PARAMS ##
################
png_resolution=800
projection=-JM20c
gmt gmtset PROJ_ELLIPSOID Clarke-1880-IGN
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



#####################################
## Create Baseline from BATHY file ##
#####################################
if [ "${HASBASELINE}" = false ] ; then
    log "notice" "No baseline found. Create a baseline with mode ${BASELINEMODE}"
    case ${BASELINEMODE} in 
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

			gmt mapproject ${work}/bathy.xyz $paramJ1 -Rg -F -S  $offset > ${work}/bathy-lamb93.xyz
			log $? "mapproject bathy"

			isobathmin=$(awk "BEGIN {print $H0-0.1; exit}")
			isobathmax=$(awk "BEGIN {print $H0+0.1; exit}")
			
			# gmt gmtselect ${envelope} -Z${isobathmin}/${isobathmax} ${work}/bathy.xyz > ${work}/isobath.xyz
			envelope=${envelopelambert93}
			gmt gmtselect ${envelope} -Z${isobathmin}/${isobathmax} ${work}/bathy-lamb93.xyz > ${work}/isobath.xyz
			log $? "iso-value point"

			gmt blockmean ${work}/isobath.xyz ${envelope} -I${BLOCKMEANINC} > ${work}/isobathblockmean.xyz
			log $? "blockmean isobath"
			sort ${work}/isobathblockmean.xyz > ${work}/isobathblockmean-sorted.xyz

			ed -s work/isobathblockmean-sorted.xyz <<< $'4m0\nw'
			log "warning" "trick to sort out the isobath points"

			awk '{print NR,$1}' ${work}/isobathblockmean-sorted.xyz > ${work}/isobath.ix
			awk '{print NR,$2}' ${work}/isobathblockmean-sorted.xyz > ${work}/isobath.iy
			max=$(( $(wc -l ${work}/isobathblockmean-sorted.xyz | awk '{print $1}') ))

			min=1
			# nb=200
			nbpersec=10
			# inc=$(echo "$max / $nb" | bc -l)
			inc=$(echo "1 / $nbpersec" | bc -l)
			envelopei=-R${min}/${max}
			gmt greenspline ${work}/isobath.iy ${envelopei} -I${inc} -St0.95 -G${work}/baseline_splineY
			gmt greenspline ${work}/isobath.ix ${envelopei} -I${inc} -St0.95 -G${work}/baseline_splineX
			join ${work}/baseline_splineX ${work}/baseline_splineY | awk '{print $2,$3}' > ${work}/baseline_spline

			# envelopelon=$(gmt minmax -I0.001 ${work}/isobath.xyz |awk 'BEGIN {FS="/"} {print $1"/"$2 }' )
			# envelopelat=$(gmt minmax -I0.001 ${work}/isobath.xyz |awk 'BEGIN {FS="/"} {print "-R"$3"/"$4 }' )

			# minX=$(gmt minmax -I0.001 ${work}/isobath.xyz |awk 'BEGIN {FS="/"} { print $1 }' |sed 's/-R//')
			# maxX=$(gmt minmax -I0.001 ${work}/isobath.xyz |awk 'BEGIN {FS="/"} { print $2 }' )
			# minY=$(gmt minmax -I0.001 ${work}/isobath.xyz |awk 'BEGIN {FS="/"} { print $3 }' )
			# maxY=$(gmt minmax -I0.001 ${work}/isobath.xyz |awk 'BEGIN {FS="/"} { print $4 }'  )

			# nb=100
			# incX=$(echo "($maxX-$minX) / $nb" | bc -l)
			# incY=$(echo "($maxY-$minY) / $nb" | bc -l)
			
			# awk '{print $2,$1}' ${work}/isobath.xyz > ${work}/isobath.yx
			# awk '{print $1,$2}' ${work}/isobath.xyz > ${work}/isobath.xy
			# gmt greenspline ${work}/isobath.xy ${envelopelon} -I${incX} -Cn15 -St0.95 -G${work}/baseline_splineX 
			# gmt greenspline ${work}/isobath.yx ${envelopelat} -I${incY} -Cn15 -St0.95 -G${work}/baseline_splineY 

			# log $? "greenspline"

			# # gmt grdtrack ${work}/isobath.xyz -Ar -C0.1/0.01/0.05 -G${bathy} ${envelope}| awk '$0 ~ /^>/ {print $7}' | sed 's/\// /g' > ${work}/track.xyz
			# # log $? "grdtrack"
			# paste ${work}/baseline_splineX ${work}/baseline_splineY | awk '{print $4,$3}' > ${work}/baseline_spline

			log $? "baseline created";;
		3 )
			echo "TODO";;
	esac
else
	if [ -f ${work}/baseline_spline ]; then
		log "notice" "Baseline found."
	else
		log -1 "Baseline not found."
	fi
fi

##########################
## Create profile files ##
##########################
if [ "${HASPROFILE}" = false ]; then
	log "notice" "No profile found."

	normale="${work}/normale"
	profiletrack="${work}/profiletrack"

	nbrprofiles=$(
		awk -v profiletrack=${profiletrack} -v normale=${normale} -v dist="$ORTHOLENGTH" -v npoints="$ORTHONBPOINTS" '
     		BEGIN { 
     		 numpi  = 3.1415926535897932384626433832795;
             basefileout = profiletrack;
             normalefileout = normale;
             numprofile = 1;
             pstep= dist/npoints;
             getline; xi=$1 ; yi=$2 ;
	         print("Info sur les normales a Baseline : numprofil xi yi xo yo mo po pxstep pystep") > normalefileout ;
           	}
           	{ 
           	 xj=$1 ; yj=$2 ; if (xi == xj && yi == yj) { next; }
             if (xi != xj) { mi = (yj-yi)/(xj-xi); } else { mi = 999999999999999; }
             pi = yi-mi*xi;
             if (yi != yj) { mo = (xi-xj)/(yj-yi); } else { mo = 999999999999999; }
             po = yi-mo*xi;
             dx = xj-xi;
             dy = yj-yi;
             dij = sqrt(dx*dx+dy*dy);
             r = dist/dij;
             xo = xi-r*dy;
             yo = yi+r*dx;
             pxstep = (xo-xi)/npoints;
             pystep = (yo-yi)/npoints;
             printf("%d %f %f %f %f %f %f %f %f\n",numprofile,xi,yi,xo,yo,mo,po,pxstep,pystep) >> normalefileout ;
             myfile= (basefileout numprofile)
             printf("PROFILE %s [ Distance between %d points: %f ]\n",numprofile,npoints,pstep) > myfile
             for (i = 0; i <= npoints; i++) {
                xcur=xi+i*pxstep;
                ycur=yi+i*pystep;
                printf("%f %f\n",xcur,ycur) >> myfile ;
                }
             close(myfile);
             numprofile++;
             xi = xj ; yi = yj ;
            } 
      END { print numprofile}
      ' ${work}/baseline_spline)
	log "notice" "Number of profiles created : $nbrprofiles"

	nbrpprof=`awk 'END { print NR }' ${work}/profiletrack1`
fi

##########################
## Shoreline extraction ##
##########################
if [ "${HASSHORELINE}" = false ]; then
	# gmt gshhg ${PATHTOGSHHGBIN} -I0 -N1 > ${input}/coast.dat
	log $? "extract shoreline to ascii file"

	awk '$1 ~ /^[0-9]/ {print $0}' ${input}/coast.dat > ${work}/coast.xy
	log $? "make it continuous file"

	gmt mapproject ${work}/coast.xy $paramJ1 -Rg -F -S  $offset > ${work}/coast-lambert93.xy
	log $? "mapproject shoreline"	
fi



###########################
## Geometry computations ##
###########################
if [ "${HASGEOMETRY}" = false ]; then

	log "notice" "Begin to work on profile geometry computations"
	shorelineheight=${work}/shorelineheights.dat # file defined as (profile number, highest dune height)
	slopeprofiles=${work}/slopes.dat # file defined as (profile number, slope)
	distances=${work}/distances.dat # file defined as (profile number, distance between (x_baseline, x_coast))

	for profile in $(ls ${work}/profiletrack*) ; do
		profilenb=$(echo ${profile} |grep -o '[0-9]\+')
		log "notice" "work on profile ${profilenb}"

		gmt grdtrack $profile -hi1 -G${bathy} -sa > ${profile}-tmp
		log $? "grdtrack profile points"

		awk '$3 < min {min=$3; minline=$0}; END {print minline}' ${profile}-tmp > ${work}/tmp_coastpoint
	    z=$(awk '{print $3}' ${work}/tmp_coastpoint)
	    echo "$profilenb,$z"
	    if [ "$z" == "" ]; then 
	    	echo "$profilenb	NaN" >> $shorelineheight
			echo "$profilenb	NaN" >> $distances
			echo "$profilenb	NaN" >> $slopeprofiles
		else 
			xd=$(awk '{print $1;}' ${work}/tmp_coastpoint)
			yd=$(awk '{print $2;}' ${work}/tmp_coastpoint)
			xi=$(awk '$1 == '$profilenb' {print $2}' ${work}/normale)
			yi=$(awk '$1 == '$profilenb' {print $3}' ${work}/normale)
			ddune=$(echo "sqrt(($xd - $xi)*($xd - $xi) + ($yd - $yi)*($yd - $yi))" |bc -l)
			slope=$(echo "100*($z+$H0)/$ddune" |bc -l)
			
			echo "$profilenb	$z" >> $shorelineheight
			echo "$profilenb 	$ddune" >> $distances
			echo "$profilenb	$slope" >> $slopeprofiles
	    fi
	done

	sort -g -k1 $shorelineheight > ${output}/shorelineheights.dat
	sort -g -k1 $distances > ${output}/distances.dat
	sort -g -k1 $slopeprofiles > ${output}/slopes.dat
fi

#########################
## Extract storms data ##
#########################
if [ "${HASEXTRACTED}" = false ]; then
	# Interpolation of vars to the (xi,yi) points, i.e. the baselines points
	storms=$(ls ${STORMSDIR}/storm*.nc)
	vars=${VARIABLES}
	# vars="hs_uplifted "
	flux=${storms[0]}

	## PROCESS GRID POINTS FROM WW3-4.18 ounf outputs ##
	#lon 
	ncks -v longitude $flux |awk '! /count/ {print $0}' |awk ' /degree_east/ {print $0}' |grep "node"  > $work/longitude.ncks
	awk '{print $2}' $work/longitude.ncks > $work/longitude.ncks2
	awk 'BEGIN { FS = "=" } ;{print $2}' $work/longitude.ncks2 > $work/longitude
	log $? "long extraction"

	#lat
	ncks -v latitude $flux |awk '! /count/ {print $0}' |awk ' /degree_north/ {print $0}' |grep "node"  > $work/latitude.ncks
	awk '{print $2}' $work/latitude.ncks > $work/latitude.ncks2
	awk 'BEGIN { FS = "=" } ;{print $2}' $work/latitude.ncks2 > $work/latitude
	log $? "lat extraction"

	paste $work/longitude $work/latitude > $work/nodes.xy
	log $? "join node in a lonlat file"

	for storm in ${storms} ; do #for each storm
		nbstorm=$(echo $storm|grep -o '[0-9]\+')
		log "notice" "work on storm ${nbstorm}"

		for var in ${vars} ; do #for each var
			outfile=${output}/storm-${nbstorm}-$var

			nbtimesteps=$(ncks -M ${storm} |grep "name = time" | awk -F= '{print $3}' |bc -l)

			touch ${outfile}
			awk '{print $1}' ${output}/shorelineheights.dat > ${work}/${var}-a

			for t in $(seq 0 $((nbtimesteps-1))); do #for all time step in the given storm
				#process data
				ncks -v $var -d time,$t $flux |awk '! /count/ {print $0}'| grep "time\[$t\]" |grep "node" |awk 'BEGIN {FS="="} {print $4} ' |sed 's/ m//g' |sed 's/_/0.0/g' > $work/$var-$t
				# log $? "extract $var-$t in 1 column"

				paste $work/nodes.xy $work/$var-$t > $work/$var-$t.xyz
				# log $? "prepare xyz file"

				gmt nearneighbor $work/$var-$t.xyz $envelope -S${NEARNEIGHBORS} -N${NEARNEIGHBORN} -I${NEARNEIGHBORINC}/${NEARNEIGHBORINC} -G$work/$var-$t.grd
				# log $? "xyz2grd"

				gmt grdtrack ${work}/baseline_spline -G$work/$var-$t.grd -sa |awk '{print $3}' > ${work}/${var}-tmp
				log $? "grdtrack var:$var t:$t"
				paste ${work}/${var}-a ${work}/${var}-tmp > ${work}/${var}-b
				mv ${work}/${var}-b ${work}/${var}-a
			done
			mv ${work}/${var}-a ${outfile}
		done
	done
fi

#################
## Wave impact ##
#################



################
## GMT PLOTS  ##
################
if [ "${GMT}" = true ]; then

	MAPPROJECT=true
	
	if [ "${MAPPROJECT}" = true ]; then
		palette="work/shallow-water.cpt"
		envelope="-R3/5/42.25/43.60"
		# gmt makecpt -Z-5000/0/50 -Csealand -G-5000/0 -M > $palette		
		# gmt makecpt -Z0/2000/50 -C -G0/2000 -M > $palette
		gmt xyz2grd ${work}/bathy-lamb93.xyz $envelopelambert93 -G${work}/bathy-lamb93.grd -I832/832
		log $? "xyz2grd"
		gmt grd2cpt ${work}/bathy-lamb93.grd -Z > $palette
		log $? "grd2cpt"

		gmt	grdimage -V ${work}/bathy-lamb93.grd $projX -C$palette -P -K > ${outfile}.ps
		log $? "grdimage"
		gmt pscoast -V ${paramJ2} ${envelope} -Df -G#d9bb7a -Cwhite -N1/0.2p,#0000FF,solid  -P -K -O >> ${outfile}.ps
		log $? "pscoast"
		gmt	psbasemap -V ${paramJ2} ${envelope} $echelle $rose -Bf0.5a1:longitude:/f0.25a0.5:latitude:/:."bla":WeSn -P -O -K >> ${outfile}.ps
		log $? "psbasemap"
		gmt	grdcontour -V ${work}/bathy-lamb93.grd $projX ${envelopelambert93} -S -A60+gwhite+f4 -Wcthinnest,black,solid -Wathinner,black,solid -P -O -K >> ${outfile}.ps
		log $? "grdcontour"

		gmt psxy -V ${work}/isobath.xyz ${projX} ${envelopelambert93} -S+0.4c -W1p,red -P -O -K >> ${outfile}.ps
		gmt psxy ${work}/baseline_spline ${projX} ${envelopelambert93} -S+0.4c -W1p,black -P -O -K >> ${outfile}.ps
		# gmt psxy ${work}/coast-lambert93.xy ${projX} ${envelopelambert93} -S+0.4c -W1p,black -P -O -K >> ${outfile}.ps

		if [ -f ${work}/isobathblockmean.xyz ] ; then 
			gmt psxy ${work}/isobathblockmean.xyz ${projX} ${envelopelambert93} -S+0.4c -W1p,blue -P -O -K >> ${outfile}.ps
		fi

		# #gmt psxy ${work}/track.xyz $projection $envelope -S+0.4c -W1p,black -P -O -K >> ${outfile}.ps
		for profile in $(ls ${work}/profiletrack*) ; do
			tail -n +2 $profile > ${work}/profiletmp
			gmt psxy ${work}/profiletmp ${projX} ${envelopelambert93} -Sp0.1c -W0.5p,brown -P -O -K >> ${outfile}.ps
		done
		log $? "psxy"

		gmt	psscale -D21/9/17.5/0.3 -C$palette -B500:"":/:"Depth(m)": -E -O  >> ${outfile}.ps
		log $? "psscale"

		gmt ps2raster -E$png_resolution -A -Tg -P ${outfile}.ps
		log $? "psraster"
		rm ${outfile}.ps
	else
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
		gmt psxy ${work}/baseline_spline $projection $envelope -S+0.4c -W1p,black -P -O -K >> ${outfile}.ps

		if [ -f ${work}/isobathblockmean.xyz ] ; then 
			gmt psxy ${work}/isobathblockmean.xyz $projection $envelope -S+0.4c -W1p,blue -P -O -K >> ${outfile}.ps
		fi

		# gmt psxy ${work}/track.xyz $projection $envelope -S+0.4c -W1p,black -P -O -K >> ${outfile}.ps
		for profile in $(ls ${work}/profiletrack*) ; do
			tail -n +2 $profile > ${work}/profiletmp
			gmt psxy ${work}/profiletmp $projection $envelope -Sp0.1c -W0.5p,brown -P -O -K >> ${outfile}.ps
		done
		log $? "psxy"

		gmt	psscale -D21/9/17.5/0.3 -C$palette -B500:"":/:"Depth(m)": -E -O  >> ${outfile}.ps
		log $? "psscale"

		gmt ps2raster -E$png_resolution -A -Tg -P ${outfile}.ps
		log $? "psraster"
		rm ${outfile}.ps
	fi
fi

rightnow
log "notice" "That's all folks !  $d"