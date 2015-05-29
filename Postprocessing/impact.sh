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

			gmt mapproject ${work}/bathy.xyz $paramJ1 ${envelope} -F -S $offset > ${work}/bathy-lamb93.xyz
			log $? "mapproject bathy"
			gmt xyz2grd ${work}/bathy-lamb93.xyz ${envelopelambert93} -G${work}/bathy-lamb93.grd -I1000/1000
			log $? "xyz2grd"

			isobathmin=$(awk "BEGIN {print $H0-0.1; exit}")
			isobathmax=$(awk "BEGIN {print $H0+0.1; exit}")
			
			gmt gmtselect ${envelopelambert93} -Z${isobathmin}/${isobathmax} ${work}/bathy-lamb93.xyz > ${work}/isobath.xyz
			log $? "iso-value point"

			gmt blockmean ${work}/isobath.xyz ${envelopelambert93} -I${BLOCKMEANINC} > ${work}/isobathblockmean.xyz
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
	         print("Info sur les normales a Baseline : numprofil xi yi xo yo mo po pxstep pystep cosalpha sinalpha") > normalefileout ;
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
             cosalpha = dy/dij;
             sinalpha = dx/dij;
             printf("%d %f %f %f %f %f %f %f %f %f %f\n",numprofile,xi,yi,xo,yo,mo,po,pxstep,pystep,cosalpha,sinalpha) >> normalefileout ;
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
	gmt gshhg ${PATHTOGSHHGBIN} -Ic -N1 > ${input}/coast.dat
	log $? "extract shoreline to ascii file"

	# awk '$1 ~ /^[0-9]/ {print $0}' ${input}/coast.dat > ${work}/coast.xy
	sed '3d' ${input}/coast.dat > ${work}/coast.xy
	log $? "make it continuous file"
	
	gmt mapproject ${work}/coast.xy $paramJ1 ${envelope} -F -S $offset > ${work}/coast-lamb93.xy 
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

	if [ -f $shorelineheight ] ; then 
	    rm $shorelineheight  $distances $slopeprofiles
	fi

	echo "profilenb c1x c1y c2x c2y :" > ${work}/coast-intersections-points-info
	echo "profilenb xc yc" > ${work}/coast-intersections-points-tmp
	log $? "prepare coast-intersections-points file"

	for profile in $(ls ${work}/profiletrack*) ; do
		profilenb=$(echo ${profile} |grep -o '[0-9]\+')
		log "notice" "work on profile ${profilenb}"

		# gmt grdtrack $profile -hi1 -G${work}/bathy-lamb93.grd -sa > ${profile}.xyz
		# log $? "grdtrack profile points"
		sed '1d' ${profile} > ${profile}.xyz
		cat -n ${profile}.xyz > ${profile}-tmp
		rm ${profile}.xyz

		xi=$(awk '$1 == '$profilenb' {print $2}' ${work}/normale)
		yi=$(awk '$1 == '$profilenb' {print $3}' ${work}/normale)
		x0=$(awk '$1 == '$profilenb' {print $4}' ${work}/normale)
		y0=$(awk '$1 == '$profilenb' {print $5}' ${work}/normale)
		cosalpha=$(awk '$1 == '$profilenb' {print $10}' ${work}/normale)
		sinalpha=$(awk '$1 == '$profilenb' {print $11}' ${work}/normale)
		
		
		log $? "gmtselect coast within polygon-i-up"
		# Construction Polygon Pi-low
		x1=$(echo "${xi}-${POLYGONINC}*${sinalpha}" |bc -l)
		y1=$(echo "${yi}-${POLYGONINC}*${cosalpha}" |bc -l)
		echo "${x1} ${y1}" > ${work}/polygon-i.xy
		x1=$(echo "${x0}-${POLYGONINC}*${sinalpha}" |bc -l)
		y1=$(echo "${y0}-${POLYGONINC}*${cosalpha}" |bc -l)
		echo "${x1} ${y1}" >> ${work}/polygon-i.xy
		x1=$(echo "${x0}" |bc -l)
		y1=$(echo "${y0}" |bc -l)
		echo "${x1} ${y1}" >> ${work}/polygon-i.xy
		x1=$(echo "${xi}" |bc -l)
		y1=$(echo "${yi}" |bc -l)
		echo "${x1} ${y1}" >> ${work}/polygon-i.xy
		gmt gmtselect ${work}/coast-lamb93.xy -F${work}/polygon-i.xy > ${work}/tmp_coast-lamb93-low.xy
		log $? "gmtselect coast within polygon-i-low"
		
		awk -v xi=${xi} -v yi=${yi} -v outfile=${work}/c1.xy '
		BEGIN {dist=9999999}
		{
			xc=$1;
			yc=$2;
			dx = xc-xi;
            dy = yc-yi;
            newdist = sqrt(dx*dx+dy*dy);
			if (newdist < dist) {
				xcur=xc;
				ycur=yc;
				dist=newdist;
			}
		}
		END {
			printf("%f %f\n",xcur,ycur) > outfile ;
		}
		' ${work}/tmp_coast-lamb93-low.xy
		# log $? "find c1_i"

		# Construction Polygon Pi-up
		x1=$(echo "${xi}" |bc -l)
		y1=$(echo "${yi}" |bc -l)
		echo "${x1} ${y1}" > ${work}/polygon-i.xy
		x1=$(echo "${x0}" |bc -l)
		y1=$(echo "${y0}" |bc -l)
		echo "${x1} ${y1}" >> ${work}/polygon-i.xy
		x1=$(echo "${x0}+${POLYGONINC}*${sinalpha}" |bc -l)
		y1=$(echo "${y0}+${POLYGONINC}*${cosalpha}" |bc -l)
		echo "${x1} ${y1}" >> ${work}/polygon-i.xy
		x1=$(echo "${xi}+${POLYGONINC}*${sinalpha}" |bc -l)
		y1=$(echo "${yi}+${POLYGONINC}*${cosalpha}" |bc -l)
		echo "${x1} ${y1}" >> ${work}/polygon-i.xy
		# Select only coast points close to the profile
		gmt gmtselect ${work}/coast-lamb93.xy -F${work}/polygon-i.xy > ${work}/tmp_coast-lamb93-up.xy

		awk -v xi=${xi} -v yi=${yi} -v outfile=${work}/c2.xy '
		BEGIN {dist=9999999}
		{
			xc=$1;
			yc=$2;
			dx = xc-xi;
            dy = yc-yi;
            newdist = sqrt(dx*dx+dy*dy);
			if (newdist < dist) {
				xcur=xc;
				ycur=yc;
				dist=newdist;
			}
		}
		END {
			printf("%f %f\n",xcur,ycur) > outfile ;
		}
		' ${work}/tmp_coast-lamb93-up.xy
		# log $? "find c2_i"

		c1x=$(awk '{print $1}' ${work}/c1.xy)
		c1y=$(awk '{print $2}' ${work}/c1.xy)
		c2x=$(awk '{print $1}' ${work}/c2.xy)
		c2y=$(awk '{print $2}' ${work}/c2.xy)

		#find intersection point ci
		a1=$(python -c "a1 = ($y0-$yi)/($x0-$xi); print a1;")
		a2=$(python -c "a2 = ($c2y-$c1y)/($c2x-$c1x); print a2;")
		b1=$(python -c "b1 = $yi - $a1*$xi ; print b1;")
		b2=$(python -c "b2 = $c1y - $a2*$c1x ; print b2;")
		
		xc=$(python -c "xc = ($b2-$b1)/($a1-$a2) ; print xc")
		yc=$(python -c "yc = $a1*($b2-$b1)/($a1-$a2) + $b1 ; print yc")

		distc1c2=$(python -c "import math as math; dist = math.sqrt( ($c2x-$c1x)**2 + ($c2y - $c1y)**2 ) ; print dist")
		delta=$(python -c "import math as math; delta = math.acos( ($c2y - $c1y)/$distc1c2 ) ; print delta")
		xt=$(python -c "import math as math; xt= $xc - 5000*math.cos($delta+2*math.pi); print xt")
		yt=$(python -c "import math as math; yt= $yc + 5000*math.sin($delta+2*math.pi); print yt")
		
		echo "${profilenb} ${xc} ${yc} ${delta} ${xt} ${yt}" >> ${work}/coast-intersections-points-tmp
	done
	
	# awk -v myfilec1=${work}/tmp_c1.xy -v myfilec2=${work}/tmp_c2.xy '{printf("%f %f\n",$2,$3) > myfilec1; printf("%f %f\n",$4,$5) > myfilec2;}' ${work}/coast-intersections-points-info
	# awk  ' {print $2,$3}' ${work}/coast-intersections-points-tmp > ${work}/coast-intersections-points
	# awk  ' {print $5,$6}' ${work}/coast-intersections-points-tmp > ${work}/testpoints

	sort -g -k1 ${work}/coast-intersections-points-tmp > ${work}/coast-intersections-points
	
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

	gmt mapproject $work/nodes.xy $paramJ1 ${envelope} -F -S $offset > $work/nodes-lamb93.xy
	log $? "projection mesh to lamb93"

	for storm in ${storms} ; do #for each storm
		nbstorm=$(echo $storm|grep -o '[0-9]\+')
		log "notice" "work on storm ${nbstorm}"

		for var in ${vars} ; do #for each var
			outfile2=${output}/storm-${nbstorm}-$var

			nbtimesteps=$(ncks -M ${storm} |grep "name = time" | awk -F= '{print $3}' |bc -l)

			touch ${outfile2}
			awk '{print $1}' ${output}/shorelineheights.dat > ${work}/${var}-a

			for t in $(seq 0 $((nbtimesteps-1))); do #for all time step in the given storm
				#process data
				ncks -v $var -d time,$t $flux |awk '! /count/ {print $0}'| grep "time\[$t\]" |grep "node" |awk 'BEGIN {FS="="} {print $4} ' |sed 's/ m//g' |sed 's/_/0.0/g' > $work/$var-$t
				# log $? "extract $var-$t in 1 column"

				paste $work/nodes-lamb93.xy $work/$var-$t > $work/$var-$t.xyz
				# log $? "prepare xyz file"

				gmt xyz2grd $work/$var-$t.xyz  ${envelopelambert93} -I2000/2000 -G$work/$var-$t.grd
				# gmt nearneighbor $work/$var-$t.xyz $envelope -S${NEARNEIGHBORS} -N${NEARNEIGHBORN} -I${NEARNEIGHBORINC}/${NEARNEIGHBORINC} -G$work/$var-$t.grd
				# log $? "xyz2grd"

				gmt grdtrack ${work}/baseline_spline -G$work/$var-$t.grd -sa |awk '{print $3}' > ${work}/${var}-tmp
				log $? "grdtrack var:$var t:$t"
				paste ${work}/${var}-a ${work}/${var}-tmp > ${work}/${var}-b
				mv ${work}/${var}-b ${work}/${var}-a
			done
			mv ${work}/${var}-a ${outfile2}
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
		gmt makecpt -Z-5000/0/50 -Csealand -G-5000/0 -M > $palette		
		# gmt grd2cpt ${work}/bathy-lamb93.grd -Z > $palette
		log $? "grd2cpt"

		# gmt	grdimage ${work}/bathy-lamb93.grd $projX -C$palette -P -K > ${outfile}.ps
		gmt	grdimage ${bathy} ${paramJ2} ${envelope} -C$palette -P -K > ${outfile}.ps
		log $? "grdimage"
		gmt pscoast ${paramJ2} ${envelope} -Df -G#d9bb7a -C#d9bb7a -N1/0.2p,#0000FF,solid  -P -K  -O >> ${outfile}.ps
		log $? "pscoast"
		gmt	psbasemap ${paramJ2} ${envelope} $echelle $rose -Bf0.5a1:longitude:/f0.25a0.5:latitude:/:."bla":WeSn -P -O -K >> ${outfile}.ps
		log $? "psbasemap"
		# gmt	grdcontour ${work}/bathy-lamb93.grd $projX ${envelopelambert93} -S -A60+gwhite+f4 -Wcthinnest,black,solid -Wathinner,black,solid -P -O -K >> ${outfile}.ps
		# log $? "grdcontour"

		# gmt psxy ${work}/isobath.xyz ${projX} ${envelopelambert93} -S+0.4c -W1p,red -P -O -K >> ${outfile}.ps
		gmt psxy ${work}/baseline_spline ${projX} ${envelopelambert93} -S+0.4c -W1p,black -P -O -K >> ${outfile}.ps
		
		# gmt psxy ${work}/coast-lamb93.xy ${projX} ${envelopelambert93} -S+0.2c -W0.2p,black -P -O -K >> ${outfile}.ps
		# gmt psxy inputs/fab-coast-lambert93.dat ${projX} ${envelopelambert93} -S+0.1c -W0.1p,red -P -O -K >> ${outfile}.ps
		# gmt psxy ${work}/tmp_coast-lamb93.xy ${projX} ${envelopelambert93}  -W1p,red -S+0.05c -P -O -K >> ${outfile}.ps

		# gmt psxy -L ${work}/polygon-i.xy ${projX} ${envelopelambert93}  -W1p,green -P -O -K >> ${outfile}.ps
		# gmt psxy ${work}/tmp_c1.xy ${projX} ${envelopelambert93}  -W1p,black -S+0.08c -P -O -K >> ${outfile}.ps
		# gmt psxy ${work}/tmp_c2.xy ${projX} ${envelopelambert93}  -W1p,red -S+0.08c -P -O -K >> ${outfile}.ps

		# gmt psxy ${work}/coast-intersections-points ${projX} ${envelopelambert93}  -W0.5p,red -S+0.1c -P -O -K >> ${outfile}.ps
		# gmt psxy ${work}/testpoints ${projX} ${envelopelambert93}  -W0.5p,black -S+0.1c -P -O -K >> ${outfile}.ps
		

		if [ -f ${work}/isobathblockmean.xyz ] ; then 
			gmt psxy ${work}/isobathblockmean.xyz ${projX} ${envelopelambert93} -S+0.4c -W1p,blue -P -O -K >> ${outfile}.ps
		fi

		for profile in $(ls ${work}/profiletrack*) ; do
			tail -n +2 $profile > ${work}/profiletmp
			gmt psxy -L ${work}/profiletmp ${projX} ${envelopelambert93} -W0.5p,brown -P -O -K >> ${outfile}.ps
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
		# gmt	grdimage ${bathy} $projection $envelope -C$palette -P -K > ${outfile}.ps
		# log $? "grdimage"
		gmt pscoast $projection $envelope -Df -G#d9bb7a -Cwhite -N1/0.2p,#0000FF,solid  -P -K  >> ${outfile}.ps
		log $? "pscoast"
		gmt	psbasemap $envelope $projection $echelle $rose -Bf0.5a1:longitude:/f0.25a0.5:latitude:/:."bla":WeSn -P -O -K >> ${outfile}.ps
		log $? "psbasemap"
		# gmt	grdcontour $bathy $envelope -S $projection -C$contourfile -A60+gwhite+f4 -Wcthinnest,black,solid -Wathinner,black,solid -P -O -K >> ${outfile}.ps
		# log $? "grdcontour"

		# gmt psxy ${work}/isobath.xyz $projection $envelope -S+0.4c -W1p,red -P -O -K >> ${outfile}.ps
		# gmt psxy ${work}/baseline_spline $projection $envelope -S+0.4c -W1p,black -P -O -K >> ${outfile}.ps
		gmt psxy ${input}/coast.dat $projection $envelope -S+0.1c -W1p,black -P -O -K >> ${outfile}.ps
		# gmt psxy ${work}/coast-lamb93-inverted.xy $projection $envelope -S+0.1c -W0.1p,red -P -O -K >> ${outfile}.ps

		# if [ -f ${work}/isobathblockmean.xyz ] ; then 
		# 	gmt psxy ${work}/isobathblockmean.xyz $projection $envelope -S+0.4c -W1p,blue -P -O -K >> ${outfile}.ps
		# fi

		# gmt psxy ${work}/track.xyz $projection $envelope -S+0.4c -W1p,black -P -O -K >> ${outfile}.ps
		# # for profile in $(ls ${work}/profiletrack*) ; do
		# # 	tail -n +2 $profile > ${work}/profiletmp
		# # 	gmt psxy ${work}/profiletmp $projection $envelope -Sp0.1c -W0.5p,brown -P -O -K >> ${outfile}.ps
		# # done
		# log $? "psxy"

		gmt	psscale -D21/9/17.5/0.3 -C$palette -B500:"":/:"Depth(m)": -E -O  >> ${outfile}.ps
		log $? "psscale"

		gmt ps2raster -E$png_resolution -A -Tg -P ${outfile}.ps
		log $? "psraster"
		rm ${outfile}.ps
	fi
fi

rightnow
log "notice" "That's all folks !  $d"