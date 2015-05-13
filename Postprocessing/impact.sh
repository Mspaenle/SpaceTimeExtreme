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

outfile="outputs/impact"

# GMT PARAMS #
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

			isobathmin=$(awk "BEGIN {print $H0-0.1; exit}")
			isobathmax=$(awk "BEGIN {print $H0+0.1; exit}")
			gmt gmtselect ${envelope} -Z${isobathmin}/${isobathmax} ${work}/bathy.xyz > ${work}/isobath.xyz
			log $? "iso-value point"

			envelopelon=$(gmt minmax -I0.001 ${work}/isobath.xyz |awk 'BEGIN {FS="/"} {print $1"/"$2 }' )
			envelopelat=$(gmt minmax -I0.001 ${work}/isobath.xyz |awk 'BEGIN {FS="/"} {print "-R"$3"/"$4 }' )
			# awk '{print $2,$1}' ${work}/isobath.xyz > ${work}/isobath.yx
			# awk '{print $1,$2}' ${work}/isobath.xyz > ${work}/isobath.xy
			# gmt greenspline ${work}/isobath.xy ${envelopelon} -I${EQUIDIST} -Cn15 -St0.95 -G${work}/baseline_splineX 
			# gmt greenspline ${work}/isobath.yx ${envelopelat} -I${EQUIDIST} -Cn15 -St0.95 -G${work}/baseline_splineY 

			awk '{print NR,$1}' ${work}/isobath.xyz > ${work}/isobath.ix
			awk '{print NR,$2}' ${work}/isobath.xyz > ${work}/isobath.iy
			max=82
			min=1
			nb=100
			inc=$(echo "($max-$min) / $nb" | bc -l)
			# envelopei=$(gmt minmax -I0.1 ${work}/isobath.ix |awk 'BEGIN {FS="/"} {print $1"/"$2 }' )
			envelopei=-R1/82
			gmt greenspline ${work}/isobath.iy ${envelopei} -I${inc} -Cn5 -St0.95 -G${work}/baseline_splineY 
			gmt greenspline ${work}/isobath.ix ${envelopei} -I${inc} -Cn5 -St0.95 -G${work}/baseline_splineX 

			# gmt greenspline ${work}/isobath.xyz ${envelopelon} -I0.01 -Cn10 -St0.95 -G${work}/baseline_spline -V
			log $? "greenspline"

			# gmt grdtrack ${work}/isobath.xyz -Ar -C0.1/0.01/0.05 -G${bathy} ${envelope}| awk '$0 ~ /^>/ {print $7}' | sed 's/\// /g' > ${work}/track.xyz
			# log $? "grdtrack"
			paste ${work}/baseline_splineX  ${work}/baseline_splineY | awk '{print $2,$4}' > ${work}/baseline_spline
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
	# gmt grdtrack ${work}/baseline_spline -C${PROFILELENGTH}/${EQUIPROFILEDS}/${EQUIPROFILESPACE} -Ar  -G${bathy} ${envelope} > ${work}/profiles
	# log $? "grdtrack"

	# cat  ${work}/profiles |awk '$1 ~ /^>/ { print $0 }' > ${work}/profilesheader
	# log $? "save profiles header"

	# echo ">"  >> ${work}/profiles
	# countfull=$(wc -l ${work}/profilesheader |awk '{print $1}')
	# count=$(( ${countfull} - 1 ))
	# csplit -s -n 5 -f ${work}/profiletrack  ${work}/profiles "/>/" "{${count}}" 
	# log $? "csplit files"
	# rm ${work}/profiletrack00${countfull}

	# for profile in $(ls ${work}/profiletrack*) ; do
	# 	awk '$1 ~ !/^>/ {if ($5 >= -60.1 && $5 <= 1 && $5 != NaN ) print $1,$2,$5 }' ${profile}  > ${profile}-points
	# 	rm ${profile}
	# done

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

		awk '$3 > max {max=$3; maxline=$0}; END {print maxline}' ${profile}-tmp > ${work}/tmp_coastpoint
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
if [ "${HASTOEXTRACT}" = true ]; then
	# Interpolation of vars to the (xi,yi) points, i.e. the baselines points
	storms=$(ls /storm*.nc)
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







# On calcul l'impact en (xi,yi) sur baseline

#                             #
# Calcul des indices d'impact #
#                             #

# nbrbouees=`awk 'END {print NR}' ./input/bouee.txt`
# nbrregimes=`awk 'END {print NR}' ./input/regime_vague.txt`
# nbbouee=$(($nbrbouees-1))
# nbregime=$(($nbrregimes-1))


# myimpactelementaire="./garbage/impact_elementaire" # numzone, numbouee, numregime, Ikji
# myimpactbouee="./garbage/impact_bouee"             # numzone, numbouee, Iji
# myimpactregime="./garbage/impact_regime"           # numzone, numregime, Iki
# myimpactglobal="./output/impact_global"            # numzone, Ii

# for i in `seq 2 $nbrprofiles`
#     do Sj=0
#        Sdbc=0
#        xi=`awk 'NR == '$i' {print $2}' $normale`
#        yi=`awk 'NR == '$i' {print $3}' $normale`
#        xo=`awk 'NR == '$i' {print $4}' $normale`
#        yo=`awk 'NR == '$i' {print $5}' $normale`
#        dio=`echo "sqrt(($xi - $xo)*($xi - $xo) + ($yi - $yo)*($yi - $yo))" |bc -l`
#        coslambda=`echo "($yo - $yi)/$dio" |bc -l`
#        sinlambda=`echo "($xi - $xo)/$dio" |bc -l`
#        numzone=$(($i-1))
#     echo $numzone
#        mytemp="./garbage/impact_temp_$numzone"
#        for j in `seq 2 $nbrbouees`
#            do Sk=0
#               xb=`awk 'NR == '$j' {print $1}' ./garbage/bouee_lambert`
#               yb=`awk 'NR == '$j' {print $2}' ./garbage/bouee_lambert`
#               dbc=`echo "sqrt(($xi - $xb)*($xi - $xb) + ($yi - $yb)*($yi - $yb))" |bc -l` # distance bouée-côte
#               cosgamma=`echo "($yi - $yb)/$dbc" |bc -l`
#               singamma=`echo "($xb - $xi)/$dbc" |bc -l`
#               numbouee=$(($j-1))
#               nregime=0
#               for k in `seq 2 $nbrregimes`
#                   do betameteo=`awk 'NR == '$k' {print $1}' ./input/regime_vague.txt`
#                      numpi=3.1415926535897932384626433832795
#                      cosbeta=`echo "c((180 - $betameteo)*$numpi/180)" |bc -l`
#                      sinbeta=`echo "s((180 - $betameteo)*$numpi/180)" |bc -l`
#                      numregime=$(($k-1))
#                      cbouee=$(($numbouee+1))
#                      pbouee=`awk 'NR == '$k' {print$'$cbouee'}' ./input/regime_vague.txt`
#                      if [ $pbouee == 0 ]
#                         then cosalpha=0
#                              cosomega=0
#                         else cosalpha=`echo "$singamma*$sinbeta + $cosgamma*$cosbeta" |bc -l`
#                              cosomega=`echo "$sinlambda*$sinbeta + $coslambda*$cosbeta" |bc -l`
#                              nregime=$(($nregime+1))
#                      fi
#                      if [ $(echo "$cosalpha < 0" |bc) = 1 ]
#                         then cosalpha=0
#                      fi
#                      if [ $(echo "$cosomega < 0" |bc) = 1 ]
#                         then cosomega=0
#                      fi
#                      Ikji=`echo "100*$cosalpha*$cosomega" |bc -l` # Impact élémentaire du régime k à la bouée j sur la zone i
#                      awk 'BEGIN { printf("%d %d %d %f\n",'$numzone','$numbouee','$numregime','$Ikji') }' >> $myimpactelementaire
#                      Sk=`echo "$Sk + $Ikji" |bc -l` # somme sur les k de Ikji
#                      Itemp=`echo "$Ikji/$dbc" |bc -l`
#                      awk 'BEGIN { printf("%d %d %d %f %f %d\n",'$numzone','$numbouee','$numregime','$Itemp','$dbc','$pbouee') }' >> $mytemp
#               done
#               Iji=`echo "$Sk/$nregime" |bc -l` # impact de tous les régime à la bouée j sur i
#               awk 'BEGIN { printf("%d %d %f\n",'$numzone','$numbouee','$Iji') }' >> $myimpactbouee
#               Sj=`echo "$Sj + $Iji/$dbc" |bc -l` # somme sur les j de Iji pondérée par 1/dbc
#               Sdbc=`echo "$Sdbc + (1/$dbc)" |bc -l` # somme sur les j des 1/dbc
#        done
#        Ii=`echo "$Sj/$Sdbc" |bc -l` # impact de tous les régime à toutes les bouées sur i
#        awk 'BEGIN { printf("%d %f\n",'$numzone','$Ii') }' >> $myimpactglobal
#        for k in `seq 1 $nbregime`
#            do awk 'BEGIN { I=0 ; D=0 ;
#                          }
#                    $3 == '$k' { I = I + $4 ; D = D + $6/$5 } 
#                    END { Iki = I/D ;
#                          printf("%d %d %f\n",'$numzone','$k',Iki) 
#                        }' $mytemp >> $myimpactregime
#        done
# done

################
## GMT PLOTS  ##
################
if [ "${GMT}" = true ]; then

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
	# gmt psxy ${work}/track.xyz $projection $envelope -S+0.4c -W1p,black -P -O -K >> ${outfile}.ps
	for profile in $(ls ${work}/profiletrack*) ; do
		tail -n +2 $profile > ${work}/profiletmp
		gmt psxy ${work}/profiletmp $projection $envelope -Sp0.1c -W0.5p,black -P -O -K >> ${outfile}.ps
	done
	log $? "psxy"

	gmt	psscale -D21/9/17.5/0.3 -C$palette -B500:"":/:"Depth(m)": -E -O  >> ${outfile}.ps
	log $? "psscale"

	gmt ps2raster -E$png_resolution -A -Tg -P ${outfile}.ps
	log $? "psraster"
	rm ${outfile}.ps
fi

rightnow
log "notice" "That's all folks !  $d"