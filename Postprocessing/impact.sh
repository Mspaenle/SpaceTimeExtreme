#! /bin/bash

source "properties.sh"
source "resources/sh/utility.sh"
rightnow
log "notice" "STARTING... $d"

## PARAMETERS ##
bathy="inputs/sirocco.europe.grd"

## Create Baseline from BATHY file ##
if [ "$HASBASELINE" = false ] ; then
    log "notice" "No baseline found. Create a baseline with mode $BASELINEMODE"
    case $BASELINEMODE in 
    	1 )
			echo "TODO";;
		2 ) 
			log 0 "H0 value is $H0"
			log $? "baseline created";;
		3 )
			echo "TODO";;
	esac
fi


rightnow
log "notice" "That's all folks !  $d"