#! /bin/bash

source "properties.sh"
source "resources/sh/utility.sh"
rightnow
log "notice" "STARTING... $d"


for rp in 25 50 75 100 125 150; do
	export STORMSDIR=inputs/storms/NO-REF/T${rp}
	log "notice" "work on Return Period ${rp}"
	for st in 1 2 3 4 5 6 7 8 9 0; do
		export STORM=inputs/storms/NO-REF/T${rp}/storm-${st}.nc
		./impact.sh
		mv outputs/impacts-storm-$st "/Users/rchailan/Desktop/OnGoing/R5-Thesis/R/data/space-time-impacts/T0I_MODE3-NOREF/T${rp}/."
		log $? "final move" 
	done
done