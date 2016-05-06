#!/bin/bash 

############################
## BASELINE CONFIGURATION ##
############################
export HASBASELINE=true
# mode 1 : h0 max from  => max[Tp(t)] => T_0^max => L_0^max => h0^max
# mode 2 : h0 is fixed at a reference value
# mode 3 : h0 ~ 10-20m assuming WW3 is correct close to the littoral
export BASELINEMODE=2
export H0=60
export BLOCKMEANINC=18000

######################################
## NORMAL to BASELINE CONFIGURATION ##
######################################
export HASPROFILE=true
export ORTHOLENGTH=25000
export ORTHONBPOINTS=10

##########################
## Shoreline extraction ##
##########################
export HASSHORELINE=true
export PATHTOGSHHGBIN=/Users/rchailan/Applications/gshhg/gshhg-bin-2.3.4/gshhs_f.b

###########################
## Geometry computations ##
###########################
export HASGEOMETRY=true
export POLYGONINC=3000

#########################
## Extract storms data ##
#########################
export HASEXTRACTED=false
export NEARNEIGHBORINC=1000
export NEARNEIGHBORS=10000
export NEARNEIGHBORN=4/2

export STORMSDIR=inputs/storms/NO-REF/T100
export STORM=inputs/storms/NO-REF/T100/storm-0.nc
# export STORM=inputs/storms/T75/storm-9.nc
export VARIABLES="hs_uplifted dir tp_uplifted"
# export VARIABLES="hs dir tp"


#################
## Wave impact ##
#################
export DOIMPACT=true

################
## GMT PLOTS  ##
################
export GMT=false