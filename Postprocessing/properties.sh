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
export HASPROFILE=false
export ORTHOLENGTH=25000
export ORTHONBPOINTS=80

##########################
## Shoreline extraction ##
##########################
export HASSHORELINE=true
export PATHTOGSHHGBIN=/Users/rchailan/Applications/gshhg/gshhg-bin-2.3.4/gshhs_f.b

###########################
## Geometry computations ##
###########################
export HASGEOMETRY=true

#########################
## Extract storms data ##
#########################
export HASEXTRACTED=true
export NEARNEIGHBORINC=0.01
export NEARNEIGHBORS=5k
export NEARNEIGHBORN=4/2

export STORMSDIR=inputs/storms
export VARIABLES="hs_uplifted dir fp"

################
## GMT PLOTS  ##
################
export GMT=false