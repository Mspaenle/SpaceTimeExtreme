#!/bin/bash 

############################
## BASELINE CONFIGURATION ##
############################
export HASBASELINE=false
# mode 1 : h0 max from  => max[Tp(t)] => T_0^max => L_0^max => h0^max
# mode 2 : h0 is fixed at a reference value
# mode 3 : h0 ~ 10-20m assuming WW3 is correct close to the littoral
export BASELINEMODE=2
export H0=60
# In degrees for the moment
export EQUIDIST=0.01 

######################################
## NORMAL to BASELINE CONFIGURATION ##
######################################
export HASPROFILE=true

# following three dist needs to be in the same unit (see -C grdtrack)
#export PROFILELENGTH=70k
#export EQUIPROFILEDS=0.1k
#export EQUIPROFILESPACE=2k

export ORTHOLENGTH=0.5
export ORTHONBPOINTS=30

###########################
## Geometry computations ##
###########################
export HASGEOMETRY=true

#########################
## Extract storms data ##
#########################
export HASTOEXTRACT=false
export NEARNEIGHBORINC=0.01
export NEARNEIGHBORS=5k
export NEARNEIGHBORN=4/2

export STORMSDIR=inputs/storms
export VARIABLES="hs_uplifted dir fp"

################
## GMT PLOTS  ##
################
export GMT=true