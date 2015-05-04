#!/bin/bash 

# BASELINE CONFIGURATION #
export HASBASELINE=false
# mode 1 : h0 max from  => max[Tp(t)] => T_0^max => L_0^max => h0^max
# mode 2 : h0 is fixed at a reference value
# mode 3 : h0 ~ 10-20m assuming WW3 is correct close to the littoral
export BASELINEMODE=2
export H0=60