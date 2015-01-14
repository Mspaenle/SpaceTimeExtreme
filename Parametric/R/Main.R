dyn.load(paste(getwd(),"/../resources/spt/spt/spt/","sptlib",.Platform$dynlib.ext,sep=""))
Sys.setenv(OMP_NUM_THREADS=2)

#-- Test --#
source("SptTestorOpenmp.R")