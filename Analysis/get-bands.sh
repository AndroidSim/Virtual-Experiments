#!/bin/bash
###
## Executes the data reduction scripts for a set of bands 
## on each experiment listed as an argument.
## This is specifically used to obtain the spatial bands for a 3D plot.
##
## Time-stamp: 
###
shopt -s extglob

## first argument is max distance for dCV or dPV obtained from the intersection point (xpt) 
## this is used to divide the space into 2 sections: from CV to xpt and from PV to xpt
## the second argument is the step size or band width
## e.g. set of bands in intervals [0,STEP) until [MAX-1,MAX)

# arguments 3, 4, etc. are the experiment names.

if [ $# -lt 4 ]
then
   echo "  Usage: get-bands.sh min max step <experiment 1> <experiment 2> ..."
   echo "      * Need at least 1 argument being experiment name(s)"
   exit 1
else
   nexp=$(($#-3))
   echo "number of experiments = ${nexp}"
fi

MIN=$1
MAX=$2
STEP=$3
shift 3

SRC_DIR=$(dirname ${BASH_SOURCE[0]})

# get the banded measurements for each experiment

EXPs=$@
for exp in ${EXPs}
do
	START=${MIN}
	STOP=$((${MIN}+${STEP}))
	while [ ${STOP} -le ${MAX} ]
	do
		echo "reducing cell state event data ..."
		${SRC_DIR}/reduce-event-data-inband.r dCV ${START} ${STOP} ${exp}
		echo "metabolism measurements..."
		${SRC_DIR}/inextra-inband.r dCV ${START} ${STOP} ${exp}
		echo "dataperH-inband.r ${DEFMIN} ${DEFMAX} ${exp}"
		${SRC_DIR}/dataperH-inband.r dCV ${START} ${STOP} ${exp}
		echo "calculating exposure with ma = 51 #vHPC = 1000 "
		#${SRC_DIR}/calc-exposure.r dPV ${START} ${STOP} ${exp}
		${SRC_DIR}/calc-exposure.r dCV ${START} ${STOP} ${exp}
		START=$((${START}+${STEP}))
		STOP=$((${STOP}+${STEP}))
	done
	echo "moving band files to -reduced directory"
	mkdir -p ${exp}-reduced
	mv ${exp}*.csv ${exp}-reduced
done

exit 0
