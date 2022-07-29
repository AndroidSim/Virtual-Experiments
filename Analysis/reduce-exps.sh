#!/bin/bash
###
## Executes the data reduction scripts on each experiment listed as an argument.
##
## Time-stamp: <2019-09-03 10:57:32 gepr>
###
shopt -s extglob

DEFMIN=0
DEFMAX=100
PPMIN=5
PPMAX=10
#PCMIN=0
#PCMAX=8
PCMIN=0
PCMAX=4
MIDMIN=12
MIDMAX=16

SRC_DIR=$(dirname ${BASH_SOURCE[0]})
source ${SRC_DIR}/redlib.sh

function usage() {
   echo "  Usage: reduce-exps.sh <experiment 1> <experiment 2> ..."
   echo "      * Need at least 1 argument being experiment name(s)"
   exit 1
}

# arguments 1, 2, etc. are the experiment names.
nexp=$#
if [ ${nexp} -lt 1 ]; then usage; fi

echo "number of experiments = $#"

# data reduction and extraction

EXPs=$@
for exp in ${EXPs}
do
  echo "working on experiment ${exp}:"
  cleanexp=$(dirname ${exp})/$(basename ${exp}) # clean up trailing slash for file names
  
  pre_mv_whole ${exp}
  pre_mv_banded dPV ${DEFMIN} ${DEFMAX} ${cleanexp}
  pre_mv_banded dPV ${PPMIN} ${PPMAX} ${cleanexp}
  pre_mv_banded dCV ${PCMIN} ${PCMAX} ${cleanexp}
  pre_mv_banded dCV ${MIDMIN} ${MIDMAX} ${cleanexp}

  echo "moving data reduction files to analysis directory"
  mkdir -p ${cleanexp}-reduced
  mv ${cleanexp}*.csv ${cleanexp}-reduced

  # requires files from ${exp}-reduced/Hcounts-avgsd.r
  post_mv_whole ${cleanexp}

  # requires files from ${exp}-reduced
  ${SRC_DIR}/dcvdpv.r ${PPMIN} ${PPMAX} ${PCMIN} ${PCMAX} ${exp}
  if [ $? -ne 0 ]; then error "dCV/dPV dPV∈[${PPMIN},${PPMAX}), dCV∈[${PCMIN},${PCMAX}) script failed."; fi

  post_mv_banded dPV ${DEFMIN} ${DEFMAX} ${cleanexp}
  post_mv_banded dPV ${PPMIN} ${PPMAX} ${cleanexp}
  post_mv_banded dCV ${PCMIN} ${PCMAX} ${cleanexp}
  post_mv_banded dCV ${MIDMIN} ${MIDMAX} ${exp}
   
done

exit 0
