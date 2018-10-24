#!/bin/bash
###
## Executes the data reduction scripts on each experiment listed as an argument.
##
## Time-stamp: <2018-10-16 10:28:34 gepr>
###
shopt -s extglob

DEFMIN=0
DEFMAX=100
PPMIN=5
PPMAX=10
PCMIN=0
PCMAX=8

# arguments 1, 2, etc. are the experiment names.

if [ $# -lt 1 ]
then
   echo "  Usage: reduce-exps.sh <experiment 1> <experiment 2> ..."
   echo "      * Need at least 1 argument being experiment name(s)"
   exit 1
else
   nexp=$#
   echo "number of experiments = $#"
fi

SRC_DIR=$(dirname ${BASH_SOURCE[0]})

# data reduction and extraction

EXPs=$@
for exp in ${EXPs}
do
   echo "working on experiment ${exp}:"
   echo "solute measurements..."
   ${SRC_DIR}/body-avg.r ${exp}
   ${SRC_DIR}/extra-avg.r ${exp}
   ${SRC_DIR}/fextract.r ${exp}
   echo "reducing cell state event data ..."
   ${SRC_DIR}/reduce-event-data-inband.r ${DEFMIN} ${DEFMAX} ${exp}
   ${SRC_DIR}/reduce-event-data-inband.r ${PPMIN} ${PPMAX} ${exp}
   ${SRC_DIR}/reduce-event-data-inband.r ${PCMIN} ${PCMAX} ${exp}
   echo "enzyme measurements..."
   ${SRC_DIR}/eg-dPV.r ${exp}
   ${SRC_DIR}/eg-inband.r ${DEFMIN} ${DEFMAX} ${exp}_enzymes.csv
   ${SRC_DIR}/eg-inband.r ${PPMIN} ${PPMAX} ${exp}_enzymes.csv
   ${SRC_DIR}/ei-groups.r ${exp}
   echo "metabolism measurements..."
   ${SRC_DIR}/inextra-inband.r ${DEFMIN} ${DEFMAX} ${exp}
   ${SRC_DIR}/inextra-inband.r ${PPMIN} ${PPMAX} ${exp}
   ${SRC_DIR}/inextra-inband.r ${PCMIN} ${PCMAX} ${exp}
   ${SRC_DIR}/rxnfield.r ${exp}
   echo "averaging #Hepatocytes ∈[${DEFMIN},${DEFMAX})"
   ${SRC_DIR}/hcounts-inband.r dPV ${DEFMIN} ${DEFMAX} ${exp} > ${exp}_hcounts-dCV∈\[${DEFMIN}\,${DEFMAX}\).csv
   ${SRC_DIR}/hcounts-inband.r dCV ${PCMIN} ${PCMAX} ${exp} > ${exp}_hcounts-dCV∈\[${PCMIN}\,${PCMAX}\).csv
   ${SRC_DIR}/hcounts-inband.r dPV ${PPMIN} ${PPMAX} ${exp} > ${exp}_hcounts-dPV∈\[${PPMIN}\,${PPMAX}\).csv
   echo "dataperH-inband.r ${DEFMIN} ${DEFMAX} ${exp}"
   ${SRC_DIR}/dataperH-inband.r ${DEFMIN} ${DEFMAX} ${exp}
   ${SRC_DIR}/dataperH-inband.r ${PPMIN} ${PPMAX} ${exp}
   ${SRC_DIR}/dataperH-inband.r ${PCMIN} ${PCMAX} ${exp}
   echo "calc-exposure.r 51 1000 ${exp}_entries-dPV-avg-pHPC-pMC∈\[${DEFMIN}\,${DEFMAX}\).csv"
   ${SRC_DIR}/calc-exposure.r 51 1000 ${exp}_entries-avg-pHPC-pMC-dPV∈\[${DEFMIN}\,${DEFMAX}\).csv
   ${SRC_DIR}/calc-exposure.r 51 1000 ${exp}_entries-avg-pHPC-pMC-dPV∈\[${PPMIN}\,${PPMAX}\).csv
   ${SRC_DIR}/calc-exposure.r 51 1000 ${exp}_entries-avg-pHPC-pMC-dCV∈\[${PCMIN}\,${PCMAX}\).csv
   echo "Hsums-Dist-MC.r ${DEFMIN} ${DEFMAX} ${exp}"
   ${SRC_DIR}/Hsums-Dist-MC.r ${DEFMIN} ${DEFMAX} ${exp}
   ${SRC_DIR}/Hsums-Dist-MC.r ${PPMIN} ${PPMAX} ${exp}
   ${SRC_DIR}/Hsums-Dist-MC.r ${PCMIN} ${PCMAX} ${exp}
   echo "Hcounts-avgsd.r ${DEFMIN} ${DEFMAX} ${exp}"
   ${SRC_DIR}/Hcounts-avgsd.r ${DEFMIN} ${DEFMAX} ${exp}
   ${SRC_DIR}/Hcounts-avgsd.r ${PPMIN} ${PPMAX} ${exp}
   ${SRC_DIR}/Hcounts-avgsd.r ${PCMIN} ${PCMAX} ${exp}
   echo "moving data reduction files to analysis directory"
   mkdir -p ${exp}-reduced
   mv ${exp}*.csv ${exp}-reduced

   ## the following derived data assumes inputs are in ${exp}-reduced
   echo "dCV/dPV derived measures."
   ${SRC_DIR}/dcvdpv.r ${exp}
done

exit 0
