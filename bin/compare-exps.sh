#!/bin/bash
###
## Executes the comparison scripts cmp-by-col.r and cmp-movavg.r
## that compare the output between a set of experiments given as arguments.
##
## Time-stamp: <2018-09-05 09:36:14 gepr>
###
shopt -s extglob

DEFBAND="[0,100)"
PPBAND="[5,10)"
PCBAND="[0,8)"

move() {
  tgt_dir="g-$1"
  if ! test -e ${tgt_dir}; then mkdir ${tgt_dir}; fi
  if ! test -z "$(ls -A graphics/)" ; then mv graphics/* ${tgt_dir}/ ; fi
}

# arguments 1, 2, etc. are the experiment names.

if [ $# -lt 2 ]
then
	echo "  Usage: compare-exps.sh <experiment 1> <experiment 2> ..."
	echo "		Need at least 2 experiments to compare"
	echo "    	The data is plotted in graphics directory."
	exit 1
else
	nexp=$#
	echo "number of experiments = $#"
fi

SRC_DIR=$(dirname ${BASH_SOURCE[0]}) # this script's location
#GH_DIR="${HOME}/Research/BioSystemsGroup-Github/scripts" # github scripts location
GH_DIR="${HOME}/local/scripts" # github scripts location
plot() {
  type=$1; shift 1
  if [[ $1 == "" ]]; then echo "No files to plot."; return; fi
  if [[ ${type} == "raw" ]]; then ${SRC_DIR}/cmp-by-col.r $@
  else ${GH_DIR}/cmp-movavg.r $@
  fi
}

## non-banded data
n=1
for exp in $@
do
  file="${exp}-reduced/${exp}_body-avg.csv"
  if test -e ${file}; then body_files[${n}]=${file}; fi
  file="${exp}-reduced/${exp}_outFract.csv"
  if test -e ${file}; then outFract_files[${n}]=${file}; fi
  file="${exp}-reduced/${exp}_extRatio.csv"
  if test -e ${file}; then extRatio_files[${n}]=${file}; fi
  file="${exp}-reduced/${exp}_extra.csv"
  if test -e ${file}; then extra_files[${n}]=${file}; fi
  (( n++ ))
done
plot raw ${body_files[@]}
move body
plot avg ${outFract_files[@]}
move outFract
plot avg ${extRatio_files[@]}
move extRatio
plot raw ${extra_files[@]}
move extra

## dPV data
for band in $DEFBAND $PPBAND; do
  n=1
  for exp in $@; do
    file="${exp}-reduced/${exp}_necrotic-dPV∈${band}.csv"
    if test -e ${file}; then necroticdPV_files[${n}]=${file}; fi
    file="${exp}-reduced/${exp}_nectrig-dPV∈${band}.csv"
    if test -e ${file}; then nectrigdPV_files[${n}]=${file}; fi
    file="${exp}-reduced/${exp}_hsolute-dPV∈${band}.csv"
    if test -e ${file}; then hsoldPV_files[${n}]=${file}; fi
    file="${exp}-reduced/${exp}_celladj-dPV∈${band}.csv"
    if test -e ${file}; then celladjdPV_files[${n}]=${file}; fi
    file="${exp}-reduced/${exp}_enzymes∈${band}.csv"
    if test -e ${file}; then eg_files[${n}]=${file}; fi
    file="${exp}-reduced/${exp}_entries-dPV∈${band}.csv"
    if test -e ${file}; then entriesdPV_files[${n}]=${file}; fi
    file="${exp}-reduced/${exp}_exits-dPV∈${band}.csv"
    if test -e ${file}; then exitsdPV_files[${n}]=${file}; fi
    file="${exp}-reduced/${exp}_rejects-dPV∈${band}.csv"
    if test -e ${file}; then rejectsdPV_files[${n}]=${file}; fi
    file="${exp}-reduced/${exp}_traps-dPV∈${band}.csv"
    if test -e ${file}; then trapsdPV_files[${n}]=${file}; fi
    file="${exp}-reduced/${exp}_entries-avg-pHPC-pMC-dPV∈${band}.csv"
    if test -e ${file}; then exposuredPV_files[${n}]=${file}; fi
    (( n++ ))
  done

  echo "Plotting Necrosed dPV∈${band}"
  plot avg ${necroticdPV_files[@]}
  move necrotic
  echo "Plotting Triggered dPV∈${band}"
  plot avg ${nectrigdPV_files[@]}
  move nectrig
  echo "Plotting intra dPV∈${band}"
  plot avg ${hsoldPV_files[@]}
  move intra
  echo "Plotting celladj dPV∈${band}"
  plot avg ${celladjdPV_files[@]}
  move celladj
  echo "Plotting enzymes dPV∈${band}"
  plot avg ${eg_files[@]}
  move enzymes
  echo "Plotting MITs dPV∈${band}"
  plot raw ${entriesdPV_files[@]}
  plot raw ${exitsdPV_files[@]}
  plot raw ${rejectsdPV_files[@]}
  plot raw ${trapsdPV_files[@]}
  move mits
  echo "Plotting exposure dPV∈${band}"
  plot raw ${exposuredPV_files[@]}
  move exposure
done

## dCV data
for band in $DEFBAND $PCBAND; do
  n=1
  for exp in $@; do
    file="${exp}-reduced/${exp}_necrotic-dCV∈${band}.csv"
    if test -e ${file}; then necroticdCV_files[${n}]=${file}; fi
    file="${exp}-reduced/${exp}_nectrig-dCV∈${band}.csv"
    if test -e ${file}; then nectrigdCV_files[${n}]=${file}; fi
    file="${exp}-reduced/${exp}_hsolute-dCV∈${band}.csv"
    if test -e ${file}; then hsoldCV_files[${n}]=${file}; fi
    file="${exp}-reduced/${exp}_celladj-dCV∈${band}.csv"
    if test -e ${file}; then celladjdCV_files[${n}]=${file}; fi
    file="${exp}-reduced/${exp}_entries-dCV∈${band}.csv"
    if test -e ${file}; then entriesdCV_files[${n}]=${file}; fi
    file="${exp}-reduced/${exp}_exits-dCV∈${band}.csv"
    if test -e ${file}; then exitsdCV_files[${n}]=${file}; fi
    file="${exp}-reduced/${exp}_rejects-dCV∈${band}.csv"
    if test -e ${file}; then rejectsdCV_files[${n}]=${file}; fi
    file="${exp}-reduced/${exp}_traps-dCV∈${band}.csv"
    if test -e ${file}; then trapsdCV_files[${n}]=${file}; fi
    file="${exp}-reduced/${exp}_entries-avg-pHPC-pMC-dCV∈${band}.csv"
    if test -e ${file}; then exposuredCV_files[${n}]=${file}; fi
    (( n++ ))
  done

  echo "Plotting Necrosed dCV∈${band}"
  plot avg ${necroticdCV_files[@]}
  move necrotic
  echo "Plotting Triggered dCV∈${band}"
  plot avg ${nectrigdCV_files[@]}
  move nectrig
  echo "Plotting intra dCV∈${band}"
  plot avg ${hsoldCV_files[@]}
  move intra
  echo "Plotting celladj dCV∈${band}"
  plot avg ${celladjdCV_files[@]}
  move celladj
  echo "Plotting MITs dCV∈${band}"
  plot raw ${entriesdCV_files[@]}
  plot raw ${exitsdCV_files[@]}
  plot raw ${rejectsdCV_files[@]}
  plot raw ${trapsdCV_files[@]}
  move mits
  echo "Plotting exposure dCV∈${band}"
  plot raw ${exposuredCV_files[@]}
  move exposure
done

exit 0
