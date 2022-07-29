#!/bin/bash
###
## Executes the single data plotting script plot-cols.r
## for a single experiment each from a set of experiments given as arguments.
##
## Time-stamp: <2020-07-13 aks>
###
shopt -s extglob

function error() { echo $'\n'"!!!Error: $@"$'\n'; exit -1; }
function usage() {
	echo "  Usage: plot-exps.sh <experiment 1> <experiment 2> ..."
	echo "    	The data is plotted in graphics directory."
   exit 1
}

# arguments 1, 2, etc. are the experiment names.
nexp=$#
if [ ${nexp} -lt 1 ]; then usage; fi

echo "number of experiments = $#"
exps=$*

SRC_DIR=$(dirname ${BASH_SOURCE[0]}) # this script's location
GH_DIR="${HOME}/Research/BioSystemsGroup-Github/scripts" # github scripts location
#GH_DIR="${HOME}/local/scripts" # github scripts location


function warn() { echo "Warning: $@"; }

move() {
  # 1st arg = data name, 2nd arg = experiment name 
  # move plot files to target directory	
  tgt_dir="g-$1"
  if ! test -e ${tgt_dir}; then mkdir ${tgt_dir}; fi
  if ! test -z "$(ls -A graphics/)" ; then mv graphics/* ${tgt_dir}/ ; fi
}

plot() {
  type=$1; shift 1
  if [[ $1 == "" ]]; then return; fi
  if [[ ${type} == "raw" ]]
  then
    ${GH_DIR}/plot-cols.r raw $@
    if [ $? -ne 0 ]; then error "${GH_DIR}/plot-cols.r raw $@ failed."; fi
  else
    ${GH_DIR}/plot-cols.r ${type} $@
    if [ $? -ne 0 ]; then error "${GH_DIR}/plot-cols.r ${type} $@ failed."; fi
  fi
}

function getBands() {
  dir=$1; exp=$2
  regexstr=".*${dir}∈\[[[:digit:]]*,[[:digit:]]*).*.csv"
  directory="$exp-reduced"
  # "grep -v" needed to remove the ratio files
  files=$(find -L ${directory} -regextype posix-basic -regex ${regexstr} | grep -v "\-to\-")
  bandNdx=1
  for f in $files; do
    bands="$bands "$(echo $f | sed 's/.*∈//' | sed 's/).*/)/')
    bandNdx=$(( bandNdx + 1))
  done
  echo $( set -f; printf "%s\n" $bands | sort -u | paste -sd" " )
}

function plotBandsInDir() {
  dir=$1; shift; exp=$1; shift; bands=$@;
  
  for band in ${bands}; do
    echo "Plotting Necrosed ${dir}∈${band}"
    file="${exp}-reduced/${exp}_necrotic-${dir}∈${band}.csv"
    if test -e ${file}; then plot data ${file}; move necrotic 
    else warn "No necrotic-${dir}∈${band} files to plot."; fi
    
    echo "Plotting Triggered ${dir}∈${band}"
    file="${exp}-reduced/${exp}_nectrig-${dir}∈${band}.csv"
    if test -e ${file}; then plot data ${file}; move nectrig 
    else warn "No nectrig-${dir}∈${band} files to plot."; fi
    
    echo "Plotting stressed ${dir}∈${band}"
    file="${exp}-reduced/${exp}_stressed-${dir}∈${band}.csv"
    if test -e ${file}; then plot data ${file}; move stressed 
    else warn "No stressed-${dir}∈${band} files to plot."; fi
    
    echo "Plotting intra ${dir}∈${band}"
    file="${exp}-reduced/${exp}_mobileObject-${dir}∈${band}.csv"
    if test -e ${file}; then plot data ${file}; move intra 
    else warn "No mobileObject-${dir}∈${band} files to plot."; fi
    
    echo "Plotting celladj ${dir}∈${band}"
    file="${exp}-reduced/${exp}_celladj-${dir}∈${band}.csv"
    if test -e ${file}; then plot data ${file}; move celladj 
    else warn "No celladj-${dir}∈${band} files to plot."; fi
    
    echo "Plotting enzymes ${dir}∈${band}"
    file="${exp}-reduced/${exp}_enzymes-${dir}∈${band}.csv"
    if test -e ${file}; then plot data ${file}; move enzymes 
    else warn "No enzymes-∈${band} files to plot."; fi
    
    echo "Plotting MITs ${dir}∈${band}" 
    file="${exp}-reduced/${exp}_entries-${dir}∈${band}.csv"
    if test -e ${file}; then plot rawlines ${file}; move mits 
    else warn "No entries-${dir}∈${band} files to plot."; fi
      
    file="${exp}-reduced/${exp}_exits-${dir}∈${band}.csv"
    if test -e ${file}; then plot rawlines ${file}; move mits  
    else warn "No exits-${dir}∈${band} files to plot."; fi
      
    file="${exp}-reduced/${exp}_rejects-${dir}∈${band}.csv"
    if test -e ${file}; then plot rawlines ${file}; move mits
    else warn "No rejects-${dir}∈${band} files to plot."; fi
      
    file="${exp}-reduced/${exp}_traps-${dir}∈${band}.csv"
    if test -e ${file}; then plot rawlines ${file}; move mits
    else warn "No traps-${dir}∈${band} files to plot."; fi
    
    echo "Plotting exposure ${dir}∈${band}"
    file="${exp}-reduced/${exp}_exposure-entries-${dir}∈${band}.csv"
    if test -e ${file}; then plot data ${file}; move exposure 
    else warn "No exposure-entries-${dir}∈${band} files to plot."; fi
    
    echo "Plotting clearance ${dir}∈${band}"
    file="${exp}-reduced/${exp}_clearance-[entries-exits]-${dir}∈${band}.csv"
    if test -e ${file}; then plot data ${file}; move clearance 
    else warn "No clearance-[entries-exits]-${dir}∈${band} files to plot."; fi
    
  done 
}

## for each experiment, plot the data files, then move files to directory
for exp in ${exps}
do
  ## non-banded data
  echo "Plotting body file"
  file="${exp}-reduced/${exp}_body-avg.csv"
  if test -e ${file}; then plot rawlines ${file}; move body 
  else warn "No body files to plot."; fi
  
  echo "Plotting outfract file"
  file="${exp}-reduced/${exp}_outFract.csv"
  if test -e ${file}; then plot data ${file}; move outFract 
  else warn "No OutFract files to plot."; fi
  
  echo "Plotting extratio file"
  file="${exp}-reduced/${exp}_extRatio.csv"
  if test -e ${file}; then plot data ${file}; move extRatio 
  else warn "No ExtRatio files to plot.";  fi
  
  echo "Plotting extra file"
  file="${exp}-reduced/${exp}_extra-avg.csv"
  if test -e ${file}; then plot data ${file}; move extra 
  else warn "No extraCellular files to plot."; fi
  
  ## banded data
  # find directions and bands to plot
  echo "Plotting dPV bands"
  dPVbands=$(getBands dPV ${exp})
  plotBandsInDir "dPV" ${exp} ${dPVbands} 

  echo "Plotting dCV bands"
  dCVbands=$(getBands dCV ${exp})
  plotBandsInDir "dCV" ${exp} ${dCVbands} 
  
  # move g- directories whole experiment directory, i.e. "exp-plots"
  exp_dir="${exp}-plots"
  if test -e ${exp_dir}; 
  then rm -rf ${exp_dir}; mkdir ${exp_dir}; mv g-* ${exp_dir};
  else mkdir ${exp_dir}; mv g-* ${exp_dir};
  fi
done

exit 0
