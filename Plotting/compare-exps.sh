#!/bin/bash
###
## Executes the comparison scripts cmp-by-col.r and cmp-movavg.r
## that compare the output between a set of experiments given as arguments.
##
## Time-stamp: <2019-11-29 20:37:38 gepr>
###
shopt -s extglob

SRC_DIR=$(dirname ${BASH_SOURCE[0]}) # this script's location
GH_DIR="${HOME}/Research/BioSystemsGroup-Github/scripts" # github scripts location
#GH_DIR="${HOME}/local/scripts" # github scripts location

function usage() {
	echo "  Usage: compare-exps.sh <experiment 1> <experiment 2> ..."
	echo "		Need at least 2 experiments to compare"
	echo "    	The data is plotted in graphics directory."
   exit 1
}

# arguments 1, 2, etc. are the experiment names.
nexp=$#
if [ ${nexp} -lt 2 ]; then usage; fi

echo "number of experiments = $#"
exps=("$@")

function error() { echo $'\n'"!!!Error: $@"$'\n'; exit -1; }
function warn() { echo "Warning: $@"; }
## check if experiment directories are empty
## arguments are experiment output directory
array=${exps[@]}
#array=("${exps[@]}")
index=0
for a in $array
do
	ar="$a-reduced"
	if [ ! -d ${ar} ];
	then
		warn "the experiment directory $ar does not exist; therefore, removing it";
		#array=( ${array[@]/$delete/} )
		exps=( ${exps[@]/%"$a"/} )
		#unset array[$index]?
	fi
	((index++))
done
array=${exps[@]}
for a in $array
do
	ar="$a-reduced"
	if [ ! "$(ls -A $ar)" ];
	then
		warn "the experiment directory $ar is empty; therefore, removing it";
		#array=( ${array[@]/$delete/} )
		exps=( ${exps[@]/%"$a"/} )
	fi
done
exps=${exps[@]}
nexp=${#exps[@]}
if [ ${nexp} -lt 2 ]; then usage; fi

#if [ ${#array[@]} -eq 0 ]; then
#		echo "array is empty"
#fi
	
function move() {
  tgt_dir="g-$1"
  if ! test -e ${tgt_dir}; then mkdir ${tgt_dir}; fi
  if ! test -z "$(ls -A graphics/)" ; then mv graphics/* ${tgt_dir}/ ; fi
}

function plot() {
  type=$1; shift 1
  if [[ $1 == "" ]]; then return; fi
  if [[ ${type} == "raw" ]]
  then
    ${GH_DIR}/cmp-plot.r raw $@
    if [ $? -ne 0 ]; then error "${GH_DIR}/cmp-plot.r raw $@ failed."; fi
  else
    ${GH_DIR}/cmp-plot.r ${type} $@
    if [ $? -ne 0 ]; then error "${GH_DIR}/cmp-plot.r ${type} $@ failed."; fi
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

function commonBands() {
	array1=$1
	array2=$2

	for item1 in "${array1[@]}"; do
		for item2 in "${array2[@]}"; do
			if [[ $item1 == "$item2" ]]; then
				intersections+=( "$item1" )
			fi
		done
	done

	echo  ${intersections[@]}
}

function cmpBandsInDir() {
  dir=$1; shift; band=$1; shift; exps=$@;
  
  n=1
  for exp in ${exps}; do
    # gather banded necrotic files
    file="${exp}-reduced/${exp}_necrotic-${dir}∈${band}.csv"
    if test -e ${file}; then necrotic_files[${n}]=${file}
    else warn "No necrotic-${dir}∈${band} files to plot."; fi
    
    # gather banded nectrig files
    file="${exp}-reduced/${exp}_nectrig-${dir}∈${band}.csv"
    if test -e ${file}; then nectrig_files[${n}]=${file}
    else warn "No nectrig-${dir}∈${band} files to plot."; fi
    
    # gather banded stressed files
    file="${exp}-reduced/${exp}_stressed-${dir}∈${band}.csv"
    if test -e ${file}; then stressed_files[${n}]=${file}
    else warn "No stressed-${dir}∈${band} files to plot."; fi
    
    # gather banded mobile object files
    file="${exp}-reduced/${exp}_mobileObject-${dir}∈${band}.csv"
    if test -e ${file}; then mobileObject_files[${n}]=${file}
    else warn "No mobileObject-${dir}∈${band} files to plot."; fi
    
    # gather banded mobile object per Hepatocyte files
    file="${exp}-reduced/${exp}_mobileObject-avg-pHPC-pMC-${dir}∈${band}.csv"
    if test -e ${file}; then mobileObject_pH_files[${n}]=${file}
    else warn "No mobileObject-avg-pHPC-pMC-${dir}∈${band} files to plot."; fi
    
    # gather banded bound mobile object files
    file="${exp}-reduced/${exp}_boundMObject-${dir}∈${band}.csv"
    if test -e ${file}; then boundMObject_files[${n}]=${file}
    else warn "No boundMObject-${dir}∈${band} files to plot."; fi
    
    # gather banded bound mobile object per Hepatocyte files
    file="${exp}-reduced/${exp}_boundMObject-avg-pHPC-pMC-${dir}∈${band}.csv"
    if test -e ${file}; then boundMObject_pH_files[${n}]=${file}
    else warn "No boundMObject-avg-pHPC-pMC-${dir}∈${band} files to plot."; fi
    
    # gather banded celladj files
    file="${exp}-reduced/${exp}_celladj-${dir}∈${band}.csv"
    if test -e ${file}; then celladj_files[${n}]=${file}
    else warn "No celladj-${dir}∈${band} files to plot."; fi
    
    # gather banded celladj per Hepatocyte files
    file="${exp}-reduced/${exp}_celladj-avg-pHPC-pMC-${dir}∈${band}.csv"
    if test -e ${file}; then celladj_pH_files[${n}]=${file}
    else warn "No celladj-avg-pHPC-pMC-${dir}∈${band} files to plot."; fi
    
    # gather banded enzyme files
    file="${exp}-reduced/${exp}_enzymes-${dir}∈${band}.csv"
    if test -e ${file}; then eg_files[${n}]=${file}
    else warn "No enzymes-∈${band} files to plot."; fi
    
    # gather banded enzymes per Hepatocyte files
    file="${exp}-reduced/${exp}_enzymes-avg-pHPC-pMC-${dir}∈${band}.csv"
    if test -e ${file}; then eg_pH_files[${n}]=${file}
    else warn "No enzymes-avg-pHPC-pMC-${dir}∈${band} files to plot."; fi
    
    # gather banded entries files
    file="${exp}-reduced/${exp}_entries-${dir}∈${band}.csv"
    if test -e ${file}; then entries_files[${n}]=${file}
    else warn "No entries-${dir}∈${band} files to plot."; fi
    
    # gather banded entries per Hepatocyte files
    file="${exp}-reduced/${exp}_entries-avg-pHPC-pMC-${dir}∈${band}.csv"
    if test -e ${file}; then entries_pH_files[${n}]=${file}
    else warn "No entries-avg-pHPC-pMC-${dir}∈${band} files to plot."; fi
    
    # gather banded exits files
    file="${exp}-reduced/${exp}_exits-${dir}∈${band}.csv"
    if test -e ${file}; then exits_files[${n}]=${file}
    else warn "No exits-${dir}∈${band} files to plot."; fi
    
    # gather banded exits per Hepatocyte files
    file="${exp}-reduced/${exp}_exits-avg-pHPC-pMC-${dir}∈${band}.csv"
    if test -e ${file}; then exits_pH_files[${n}]=${file}
    else warn "No exits-avg-pHPC-pMC-${dir}∈${band} files to plot."; fi
    
    # gather banded rejects files
    file="${exp}-reduced/${exp}_rejects-${dir}∈${band}.csv"
    if test -e ${file}; then rejects_files[${n}]=${file}
    else warn "No rejects-${dir}∈${band} files to plot."; fi
    
    # gather banded rejects per Hepatocyte files
    file="${exp}-reduced/${exp}_rejects-avg-pHPC-pMC-${dir}∈${band}.csv"
    if test -e ${file}; then rejects_pH_files[${n}]=${file}
    else warn "No rejects-avg-pHPC-pMC-${dir}∈${band} files to plot."; fi
    
    # gather banded traps files
    file="${exp}-reduced/${exp}_traps-${dir}∈${band}.csv"
    if test -e ${file}; then traps_files[${n}]=${file}
    else warn "No traps-${dir}∈${band} files to plot."; fi
    
    # gather banded traps per Hepatocyte files
    file="${exp}-reduced/${exp}_traps-avg-pHPC-pMC-${dir}∈${band}.csv"
    if test -e ${file}; then traps_pH_files[${n}]=${file}
    else warn "No traps-avg-pHPC-pMC-${dir}∈${band} files to plot."; fi
    
    # gather banded exposure files
    file="${exp}-reduced/${exp}_exposure-entries-${dir}∈${band}.csv"
    if test -e ${file}; then exposure_files[${n}]=${file}
    else warn "No exposure-entries-${dir}∈${band} files to plot."; fi
    
    # gather banded exposure per Hepatocyte files
    file="${exp}-reduced/${exp}_exposure-entries-perH-${dir}∈${band}.csv"
    if test -e ${file}; then exposure_pH_files[${n}]=${file}
    else warn "No exposure-entries-perH-${dir}∈${band} files to plot."; fi
    
    # gather banded clearance files
    file="${exp}-reduced/${exp}_clearance-[entries-exits]-${dir}∈${band}.csv"
    if test -e ${file}; then clearance_files[${n}]=${file}
    else warn "No clearance-[entries-exits]-${dir}∈${band} files to plot."; fi
    
    # gather banded clearance per Hepatocyte files
    file="${exp}-reduced/${exp}_clearance-[entries-exits]-perH-${dir}∈${band}.csv"
    if test -e ${file}; then clearance_pH_files[${n}]=${file}
    else warn "No clearance-[entries-exits]-perH-${dir}∈${band} files to plot."; fi
    
    # gather banded total amount files
    file="${exp}-reduced/${exp}_totalamt-${dir}∈${band}.csv"
    if test -e ${file}; then totalamt_files[${n}]=${file}
    else warn "No totalamt-${dir}∈${band} files to plot."; fi
    
    # gather banded total amount per Hepatocyte files
    file="${exp}-reduced/${exp}_totalamt-perH-${dir}∈${band}.csv"
    if test -e ${file}; then totalamt_pH_files[${n}]=${file}
    else warn "No totalamt-perH-${dir}∈${band} files to plot."; fi
    
    (( n++ ))
  done

    echo "Plotting Necrosed ${dir}∈${band}"
    plot data ${necrotic_files[@]}
    move necrotic
    echo "Plotting Triggered ${dir}∈${band}"
    plot data ${nectrig_files[@]}
    move nectrig
    echo "Plotting stressed ${dir}∈${band}"
    plot data ${stressed_files[@]}
    move stressed
    echo "Plotting intra ${dir}∈${band}"
    plot data ${mobileObject_files[@]}
    plot data ${mobileObject_pH_files[@]}
    move intra
    echo "Plotting bound intra ${dir}∈${band}"
    plot data ${boundMObject_files[@]}
    plot data ${boundMObject_pH_files[@]}
    move intra-bound
    echo "Plotting celladj ${dir}∈${band}"
    plot data ${celladj_files[@]}
    plot data ${celladj_pH_files[@]}
    move celladj
    echo "Plotting enzymes ${dir}∈${band}"
    plot data ${eg_files[@]}
    plot data ${eg_pH_files[@]}
    move enzymes
    echo "Plotting MITs ${dir}∈${band}"
    plot raw ${entries_files[@]}
    plot raw ${exits_files[@]}
    plot raw ${rejects_files[@]}
    plot raw ${traps_files[@]}
    plot raw ${entries_pH_files[@]}
    plot raw ${exits_pH_files[@]}
    plot raw ${rejects_pH_files[@]}
    plot raw ${traps_pH_files[@]}
    move mits
    echo "Plotting exposure ${dir}∈${band}"
    plot data ${exposure_files[@]}
    plot data ${exposure_pH_files[@]}
    move exposure
    echo "Plotting clearance ${dir}∈${band}"
    plot nodata ${clearance_files[@]}
    plot nodata ${clearance_pH_files[@]}
    move clearance
    echo "Plotting total amounts ${dir}∈${band}"
    plot data ${totalamt_files[@]}
    plot data ${totalamt_pH_files[@]}
    move totalamt

# cleanup
    unset necrotic_files
    unset nectrig_files
    unset stressed_files
    unset mobileObject_files
    unset mobileObject_pH_files
    unset boundMObject_files
    unset boundMObject_pH_files
    unset celladj_files
    unset celladj_pH_files
    unset eg_files
    unset eg_pH_files
    unset entries_files
    unset entries_pH_files
    unset exits_files
    unset exits_pH_files
    unset rejects_files
    unset rejects_pH_files
    unset traps_files
    unset traps_pH_files
    unset exposure_files
    unset exposure_pH_files
    unset clearance_files
    unset clearance_pH_files
    unset totalamt_files
    unset totalamt_pH_files
}


# non-banded data
n=1
for exp in ${exps}
do
  file="${exp}-reduced/${exp}_body-avg.csv"
  if test -e ${file}; then body_files[${n}]=${file}
  else warn "No body files to plot."; fi
  file="${exp}-reduced/${exp}_outFract.csv"
  if test -e ${file}; then outFract_files[${n}]=${file}
  else warn "No OutFract files to plot."; fi
  file="${exp}-reduced/${exp}_extRatio.csv"
  if test -e ${file}; then extRatio_files[${n}]=${file}
  else warn "No ExtRatio files to plot.";  fi
  file="${exp}-reduced/${exp}_extra-avg.csv"
  if test -e ${file}; then extra_files[${n}]=${file}
  else warn "No extraCellular files to plot."; fi
  (( n++ ))
done
echo "Plotting body file"
plot raw ${body_files[@]}
move body
echo "Plotting outfract file"
plot data ${outFract_files[@]}
move outFract
echo "Plotting extratio file"
plot data ${extRatio_files[@]}
move extRatio
echo "Plotting extra file"
plot raw ${extra_files[@]}
move extra

## find directions and bands to compare
dPVbands=$(getBands dPV $1)
echo " exp 1 dPV bands = {${dPVbands}}"
str1="$dPVbands "
array1=($str1)
ndx=1
for exp in ${exps}; do
  subseqbands=$(getBands dPV ${exp})
  str2="$subseqbands "
  array2=($str2)
  if [[ $ndx == 1 ]]; then
	first=$(commonBands ${array1} ${array2})
	fstr="$first "
	farray=($fstr)
  else
	next=$(commonBands ${farray} ${array2})
	nstr="$next "
	narray=($nstr)
  fi
  ((ndx++))
done
cdPVbands=$(echo "${narray[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' ' )
if [[ $cdPVbands == " " ]]; then
	echo "there are no common dPV bands across all experiments"
fi
for band in ${cdPVbands}; do
  cmpBandsInDir "dPV" ${band} ${exps}
done

dCVbands=$(getBands dCV $1)
echo "exp 1 dCV bands = {${dCVbands}}"
str1="$dCVbands "
array1=($str1)
ndx=1
for exp in ${exps}; do
  subseqbands=$(getBands dCV ${exp})
  str2="$subseqbands "
  array2=($str2)
  if [[ $ndx == 1 ]]; then
	first=$(commonBands ${array1} ${array2})
	fstr="$first "
	farray=($fstr)
  else
	next=$(commonBands ${farray} ${array2})
	nstr="$next "
	narray=($nstr)
  fi
  ((ndx++))
done
cdCVbands=$(echo "${narray[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' ' )
if [[ $cdCVbands == " " ]]; then
	echo "there are no common dCV bands across all experiments"
fi
for band in ${cdCVbands}; do
  cmpBandsInDir "dCV" ${band} ${exps}
done

exit 0
