###
## utility functions to be sourced by reduction script
##
## Time-stamp: <2019-11-25 19:31:24 gepr>
###
function warn() { echo $'\n'"!!!Warn: $@"$'\n'; }
function error() { echo $'\n'"!!!Error: $@"$'\n'; exit -1; }

## Whole HepStruct graph processing BEFORE the move to the -reduced subdirectory
function pre_mv_whole() {
  exp=$1
  
  ${SRC_DIR}/coarse-compartment.r ${exp}
  if [ $? -ne 0 ]; then error "Coarse Compartment script failed."; fi

  ${SRC_DIR}/fextract.r ${exp}
  if [ $? -ne 0 ]; then error "Fraction and Extraction Ratio script failed."; fi

  ${SRC_DIR}/reduce-event-data.r ${exp}
  if [ $? -ne 0 ]; then error "Unbanded Event Reduction script failed."; fi
  
  ${SRC_DIR}/eg-preproc.r ${exp}
  if [ $? -ne 0 ]; then error "Unbanded Enzyme Group Reduction script failed."; fi
   
  ${SRC_DIR}/rxnfield.r ${exp}
  if [ $? -ne 0 ]; then error "Reaction Field script failed."; fi

  ${SRC_DIR}/ssgeom.r ${exp}
  if [ $? -ne 0 ]; then error "SS Geometry Reduction script failed."; fi

}

## Whole HepStruct graph processing AFTER the move to the -reduced subdirectory
function post_mv_whole() {
   ${SRC_DIR}/find-xpt.r ${exp}
   if [ $? -ne 0 ]; then warn "dPV ∩ dCV script failed."; fi
}

## Banded HepStruct graph processing BEFORE the move to the -reduced subdirectory
function pre_mv_banded() {
  dir=$1; min=$2; max=$3; exp=$4
  
  ${SRC_DIR}/reduce-event-data-inband.r ${dir} ${min} ${max} ${exp}
  if [ $? -ne 0 ]; then error "Event Reduction ${dir}∈[${min},${max}) script failed."; fi
  
  ${SRC_DIR}/eg-inband.r ${dir} ${min} ${max} ${exp}
  if [ $? -ne 0 ]; then error "Enzyme Group Reduction ${dir}∈[${min},${max}) script failed."; fi
  
  ${SRC_DIR}/inextra-inband.r ${dir} ${min} ${max} ${exp}
  if [ $? -ne 0 ]; then error "Banded Solute & MIT Reduction ${dir}∈[${min},${max}) script failed."; fi

  # note the output redirect
  ${SRC_DIR}/hcounts-inband.r ${dir} ${min} ${max} ${exp} > ${exp}_hcounts-${dir}∈\[${min}\,${max}\).csv
  if [ $? -ne 0 ]; then error "Hepatocyte Counts ${dir}∈[${min},${max}) script failed."; fi
  ## remove the file if the Lobule band has zero Hepatocytes
  file=${exp}_hcounts-${dir}∈\[${min}\,${max}\).csv
  search="${dir}"
  ftext=$(grep $search $file)
  output=${ftext#*)} 
  output=${output:3}
  if [ ${output} == "0" ]; then rm $file; fi
  
  ${SRC_DIR}/dataperH-inband.r ${dir} ${min} ${max} ${exp}
  if [ $? -ne 0 ]; then error "Per-Hepatocyte derivation ${dir}∈[${min},${max}) script failed."; fi

  ${SRC_DIR}/Hcounts-avgsd.r ${dir} ${min} ${max} ${exp}
  if [ $? -ne 0 ]; then error "Hepatocyte Averages ${dir}∈[${min},${max}) script failed."; fi

}

## Banded HepStruct graph processing AFTER the move to the -reduced subdirectory
function post_mv_banded() {
  dir=$1; min=$2; max=$3; exp=$4
  ${SRC_DIR}/calc-clearance.r ${dir} ${min} ${max} ${exp}
  if [ $? -ne 0 ]; then error "Clearance Derivation ${dir}∈[${min},${max}) script failed."; fi
  
  ${SRC_DIR}/calc-exposure.r ${dir} ${min} ${max} ${exp}
  if [ $? -ne 0 ]; then error "Exposure Derivation ${dir}∈[${min},${max}) script failed."; fi
  
  ${SRC_DIR}/calc-totalamt.r ${dir} ${min} ${max} ${exp}
  if [ $? -ne 0 ]; then error "Total Mobile Object Amounts ${dir}∈[${min},${max}) script failed."; fi
  
}
