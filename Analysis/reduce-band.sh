#!/bin/bash
###
## Executes the data reduction functions on the given band.
##
## Time-stamp: <2019-11-25 19:31:16 gepr>
###
shopt -s extglob

SRC_DIR=$(dirname ${BASH_SOURCE[0]})
source ${SRC_DIR}/redlib.sh

# tests to see if the argument is an integer
function is_int() { return $(test "$@" -eq "$@" > /dev/null 2>&1); }
function usage() {
   echo "  Usage: reduce-band.sh d[CP]V min max <experiment 1> <experiment 2> ..."
   exit 1
}

# arguments 1, 2, etc. are the experiment names.

if [ $# -lt 4 ]; then usage; else
  dir=$1; shift;
  if $(is_int $1); then min=$1; shift; else echo "$1 is not a valid minimum."; usage; fi
  if $(is_int $1); then max=$1; shift; else echo "$1 is not a valid maximum."; usage; fi
  EXPs=$@
  echo "Reducing ${dir}∈[${min},${max}) for ${EXPs}"
fi

# data reduction and extraction

EXPs=$@
for exp in ${EXPs}
do
  cleanexp=$(dirname ${exp})/$(basename ${exp}) # clean up trailing slash for file names
  echo "working on experiment ${cleanexp}:"

  ## eg-inband.r needs eg-preproc.r
  if ! test -e "./${cleanexp}_enzymes-${dir}.csv"; then ${SRC_DIR}/eg-preproc.r ${cleanexp}; fi
  pre_mv_banded ${dir} ${min} ${max} ${cleanexp}
  
  ## test for -reduced directory and if its not present, don't run post_mv_banded()
  runpost=1
  if ! test -e "./${cleanexp}-reduced"; then
    runpost=0
    mkdir -p ${cleanexp}-reduced
  fi
  mv ${cleanexp}*.csv ${cleanexp}-reduced
  if (( $runpost == 1 )); then post_mv_banded ${dir} ${min} ${max} ${exp}; fi
done

exit 0
