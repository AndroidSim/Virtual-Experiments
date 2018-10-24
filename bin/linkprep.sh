#!/bin/bash
###
## prepares the links FROM the simulation directory TO the ../exp/???? experiment directory
##
## Time-stamp: <2017-09-13 08:34:57 gepr>
###
report() {
  linked=`find ${PWD} -name batch_control.properties`
  linked=`ls -ltra $linked`
  linked=${linked%/*}
  linked=${linked##*/}
  echo "Currently linked to ${linked}"
}
linkup() {
  fullfile=$1
  tgt_dir=$2
  basefile=${fullfile##*/}
  linked=$(find $tgt_dir -name $basefile)
  if [[ $linked == "" ]]
  then
    echo "WARN: $fullfile unneccessary."
    return
  fi
  tgt=$fullfile
  if ! test -e ${tgt}
  then
    echo "WARN: ${tgt} not present in experiment directory."
  else
    rm $linked
    ln -s $tgt $linked
  fi
}

if [ "$#" -lt "1" ]; then report; exit 0; fi

exp=$1

TGT_DIR=./build/classes
if ! test -e ${TGT_DIR}; then echo "You must compile first!"; exit 1; fi
SRC_DIR=${PWD%/*}/exp/${exp}

# properties files
for p in $SRC_DIR/*.properties; do linkup $p $TGT_DIR; done
# javascripts
for s in $SRC_DIR/*.js; do linkup $s $TGT_DIR; done
# hepinit json files
for j in $SRC_DIR/*.json; do linkup $j $TGT_DIR; done

exit 0
