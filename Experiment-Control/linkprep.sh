#!/bin/bash
###
## prepares the links FROM the simulation directory TO the ../exp/???? experiment directory
##
## Time-stamp: <2020-07-21 17:35:37 gepr>
###
function report() {
  linked=`find ${PWD} -name batch_control.properties`
  linked=`ls -ltra $linked`
  linked=${linked%/*}
  linked=${linked##*/}
  echo "Currently linked to ${linked}"
}
function linkup() {
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
if [ ! -d ${SRC_DIR} ]; 
then 
	echo "ERROR: the experiment directory ${SRC_DIR} does not exist"; 
	exit 1; 
else
	if [ ! "$(ls -A $SRC_DIR)" ];
	then
		echo "ERROR: the experiment directory ${SRC_DIR} is empty";
		exit 1;
	else
		if [ ! "$(ls $SRC_DIR/*.properties)" ];
		then
			echo "ERROR: the experiment directory ${SRC_DIR} contains no .properties files";
			exit 1;
		fi
	fi
fi

# properties files
for p in $SRC_DIR/*.properties; do linkup $p $TGT_DIR; done
# hepinit json files
for j in $SRC_DIR/*.json; do linkup $j $TGT_DIR; done

exit 0
