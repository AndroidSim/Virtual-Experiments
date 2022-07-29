#!/bin/bash
###
## Executes the linkprep.sh and run.sh scripts
##
## Time-stamp: <2018-01-05 09:44:14 gepr>
###
shopt -s extglob

# arguments 1, 2, etc. are the experiment names in the "exp" directory.

if [ "$#" -lt "1" ]
then
	echo "  Usage: start-exps.sh <experiment 1> <experiment 2> ..."
	echo "		* Need at least 1 experiment to run, where the"
	echo " 		  configuration files for the experiments are in the 'exp' directory"	  
	echo "    	* Uses linkprep.sh and run.sh" 
	echo "		* Run this script in islj directory containing bin"
	exit 1
else
	nexp=$#
	echo "number of experiments = $#"
fi

EXPs=$@
for exp in ${EXPs}
do
	if ./bin/linkprep.sh $exp
	then
		echo "working on experiment $exp"
		nohup time ./bin/run.sh >> keep-nohup.out 2>&1
	else
		echo "error in linkprep.sh, run command for ${exp} not executed"
	fi 
done

exit 0
