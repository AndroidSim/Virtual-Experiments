#!/bin/bash 

## This script is used to sample a certain number of MC trials across
# a set of experiment directories.
#
# command-line arguments:	1 = number of MC trials in experiment dirs
#							2 = the file containing the set of 
#								experiment directories
##

# Program logic:
# 1) user inputs nMC and file containing exp dirs
# 2) randomly select one exp dir (use dir only once?)
# 3) from that exp dir, randomly select one MC trial, which contains
#    6 files (i.e. body, dead, enzymes, run, mobileObject_zone_1_2, and 
#    mobileObject_zone_3) (use MC trial # only once?)
# 4) repeat "number of MC trials" times
# 5) put all sampled files into a directory

# check for the right number of arguments
if [[ $# < 2 ]]
then
	echo "at least 2 arguments are needed, the first is the number
	      of MC trials in the experiment directories, and the second
	      argument is a file that vertically contains the experiment 
	      directories from which to sample"
	exit 1
fi

nargs=$#
nMC=$1
#nx=$(($# - 1))

# shift args left to just have exp dirs in array
#shift
#exps=("$@")

# read experiments directories from file and store in array
xsfile=$2
x=0
while read line
do
	#echo $line
	exps[$x]="$line"
	((x++))
done < $xsfile
nx=${#exps[@]}
#echo $nx
#echo ${exps[@]}

# check if experiment directories are directories
for e in ${exps[@]}
do
	if [ ! -d ${e} ]
	then
		echo "experiment ${e} is not a directory"
		exit 1
	fi
done

outdir=${PWD}/"analog-sampling-output"
if [ ! -d ${outdir} ]
then
	mkdir $outdir
fi

# loop nMC times to sample experiment directories
for i in $(seq 1 $nMC) 
do
# 	randomly select an experiment directory
# 	using bash shell function "shuf", RANDOM is an alternative
	inx=$(($nx-1))
	rnx=$(shuf -i 0-$inx -n 1)
	x=${exps[$rnx]}
# 	randomly select a MC trial from that experiment directory
# 	MC trial = body, dead, enzymes, run, mobileObject_zone_1_2, and 
#			   mobileObject_zone_3 files
	int=$((nMC-1))
	rnt=$(shuf -i 0-$int -n 1)
#	format numbers less than 10
	if [[ $rnt < 10 ]]
	then
		rnt=0${rnt}
	fi

#	find the right MC trial files, copy to output dir, and rename them
	files=$(find ${PWD}/${x} -name "*${rnt}.csv")
	files="$files"
	for f in $files
	do
		cp $f $outdir
		of=${outdir}/$(basename ${f})
		rf=${of/.csv/${i}.csv}
		mv $of $rf
	done
done

exit 0
