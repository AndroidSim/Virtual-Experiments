#!/bin/bash
###
## Renames the files in the output directory from reduce-exps.sh  
## to another name and moves them to another directory with a name that
## is a concatenation of the raw output directory (date-number format) 
## and the new name.
##
## Time-stamp: <2018-04-13 aks>
###
shopt -s extglob

# arguments in pairs:
# argument 1 = raw output directory (i.e. date-number format)
# argument 2 = new name, e.g. experiment name in "exp" directory

if [ $# -lt 2 ]
then
	echo "  Usage: rename-exps.sh <output dir 4 exp 1> <new exp 1 name> ..."
	echo "		* Need 2 arguments, 1st = output dir, 2nd = experiment name"
	echo "		* Run this script in islj directory containing bin"
	exit 1
fi

# parse arguments in pairs, 1st = output dir, 2nd = new experiment name,
# copy or move reduced files from reduce-exps.sh to new experiment name
while (( $# ))
do
	rawout=$1
	expname=$2
	
	# check if reduced files exist
	files=${rawout}*.csv
	afe="F"
	for f in $files
	do
		if [ -e $f ]
		then
			afe="T"
			break
		fi
	done
	# if reduced directory exists and reduced files don't exist in pwd
	if [ -d ${rawout}"-reduced" ] && [ $afe == "F" ]
	then
		echo "reduced dir exists and reduced files don't exist in pwd for ${rawout}"
		cd ${rawout}"-reduced"
		files=${rawout}*.csv
		for f in $files
		do
			# copy raw file to new file (i.e. rename & keep original)
			outfile=${f/${rawout}/${expname}}
			cp ${f} ${outfile}
		done
		# move new files to directory
		cd ..
		mkdir -p ${rawout}_${expname}
		mv ${rawout}"-reduced"/${expname}*.csv ${rawout}_${expname}
	# if reduced directory exists and reduced files exist in pwd
	elif [ -d ${rawout}"-reduced" ] && [ $afe == "T" ]
	then
		echo "both reduced dir and reduced files exist in pwd for ${rawout}"
		for f in $files
		do
			# move raw file to new file (i.e. rename & delete original)
			outfile=${f/${rawout}/${expname}}
			mv ${f} ${outfile}
		done
		cd ${rawout}"-reduced"
		afiles=${rawout}*.csv
		for f in $afiles
		do
			# copy raw file to new file (i.e. rename & keep original)
			outfile=${f/${rawout}/${expname}}
			cp ${f} ${outfile}
		done
		# move new files to directory
		cd ..
		mkdir -p ${rawout}_${expname}
		mv ${rawout}"-reduced"/${expname}*.csv ${rawout}_${expname}
		mv ${expname}*.csv ${rawout}_${expname}
	# if reduced directory doesn't exist and reduced files exist in pwd
	elif [ ! -d ${rawout}"-reduced" ] && [ $afe == "T" ]
	then
		echo "reduced dir doesn't exist and reduced files exist in pwd for ${rawout}"
		for f in $files
		do
			# move raw file to new file (i.e. rename & delete original)
			outfile=${f/${rawout}/${expname}}
			mv ${f} ${outfile}
		done
		# move new files to directory
		mkdir -p ${rawout}_${expname}
		mv ${expname}*.csv ${rawout}_${expname}
	else
		echo "neither reduced dir or files exist in pwd for ${rawout}"
	fi
	
	shift 2
done

# using the command "rename" is an option; however, it can be different
# for different Linux distributions
#rename 's/$1/$2/g' *.csv

exit 0
