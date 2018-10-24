#!/bin/bash
###
## Creates a symbolic link between an experiment directory and a new 
## directory to be used in the compare-exps.sh script.
##
## Time-stamp: <2018-03-28 aks>
###
shopt -s extglob

# arguments in pairs:
# argument 1 = experiment directory
# argument 2 = symbolic link name

if [ $# -lt 2 ]
then
  echo "Usage: link-exps.sh <experiment directory> <link name> ..."
  echo "	* Need 2 arguments, 1st = output dir, 2nd = experiment name"
  echo "  e.g. link-exps.sh 2015-03-* exp000 ..."
  echo "  would create a symbolic link like \"exp000x -> 2015-03-12-1293\""
  exit 1
fi

while (( $# ))
do
	expname=$1
	linkname=$2
	ln -s ${expname} ${linkname}
	shift 2
done

exit 0
