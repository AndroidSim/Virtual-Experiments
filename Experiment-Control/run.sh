#!/bin/sh

###
## Used to execute the ISLJ.  E.g.
## nohup time ./bin/run.sh > keep-nohup.out 2>&1 &
##
## Time-stamp: <2018-02-06 16:29:19 gepr>
###

CLASSPATH=$CLASSPATH:./build/classes
#CLASSPATH=$CLASSPATH:./cfg
CLASSPATH=$CLASSPATH:./lib/colt-1.2.0.jar
CLASSPATH=$CLASSPATH:./lib/commons-math3-3.5.jar
CLASSPATH=$CLASSPATH:./lib/ecj-23-util.jar
CLASSPATH=$CLASSPATH:./lib/genson-1.3.jar
CLASSPATH=$CLASSPATH:./lib/iText-4.2.0.jar
CLASSPATH=$CLASSPATH:./lib/jcommon-1.0.23.jar
CLASSPATH=$CLASSPATH:./lib/jfreechart-1.0.19.jar
CLASSPATH=$CLASSPATH:./lib/logback-classic-1.1.3.jar
CLASSPATH=$CLASSPATH:./lib/logback-core-1.1.3.jar
CLASSPATH=$CLASSPATH:./lib/mason19-with-src.jar
CLASSPATH=$CLASSPATH:./lib/slf4j-api-1.7.12.jar
CLASSPATH=$CLASSPATH:./lib/junit-4.12.jar
#CLASSPATH=$CLASSPATH:./lib/bsg.jar
CLASSPATH=$CLASSPATH:./lib/bsg-r1437.jar
export CLASSPATH

DATE=$(date +"%F-%H%M")
NANO=$(( $(date +"%-N") / 100000 ))
DIR=$DATE"-"$NANO
mkdir $DIR 
time java -Duser.language=US -Duser.country=US isl.Main $* > $DIR/log 2>&1 
cp -rpL ./build/classes/cfg $DIR/
mv logs/* $DIR/
bzip2 -9 $DIR/log
exit 0

