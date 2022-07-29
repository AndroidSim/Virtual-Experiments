#!/bin/sh
CLASSPATH=$CLASSPATH:../build/classes
CLASSPATH=$CLASSPATH:../lib/ecj19.jar
CLASSPATH=$CLASSPATH:../lib/logback-classic-0.9.28.jar
CLASSPATH=$CLASSPATH:../lib/logback-core-0.9.28.jar
CLASSPATH=$CLASSPATH:../lib/slf4j-api-1.6.1.jar
export CLASSPATH

DIR=`date +"%F-%H%M"`
java isl.util.SRControl
exit 0

