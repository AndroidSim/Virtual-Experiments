#!/bin/bash
user=bsg
ws=islj-ws
id=ebs-2013-03-11.pem
machines=`cat ~/$ws/machines`

if [ "$#" -lt "1" ]; then echo "You must include an experiment directory."; exit 1; fi

exp=$1

# distribute
count=0
for machine in $machines
do
  # copy experiment
  echo "Sending $exp to $user@$machine:~/$ws/exp"
  scp -i $id -rp exp/$exp $user@$machine:~/$ws/exp/

  # link
  echo "Linking $exp to individual parameter files on $machine."
  ssh -i $id $user@$machine "cd $ws/islj; ./bin/linkprep.sh $exp"

  # launch
  echo "Executing."
  ssh -i $id $user@$machine "cd $ws/islj; nohup time ./bin/run.sh > nohup.out 2>&1 " &
  proc[$count]="$procs $!"
  count=$(($count+1))
done

# poll
finished=0
while (( finished == 0 ))
do
  stop=1
  count=0
  for machine in $machines
  do 
    ps ${proc[$count]}
    if (( $? == 0 )); then stop=0; fi
    count=$(($count+1))
  done
  if (( $stop == 1 )); then finished=1; fi
  sleep 10
done

# harvest and cleanup
echo "Harvesting."
restop=results/$exp
mkdir -p $restop
for machine in $machines
do
  # harvest
  resdir=$restop/$machine
  mkdir $resdir
  scp -i $id -rp $user@$machine:~/$ws/islj/2013-* $resdir
  # cleanup
  ssh -i $id $user@$machine "cd ~/$ws; rm -rf exp/$exp; rm -rf islj/2013-*"
done

