#!/bin/sh

if [ -z "$ARG1" ];
then
     echo "Error: dbsufix is not passed!"
     exit 1
fi

if [ -z "$ARG2" ];
then
     echo "Error: nprocesses is not passed!"
     exit 1
fi

if [ -z "$ARG3" ];
then
     echo "Error: path is not passed!"
     exit 1
fi

batchname=${ARG3}/qsub_test_advanced.sh
scriptname=${ARG3}/master_for_qstub_test_advanced.sh
perlfile=${ARG3}/vntrseek.pl

perl $perlfile 99 --dbsuffix $ARG1 --nprocesses $ARG2   # ask for what step needs to be run next

runnext=$?          # assuming die will always return 255 here and nothing in the 0-19 range ....


if [ $runnext -gt 19 ]; then
  echo "Nothing to do, done!"
elif [ $runnext -gt 15 ]; then
  echo "Doing steps 16-19!"             # single
  perl $perlfile 16 19 --dbsuffix $ARG1 --nprocesses $ARG2 &
  wait
  sh $batchname $ARG1 $ARG2
elif [ $runnext -gt 14 ]; then
  echo "Doing steps 15!"                # multiple
  perl $perlfile 15 15 --dbsuffix $ARG1 --nprocesses $ARG2 &
  wait
  sh $batchname $ARG1 $ARG2
elif [ $runnext -gt 13 ]; then
  echo "Doing steps 14!"                # single
  perl $perlfile 14 14 --dbsuffix $ARG1 --nprocesses $ARG2 &
  wait
  sh $batchname $ARG1 $ARG2
elif [ $runnext -gt 12 ]; then
  echo "Doing steps 13!"                # multiple
  perl $perlfile 13 13 --dbsuffix $ARG1 --nprocesses $ARG2 &
  wait
  sh $batchname $ARG1 $ARG2
elif [ $runnext -gt 10 ]; then
  echo "Doing steps 11-12!"             # single
  perl $perlfile 11 12 --dbsuffix $ARG1 --nprocesses $ARG2 &
  wait
  sh $batchname $ARG1 $ARG2
elif [ $runnext -gt 9 ]; then
  echo "Doing steps 10!"                # multiple
  perl $perlfile 10 10 --dbsuffix $ARG1 --nprocesses $ARG2 &
  wait
  sh $batchname $ARG1 $ARG2
elif [ $runnext -gt 4 ]; then
  echo "Doing steps 5-9!"               # single
  perl $perlfile 5 9 --dbsuffix $ARG1 --nprocesses $ARG2 &
  wait
  sh $batchname $ARG1 $ARG2
elif [ $runnext -gt 3 ]; then
  echo "Doing steps 4!"                 # multiple
  perl $perlfile 4 4 --dbsuffix $ARG1 --nprocesses $ARG2 &
  wait
  sh $batchname $ARG1 $ARG2
elif [ $runnext -gt 1 ]; then
  echo "Doing steps 2-3!"               # single
  perl $perlfile 2 3 --dbsuffix $ARG1 --nprocesses $ARG2 &
  wait
  sh $batchname $ARG1 $ARG2
else
  echo "Doing steps 0-1!"               # multiple
  perl $perlfile 0 1 --dbsuffix $ARG1 --nprocesses $ARG2 &
  wait
  sh $batchname $ARG1 $ARG2
fi


