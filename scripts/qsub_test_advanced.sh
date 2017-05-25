#!/bin/sh

if (($# != 2));
then
     echo "Error: Please pass NPROCESSORS and DBSUFFIX to this script. Aborting!"
     exit 1
fi

WORKD="$( dirname "${BASH_SOURCE[0]}" )"
scriptname=${WORKD}/master_for_qstub_test_advanced.sh
perlfile=${WORKD}/vntrseek.pl

nprocs=$1
dbsuffix="$2"

perl $perlfile 99 --dbsuffix $dbsuffix --nprocesses $nprocs   # ask for what step needs to be run next

runnext=$?          # assuming die will always return 255 here and nothing in the 0-19 range ....

case  $runnext  in
  0|1|4|10|13|15) 				# multiple
    echo "runnext state ($runnext), running in MULTI mode!"
    qsub -cwd -l mem_total=64G -l h_rt=40:00:00 -pe omp $nprocs -v DBSUFFIX=$dbsuffix,NPROCS=$nprocs,WORKD=$WORKD,runnext=$runnext $scriptname
    ;;
  2|3|5|6|7|8|9|11|12|14|16|17|18|19)      	# single
    echo "runnext state ($runnext), running in SINGLE mode!"
    qsub -cwd -l h_rt=40:00:00 -v DBSUFFIX=$dbsuffix,NPROCS=$nprocs,WORKD=$WORKD,runnext=$runnext $scriptname
    ;;
  20)
    echo "No more qsub jobs needed ($runnext), done!"
    ;;
  *)
    echo "Illegal runnext state ($runnext)!"
    ;;
esac 

