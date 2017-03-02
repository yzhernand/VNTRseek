#!/bin/bash

# Script to execute

# Nu
procs=$1
bamfile=$2
output_dir=$3
trf=$4
match=$5
mismatch=$6
delta=$7
PM=$8
PI=$9
minscore=$10
maxperiod=$11


trap 'jobs -p | xargs kill' EXIT  

# make bed files
samtools idxstats $bamfile | java MakeBedFiles $procs $output_dir

for spid in $(seq 1 $procs)
do    
	trf_out=$12
    #run subprocess
    ./subprocess.sh $spid $bamfile $output_dir/bed$spid $output_dir $trf $match $mismatch $delta $PM $PI $minscore $maxperiod $trf_out & 
done

rm -I output_dir/bed*
rm -I output_dir/o.*
rm -I output_dir/e.*

while pgrep -P "$BASHPID" > /dev/null; do
  wait
done
