#!/bin/bash
spid=$1
bamfile=$2
bedfile=$3
output_dir=$4
trf=$5
match=$6
mismatch=$7
delta=$8
PM=$9
PI=$10
minscore=$11
maxperiod=$12

trf_out=$13

umapped_template="*"
trap 'jobs -p | xargs kill' EXIT

#mkfifo $output_dir/pipe$spid
while read chr start end
do
    # region by region
    if [ "$chr" = "$umapped_template" ]; then
        samtools view -f 4 $bamfile | awk '{print ">" $1 "\n" $10}' | ./$trf - $match $mismatch $delta $PM $PI $minscore $maxperiod -h -ngs > $out_dir/$trf_out
    else
        samtools view $bamfile $chr:$start-$end | awk '{print ">" $1 "\n" $10}' | ./$trf - $match $mismatch $delta $PM $PI $minscore $maxperiod -h -ngs > $out_dir/$trf_out
    fi
done < $bedfile
#rm $output_dir/pipe$spid

while pgrep -P "$BASHPID" > /dev/null; do
    wait
    rm pipe$spid
done