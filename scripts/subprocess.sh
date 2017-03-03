#!/bin/bash
# subprocess for fromBam.sh
# reads segments in $bedfile from $bamfile and pipes to $trf
spid=${1}
bamfile=${2}
bedfile=${3}
output_dir=${4}
trf_exe=$5
trf_param=${6}
trf_out=${7}

umapped_template="*"
trap 'jobs -p | xargs kill' EXIT

if [ -e "$output_dir/e.$spid" ];
then
    rm ${output_dir}/e.${spid}
fi

if [ -e "$output_dir/o.$spid" ];
then
    rm ${output_dir}/o.${spid}
fi

if [ -e "$output_dir/$trf_out-$spid" ];
then
    rm ${output_dir}/${trf_out}-${spid}
fi

postfix=0
while read chr start end
do
    # region by region
    echo ${chr}:${start}-${end} " @ " ${spid}
    if [ "$chr" = "$umapped_template" ]; then
        samtools view -f 4 ${bamfile} | awk '{print ">" $1 "\n" $10}' | ./${trf_exe} - ${trf_param} -ngs >> ${output_dir}/${trf_out}-${spid} 2>> ${output_dir}/e.$spid
    else
        samtools view ${bamfile} ${chr}:${start}-${end} | awk '{print ">" $1 "\n" $10}'  | ./${trf_exe} - ${trf_param} -ngs >> ${output_dir}/${trf_out}-${spid} 2>> ${output_dir}/e.${spid}
    fi
    postfix=$((postfix+1))
done < ${bedfile}

while pgrep -P "$BASHPID" > /dev/null; do
    wait
done