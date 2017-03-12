#!/bin/bash
# subprocess for fromBam.sh
# reads segments in $bedfile from $bamfile and pipes to $trf
spid=${1}
bamfile=${2}
bedfile=${3}
output_dir=${4}
trf_param="${5}"
trf2proclu_param="${6}"

umapped_template="*"
trap 'jobs -p | xargs kill' EXIT

if [ -e "$output_dir"/"$spid"-"$spid" ];
then
    rm "$output_dir"/"$spid"-"$spid"
fi

# postfix=0
while read -r chr start end
do
    # region by region
    echo "$chr":"$start"-"$end" " @ " "$spid"
    if [ "$chr" = "$umapped_template" ]; then
        samtools view -f 4 "$bamfile" | awk '{print ">" $1 "\n" $10}'
    else
        samtools view "$bamfile" "$chr":"$start"-"$end" | awk '{print ">" $1 "\n" $10}'
    fi
    # postfix=$((postfix+1))
done < "${bedfile}" | "$trf_param" | "$trf2proclu_param" "$output_dir"/"$spid"-"$spid"

while pgrep -P "$BASHPID" > /dev/null; do
    wait
done