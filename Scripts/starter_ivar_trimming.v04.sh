#! /bin/bash

#© 2021 cutpatel

#v04 will remove read pairs of aberrant amplicons that cross the expected amplicon boundaries (option -f)

# ivar trim
# Usage: ivar trim -i <input.bam> -b <primers.bed> -p <prefix> [-m <min-length>] [-q <min-quality>] [-s <sliding-window-width>]
# Input Options    Description
           # -i    (Required) Sorted bam file, with aligned reads, to trim primers and quality
           # -b    BED file with primer sequences and positions. If no BED file is specified, only quality trimming wi                                        ll be done.
           # -f    Primer pair information file containing left and right primer names for the same amp                                        licon separated by a tab
                 # If provided, reads that do not fall within atleat one amplicon will be ignored prior to primer trim                                        ming.
           # -x    Primer position offset (Default: 0). Reads that occur at the specified offset positions relative to                                         primer positions will also be trimmed.
           # -m    Minimum length of read to retain after trimming (Default: 50% the average length of the first 1000 reads)
           # -q    Minimum quality threshold for sliding window to pass (Default: 20)
           # -s    Width of sliding window (Default: 4)
           # -e    Include reads with no primers. By default, reads with no primers are excluded
           # -k    Keep reads to allow for reanalysis: keep reads which would be dropped by
                 # alignment length filter or primer requirements, but mark them QCFAILInput Options    Description

# Output Options   Description
           # -p    (Required) Prefix for the output BAM file


primerregions=$1
threads=$2
primerpairs=$3

for file in *.gatk.bam
do
echo "Start primer trimming from $file"
name=$(echo "$file" | cut -f 1 -d '.')
#echo $name
ivar trim -i $file -b $primerregions -f $primerpairs -p "$name.gatk.ivar.bam" -e -q 20 -s 4
samtools sort -@ $threads -o "$name.gatk.ivar.sorted.bam" "$name.gatk.ivar.bam"
samtools index -@ $threads "$name.gatk.ivar.sorted.bam"
rm $file
rm "$name.gatk.bai"
rm "$name.gatk.ivar.bam"
mv "$name.gatk.ivar.sorted.bam" "$name.gatk.bam"
mv "$name.gatk.ivar.sorted.bam.bai" "$name.gatk.bai"
done
exit
