#! /bin/bash

#© 2021 cutpatel

# ivar trim
# Usage: ivar trim -i <input.bam> -b <primers.bed> -p <prefix> [-m <min-length>] [-q <min-quality>] [-s <sliding-window-width>]
# Input Options    Description
           # -i    (Required) Sorted bam file, with aligned reads, to trim primers and quality
           # -b    (Required) BED file with primer sequences and positions
           # -m    Minimum length of read to retain after trimming (Default: 30)
           # -q    Minimum quality threshold for sliding window to pass (Default: 20)
           # -s    Width of sliding window (Default: 4)
           # -e    Include reads with no primers. By default, reads with no primers are excluded
# Output Options   Description
           # -p    (Required) Prefix for the output BAM file


primerregions=$1
threads=$2

for file in *.gatk.bam
do
echo "Start primer trimming from $file"
name=$(echo "$file" | cut -f 1 -d '.')
#echo $name
ivar trim -i $file -b $primerregions -p "$name.gatk.ivar.bam" -e -q 20 -s 4 -m 30
samtools sort -@ $threads -o "$name.gatk.ivar.sorted.bam" "$name.gatk.ivar.bam"
samtools index -@ $threads "$name.gatk.ivar.sorted.bam"
rm $file
rm "$name.gatk.bai"
rm "$name.gatk.ivar.bam"
mv "$name.gatk.ivar.sorted.bam" "$name.gatk.bam"
mv "$name.gatk.ivar.sorted.bam.bai" "$name.gatk.bai"
done
exit
