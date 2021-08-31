#! /bin/bash

#© 2021 cutpatel
#v02: changed threshold to 0.9 and depth 20. Added -i for faster header naming. Needs ivar 1.3!

# ivar consensus
# Usage: samtools mpileup -aa -A -d 0 -Q 0 <input.bam> | ivar consensus -p <prefix> 
# Note : samtools mpileup output must be piped into ivar consensus
# Input Options    Description
           # -q    Minimum quality score threshold to count base (Default: 20)
           # -t    Minimum frequency threshold(0 - 1) to call consensus. (Default: 0)
                 # Frequently used thresholds | Description
                 # ---------------------------|------------
                                          # 0 | Majority or most common base
                                        # 0.2 | Bases that make up atleast 20% of the depth at a position
                                        # 0.5 | Strict or bases that make up atleast 50% of the depth at a position
                                        # 0.9 | Strict or bases that make up atleast 90% of the depth at a position
                                          # 1 | Identical or bases that make up 100% of the depth at a position. Will have highest ambiguities
           # -m    Minimum depth to call consensus(Default: 10)
           # -k    If '-k' flag is added, regions with depth less than minimum depth will not be added to the consensus sequence. Using '-k' will override any option specified using -n 
           # -n    (N/-) Character to print in regions with less than minimum coverage(Default: N)
# Output Options   Description
           # -p    (Required) Prefix for the output fasta file and quality fileivar consensus

#$1="0.9"
#$2="20"
#$3="10"
freq=$1
qual=$2
depth=$3

for file in *.gatk.bam
do
echo "Start creating consensus fasta from $file"
name=$(echo "$file" | cut -f 1 -d '.')
id=$(echo "$file" | cut -f 1 -d '_')

#echo $name
samtools mpileup -aa -A -d 0 -Q 0 $file | ivar consensus -p $name -i ""$id" "$name" Threshold_"$freq"_Quality_"$qual"_Depth_"$depth"" -t $freq -m $depth -q $qual
done
exit
