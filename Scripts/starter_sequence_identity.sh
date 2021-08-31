#! /bin/bash
#© 2021 cutpatel

date=$(date '+%Y-%m-%d')
#reference="coronavirus2_wuhan-hu-1_MN908947.3_2020-04-29"
#reference_fasta="$HOME/Data/Ref/"$reference".fasta"
reference_fasta=$1

echo SampleID$'\t'%Identity >""$date"_identity.tsv"
for fasta in *.fa; do
echo "Start calculating sequence identity of $fasta to $reference"
name=$(echo "$fasta" | cut -f 1 -d '.')
id=$(echo "$fasta" | cut -f 1 -d '_')
#echo $id
#echo $name
export id


###from https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
minimap2 -c $reference_fasta $fasta \
  | perl -ane 'if(/tp:A:P/&&/NM:i:(\d+)/){$n+=$1;$m+=$1 while/(\d+)M/g;$g+=$1,++$o while/(\d+)[ID]/g}END{print("$ENV{id}\t",1-(($n-$g+$o)/($m+$o)),"\n")}' >>""$date"_identity.tsv"

done
unset id
exit

