#! /bin/bash


date=$(date '+%Y-%m-%d')
for file in *.fa; do
	while read line;
		do 
		echo $line | grep -v '>' | grep -o "[N]" | sort | uniq -c | sed "s/ N/\tN/g" | paste -; echo $line | grep '>' | tr "\n" "\t"
		done < $file >> ""$date"_fasta_n_content.tsv"
done
exit
