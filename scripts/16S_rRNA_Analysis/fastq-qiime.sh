#!/bin/bash
#this script takes the non chimeric sequences and converts them into fastq sequences

for file in ./nonchimeras/*.qza
do
	name=$(basename ${file} .qza)
	names="${name#nonchimeras-}"
	echo $names
	#qiime tools export --input-path $file --output-path ./sequences_from_artifacts/ --output-format $name.fasta

done

for file in *.qza
do
	name=$(basename ${file} .qza)
	names="${name#nonchimeras-}"; echo $names
done	
