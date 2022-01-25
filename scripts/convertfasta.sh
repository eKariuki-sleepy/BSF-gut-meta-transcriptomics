#This script converts fasta file into tsv

for file in ./*.fasta;
do
	
	name=$(basename ${file} .fasta)

	awk 'BEGIN{RS=">"}{print "#"$1"\t"$2;}' $file | tail -n+2 > $name.tsv
done
