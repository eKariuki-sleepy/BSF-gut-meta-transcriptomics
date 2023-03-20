#This scripts takes in the representative sequences and generates the taxonomy

module load qiime2/2020.6

mkdir visualization1
for file in rep-seqs-OTU-16s-*.qza;
do
	name=$(basename ${file} .qza)
	names=${name#"rep-seqs-OTU-16s-"}

	qiime feature-classifier classify-sklearn \
  	--i-classifier classifier.qza \
  	--i-reads $file \
  	--o-classification $names-taxonomy.qza 

	qiime metadata tabulate \
  	--m-input-file $names-taxonomy.qza \
  	--o-visualization $names-taxonomy.qzv 
	#echo $names
done

mv *.qzv ./visualization1
