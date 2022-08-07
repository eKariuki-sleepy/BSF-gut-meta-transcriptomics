#Takes in the the rep-sequence

for file in new-ref-seqs*;
do
	name=$(basename ${file} .qza)
	
        qiime feature-table tabulate-seqs \
        --i-data $file \
        --o-visualization ./visualization1/$name.qzv
done


