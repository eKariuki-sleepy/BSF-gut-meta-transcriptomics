#This script takes in demux files and denoise.

for file in ../bac-16S/demux_sequences/*.qza;
do
        name=$(basename ${file} .qza) 
        prefix="demux_bac-16S-mapped-merg_bc07_rna"
        name=${name#"prefix"}
        qiime vsearch dereplicate-sequences \
	--i-sequences $file \
	--o-dereplicated-table 1table-vs-$name.qza \
	--o-dereplicated-sequences 1rep-seqs-vs-$name.qza	
done

#mv rep-seqs-vs-* qiime-vsearch
#mv table-vs-* qiime-vsearch

#visualization

