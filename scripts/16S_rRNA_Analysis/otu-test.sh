#This scripts performs OTU clustering.
#It takes in feature table and feature sequences from vsearch output.
#It takes in reference 16s-OTUs from greengenes database.

seq1='../../qiime-vsearch/rep-seqs-vs-demux_bac-16S-mapped-merg_bc07_rna.qza'
seq2='../../qiime-vsearch/rep-seqs-vs-demux_bac-16S-mapped-merg_bc08_rna.qza'
seq3='../../qiime-vsearch/rep-seqs-vs-demux_bac-16S-mapped-merg_bc09_rna.qza'
seq4='../../qiime-vsearch/rep-seqs-vs-demux_bac-16S-mapped-merg_bc10_rna.qza'
seq5='../../qiime-vsearch/rep-seqs-vs-demux_bac-16S-mapped-merg_bc11_rna.qza'

for file in ../../qiime-vsearch/f-tables/*.qza;
do
	name=$(basename ${file} .qza)	

	if [[ "$name" == "table-vs-demux_bac-16S-mapped-merg_bc07_rna" ]]
	then
	qiime vsearch cluster-features-open-reference \
	--i-table $file \
	--i-sequences $seq1 \
	--i-reference-sequences ../../silva-db/sample_data/silva_132_99.qza \
	--p-perc-identity 0.99 \
	--o-clustered-table table-OTU-16s-bc07.qza \
	--o-clustered-sequences rep-seqs-OTU-16s-bc07.qza \
	--o-new-reference-sequences new-ref-seqs-16s-bc07.qza \
 	--p-threads 8 #number of cores for processing. 
	
	elif [[ "$name" == "table-vs-demux_bac-16S-mapped-merg_bc08_rna" ]]
        then
	qiime vsearch cluster-features-open-reference \
        --i-table $file \
        --i-sequences $seq2 \
        --i-reference-sequences  ../../silva-db/sample_data/silva_132_99.qza \
        --p-perc-identity 0.99 \
        --o-clustered-table table-OTU-16s-bc08.qza \
        --o-clustered-sequences rep-seqs-OTU-16s-bc08.qza \
        --o-new-reference-sequences new-ref-seqs-16s-bc08.qza \
        --p-threads 8 #number of cores for processing.

	elif [[ "$name" == "table-vs-demux_bac-16S-mapped-merg_bc09_rna" ]]
        then
	qiime vsearch cluster-features-open-reference \
        --i-table $file \
        --i-sequences $seq3 \
        --i-reference-sequences  ../../silva-db/sample_data/silva_132_99.qza \
        --p-perc-identity 0.99 \
        --o-clustered-table table-OTU-16s-bc09.qza \
        --o-clustered-sequences rep-seqs-OTU-16s-bc09.qza \
        --o-new-reference-sequences new-ref-seqs-16s-bc09.qza \
        --p-threads 8 #number of cores for processing.

	elif [[ "$name" == "table-vs-demux_bac-16S-mapped-merg_bc10_rna" ]]
        then
	qiime vsearch cluster-features-open-reference \
        --i-table $file \
        --i-sequences $seq4 \
        --i-reference-sequences  ../../silva-db/sample_data/silva_132_99.qza \
        --p-perc-identity 0.99 \
        --o-clustered-table table-OTU-16s-bc10.qza \
        --o-clustered-sequences rep-seqs-OTU-16s-bc10.qza \
        --o-new-reference-sequences new-ref-seqs-16s-bc10.qza \
        --p-threads 8 #number of cores for processing.

	else [[ "$name" == "table-vs-demux_bac-16S-mapped-merg_bc11_rna" ]]
	qiime vsearch cluster-features-open-reference \
        --i-table $file \
        --i-sequences $seq5 \
        --i-reference-sequences  ../../silva-db/sample_data/silva_132_99.qza \
        --p-perc-identity 0.99 \
        --o-clustered-table table-OTU-16s-bc11.qza \
        --o-clustered-sequences rep-seqs-OTU-16s-bc11.qza \
        --o-new-reference-sequences new-ref-seqs-16s-bc11.qza \
        --p-threads 8 #number of cores for processing.

	fi
done
