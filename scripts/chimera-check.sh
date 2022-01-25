#This script performs chimera checking from OTU output.
#This script also performs abundance filtering.

seq1='../OTU-claster/rep-seqs-OTU-16s-bc07.qza'
seq2='../OTU-claster/rep-seqs-OTU-16s-bc08.qza'
seq3='../OTU-claster/rep-seqs-OTU-16s-bc09.qza'
seq4='../OTU-claster/rep-seqs-OTU-16s-bc10.qza'
seq5='../OTU-claster/rep-seqs-OTU-16s-bc11.qza'

for file in ../OTU-claster/table-OTU-16s-*;
do
	name=$(basename ${file} .qza)

	if [[ "$name" == "table-OTU-16s-bc07" ]] 
	then
	qiime vsearch uchime-denovo \
	--i-table $file \
	--i-sequences $seq1 \
	--o-chimeras chimeras-16s-bc07.qza \
	--o-nonchimeras nonchimeras-16s-bc07.qza \
	--o-stats stats-chimera-16s-bc07.qza 
	
	elif [[ "$name" == "table-OTU-16s-bc08" ]]
        then
        qiime vsearch uchime-denovo \
        --i-table $file \
        --i-sequences $seq2 \
        --o-chimeras chimeras-16s-bc08.qza \
        --o-nonchimeras nonchimeras-16s-bc08.qza \
        --o-stats stats-chimera-16s-bc08.qza 
		
	elif [[ "$name" == "table-OTU-16s-bc09" ]]
        then
        qiime vsearch uchime-denovo \
        --i-table $file \
        --i-sequences $seq3 \
        --o-chimeras chimeras-16s-bc09.qza \
        --o-nonchimeras nonchimeras-16s-bc09.qza \
        --o-stats stats-chimera-16s-bc09.qza 

	elif [[ "$name" == "table-OTU-16s-bc10" ]]
        then
        qiime vsearch uchime-denovo \
        --i-table $file \
        --i-sequences $seq4 \
        --o-chimeras chimeras-16s-bc10.qza \
        --o-nonchimeras nonchimeras-16s-bc10.qza \
        --o-stats stats-chimera-16s-bc10.qza
	
	else [[ "$name" == "table-OTU-16s-bc11" ]]
        qiime vsearch uchime-denovo \
        --i-table $file \
        --i-sequences $seq5 \
        --o-chimeras chimeras-16s-bc11.qza \
        --o-nonchimeras nonchimeras-16s-bc11.qza \
        --o-stats stats-chimera-16s-bc11.qza 

	fi
done

mv nonchimeras-16s-* ./nonchimeras

#Abundance filtering;filtering chimera and bordeline chimera
#Using tables

nonch1='./nonchimeras/nonchimeras-16s-bc07.qza'
nonch2='./nonchimeras/nonchimeras-16s-bc08.qza'
nonch3='./nonchimeras/nonchimeras-16s-bc09.qza'
nonch4='./nonchimeras/nonchimeras-16s-bc10.qza'
nonch5='./nonchimeras/nonchimeras-16s-bc11.qza'

for file in ../OTU-claster/table-OTU-16s-*;
do
	name=$(basename ${file} .qza)

	if [[ "$name" == "table-OTU-16s-bc07" ]]
	then
	qiime feature-table filter-features \
	--i-table $file \
	--m-metadata-file $nonch1 \
	--o-filtered-table table-nochimeric-noborderline-16s-bc07.qza 
	
	elif [[ "$name" == "table-OTU-16s-bc08" ]]
	then
        qiime feature-table filter-features \
        --i-table $file \
        --m-metadata-file $nonch2 \
        --o-filtered-table table-nochimeric-noborderline-16s-bc08.qza 

	elif [[ "$name" == "table-OTU-16s-bc09" ]]
        then
        qiime feature-table filter-features \
        --i-table $file \
        --m-metadata-file $nonch3 \
        --o-filtered-table table-nochimeric-noborderline-16s-bc09.qza 

	elif [[ "$name" == "table-OTU-16s-bc10" ]]
        then
        qiime feature-table filter-features \
        --i-table $file \
        --m-metadata-file $nonch4 \
        --o-filtered-table table-nochimeric-noborderline-16s-bc10.qza 

	else [[ "$name" == "table-OTU-16s-bc11" ]]
        qiime feature-table filter-features \
        --i-table $file \
        --m-metadata-file $nonch5 \
        --o-filtered-table table-nochimeric-noborderline-16s-bc11.qza 

	fi
done

#Abundance filtering;filtering chimera and bordeline chimera
#Using tables

for file in ../OTU-claster/rep-seqs-OTU-16s-*;
do
	name=$(basename ${file} .qza)

	if [[ "$name" == "rep-seqs-OTU-16s-bc07" ]]
	then
	qiime feature-table filter-seqs \
	--i-data $file \
	--m-metadata-file $nonch1 \
	--o-filtered-data rep-seqs-nochimeric-noborderline-16s-bc07.qza
	
	elif [[ "$name" == "rep-seqs-OTU-16s-bc08" ]]
	then
        qiime feature-table filter-seqs \
        --i-data $file \
        --m-metadata-file $nonch2 \
        --o-filtered-data rep-seqs-nochimeric-noborderline-16s-bc08.qza 

	elif [[ "$name" == "rep-seqs-OTU-16s-bc09" ]]
	then
        qiime feature-table filter-seqs \
        --i-data $file \
        --m-metadata-file $nonch3 \
        --o-filtered-data rep-seqs-nochimeric-noborderline-16s-bc09.qza 

	elif [[ "$name" == "rep-seqs-OTU-16s-bc10" ]]
	then
        qiime feature-table filter-seqs \
        --i-data $file \
        --m-metadata-file $nonch4 \
        --o-filtered-data rep-seqs-nochimeric-noborderline-16s-bc10.qza

	else [[ "$name" == "rep-seqs-OTU-16s-bc11" ]]
        qiime feature-table filter-seqs \
        --i-data $file \
        --m-metadata-file $nonch5 \
        --o-filtered-data rep-seqs-nochimeric-noborderline-16s-bc11.qza

	fi
done


#visualization denovo chimera

for file in ./stats-chimera-16s-*;
do
	name=$(basename ${file} .qza)

	qiime metadata tabulate \
	--m-input-file $file  \
	--o-visualization ../visualization/dvch-$name.qzv
done

#visualization abundance filtering

for file in ./table-nochimeric-noborderline-16s-*;
do 
	name=$(basename ${file} .qza)
	
	qiime feature-table summarize \
	--i-table $file \
	--o-visualization ../visualization/abf-$name.qzv
done


