#This script takes in table.qza(feature table) and rep-seqs.qza(feature sequence data) from vsearch output and outputs visualization for each sample.
#It also takes in the metadata files.
#Provides a .qzv file for visualization in qiime2view for all the files repectively.
#visualization

file1="../metadata/control7.tsv"
file2="../metadata/experimental8.tsv"
file3="../metadata/experimental9.tsv"
file4="../metadata/control10.tsv"
file5="../metadata/experimental11.tsv"

for file in ./qiime-vsearch/*.qza;
do
  	name=$(basename ${file} .qza)

	if [[ "$name" == "table-vs-demux_bac-16S-mapped-merg_bc07_rna" ]]
	then
	qiime feature-table summarize \
	--i-table $file \
	--o-visualization ./visualization/table-vs-$name.qzv \
	--m-sample-metadata-file $file1

	elif [[ "$name" == "table-vs-demux_bac-16S-mapped-merg_bc08_rna" ]]
	then
	qiime feature-table summarize \
	--i-table $file \
	--o-visualization ./visualization/table-vs-$name.qzv \
	--m-sample-metadata-file $file2  

	elif [[ "$name" == "table-vs-demux_bac-16S-mapped-merg_bc09_rna" ]]
	then
	qiime feature-table summarize \
	--i-table $file \
	--o-visualization ./visualization/table-vs-$name.qzv \
	--m-sample-metadata-file $file3  

	elif [[ "$name" == "table-vs-demux_bac-16S-mapped-merg_bc10_rna" ]]
	then
	qiime feature-table summarize \
	--i-table $file \
	--o-visualization ./visualization/table-vs-$name.qzv \
	--m-sample-metadata-file $file4  

	elif [[ "$name" == "table-vs-demux_bac-16S-mapped-merg_bc11_rna" ]]
        then
	qiime feature-table summarize \
	--i-table $file \
	--o-visualization ./visualization/table-vs-$name.qzv \
	--m-sample-metadata-file $file5
	
	#Takes in the feature sequence data.
        else [[ "$name" == "rep-seqs-vs-demux_bac-16S-mapped-merg_bc07_rna" ]]
	qiime feature-table tabulate-seqs \
	--i-data $file \
	--o-visualization ./visualization/rep-seqs-vs-$name.qzv

	fi
done

