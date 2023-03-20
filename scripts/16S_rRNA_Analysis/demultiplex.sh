#This script takes in qiime artifacts,barcode sequences and demultiplexes the artifacts

module load qiime2/2020.6	#load qiime

mkdir ../untrimmed
mkdir ../demux_sequences
mkdir ../visualization


file1="../../metadata/control7.tsv"
file2="../../metadata/experimental8.tsv"
file3="../../metadata/experimental9.tsv"
file4="../../metadata/control10.tsv"
file5="../../metadata/experimental11.tsv"

for file in ./*;
do
  	name=$(basename ${file} .qza)
       
	 #demultiplexing

        if ["$name" -eq "bac-16S-mapped-merg_bc07_rna"]
        then
            	qiime cutadapt demux-single --i-seqs $file \
                --m-barcodes-file $file1 --m-barcodes-column BARCODE \
                --p-error-rate 0 --o-per-sample-sequences demux_$name.qza \
                --o-untrimmed-sequences untrimmed_$name.qza --verbose

        elif ["$name" -eq "bac-16S-mapped-merg_bc08_rna"]
        then
            	qiime cutadapt demux-single --i-seqs $file \
                --m-barcodes-file $file2 --m-barcodes-column BARCODE \
                --p-error-rate 0 --o-per-sample-sequences demux_$name.qza \
                --o-untrimmed-sequences untrimmed_$name.qza --verbose

          elif ["$name" -eq "bac-16S-mapped-merg_bc09_rna"]
        then
            	qiime cutadapt demux-single --i-seqs $file \
                --m-barcodes-file $file3 --m-barcodes-column BARCODE \
                --p-error-rate 0 --o-per-sample-sequences demux_$name.qza \
        --o-untrimmed-sequences untrimmed_$name.qza --verbose

        elif ["$name" -eq "bac-16S-mapped-merg_bc10_rna"]
        then
            	qiime cutadapt demux-single --i-seqs $file \
                --m-barcodes-file $file4 --m-barcodes-column BARCODE \
                --p-error-rate 0 --o-per-sample-sequences demux_$name.qza \
                --o-untrimmed-sequences untrimmed_$name.qza --verbose

        else ["$name" -eq "bac-16S-mapped-merg_bc11_rna"]
                qiime cutadapt demux-single --i-seqs $file \
                --m-barcodes-file $file5 --m-barcodes-column BARCODE \
                --p-error-rate 0 --o-per-sample-sequences demux_$name.qza \
                --o-untrimmed-sequences untrimmed_$name.qza --verbose
         #echo $name

        #Visualize the demultiplexed sequences

        qiime demux summarize \
        --i-data demux_$name.qza \
        --o-visualization ../visualization/demux_$name.qzv

        fi

	# Organize the folders
        mv demux* ../demux_sequences
        mv untrimmed* ../untrimmed


done
