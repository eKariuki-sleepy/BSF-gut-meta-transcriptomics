#This script takes in a fastq file gzip it and converts it into a qiime artifact

mkdir artifacts
gzip *.fq
	
#Convert the file into a qiime artifact

qiime tools import --type MultiplexedSingleEndBarcodeInSequence --input-path ./bac-16S-mapped-merg_bc07_rna.fq.gz \
--output-path ./artifacts/bac-16S-mapped-merg_bc07_rna.qza  

qiime tools import --type MultiplexedSingleEndBarcodeInSequence --input-path ./bac-16S-mapped-merg_bc08_rna.fq.gz \
--output-path ./artifacts/bac-16S-mapped-merg_bc08_rna.qza

qiime tools import --type MultiplexedSingleEndBarcodeInSequence --input-path ./bac-16S-mapped-merg_bc10_rna.fq.gz \
--output-path ./artifacts/bac-16S-mapped-merg_bc10_rna.qza

qiime tools import --type MultiplexedSingleEndBarcodeInSequence --input-path ./bac-16S-mapped-merg_bc09_rna.fq.gz \
--output-path ./artifacts/bac-16S-mapped-merg_bc09_rna.qza

qiime tools import --type MultiplexedSingleEndBarcodeInSequence --input-path ./bac-16S-mapped-merg_bc11_rna.fq.gz \
--output-path ./artifacts/bac-16S-mapped-merg_bc11_rna.qza



