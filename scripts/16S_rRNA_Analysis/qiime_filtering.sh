#this script does quality control

module load qiime2/2020.6
mkdir filtered
mkdir denoise

for file in ./demux_sequences/*.qza;
do
	name=$(basename ${file} .qza)
	prefix="demux_bac-16S-mapped-merge_bc07_rna"
	name=${name#"$prefix"}
	qiime quality-filter q-score \
	--i-demux $file \
	--o-filtered-sequences demux-filtered-$name.qza \
	--o-filter-stats demux-filter-stats-$name.qza

	qiime metadata tabulate \
	--m-input-file demux-filter-stats-$name.qza \
	--o-visualization demux-filter-stats-$name.qzv
done

mv demux-filtered* filtered
mv *.qzv visualization

#This script takes in demux files and denoise.

for file in ./filtered/*.qza;
do
	name=$(basename ${file} .qza)
	prefix="demux-filtered-demux_bac-16S-mapped-merg_bc07_rna"
	name=${name#"prefix"}
	qiime dada2 denoise-single \
  	--p-trim-left 0 \
  	--p-trunc-len 0 \
  	--i-demultiplexed-seqs denoise-$name.qza \
  	--o-representative-sequences rep-seqs-$name.qza \
  	--o-table denoise-table-$name.qza \
  	--o-denoising-stats denoise-stats-$name.qza
done

mv rep-seqs* denoise
mv denoise-table* denoise
mv denoise-stats* denoise
