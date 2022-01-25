mkdir visual_untrimmed

for file in *.qza;
do
	name=$(basename ${file} .qza)
	qiime demux summarize --i-data $file --o-visualization ./visual_untrimmed/v_$name.qzv
done
