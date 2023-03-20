#This script takes in OTU tables,taxonomy files,metadata files to generate the taxabarplots for them

for file in *-taxonomy.qza;
do
	name=$(basename ${file} .qza)

	if [[ "$name" == "bc07-taxonomy" ]]
	then
		qiime taxa barplot --i-table table-OTU-16s-bc07.qza --i-taxonomy $file --m-metadata-file ../../../metadata/control7.tsv --o-visualization ./visualization1/taxbar-bc07.qzv
	
		#echo $name
	elif [[ "$name" == "bc08-taxonomy" ]]
	then 
		qiime taxa barplot --i-table table-OTU-16s-bc08.qza --i-taxonomy $file --m-metadata-file ../../../metadata/experimental8.tsv --o-visualization ./visualization1/taxbar-bc08.qzv
	
		#echo $name

	elif [[ "$name" == "bc09-taxonomy" ]]
	then
		qiime taxa barplot --i-table table-OTU-16s-bc09.qza --i-taxonomy $file --m-metadata-file ../../../metadata/experimental9.tsv --o-visualization ./visualization1/taxbar-bc09.qzv

		#echo $name

	elif [[ "$name" == "bc10-taxonomy" ]]
	then
		qiime taxa barplot --i-table table-OTU-16s-bc10.qza --i-taxonomy $file --m-metadata-file ../../../metadata/control10.tsv --o-visualization ./visualization1/taxbar-bc10.qzv

		#echo $name

	else
		qiime taxa barplot --i-table table-OTU-16s-bc11.qza --i-taxonomy $file --m-metadata-file ../../../metadata/experimental11.tsv --o-visualization ./visualization1/taxbar-bc11.qzv

		#echo $name

fi

done
