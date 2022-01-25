#This script converts the otu tables into visualizable format

for file in ./table-OTU-16s-*;
do
	name=$(basename ${file} .qza)
	
	qiime tools export --input-path $file --output-path ./visualization1/ 	
	cd ./visualization1/
	biom convert -i *biom -o $name.txt --to-tsv
	rm *.biom
	#echo $name
	cd ..

	#echo "done"
done
