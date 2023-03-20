
#This script takes in our fastq files and sorts the rRNA based on ribosomal subunits.

mkdir bac-16S
mkdir bac-23S
mkdir arc-16S
mkdir arc-23S
mkdir euk-18S
mkdir euk-28S
mkdir rfam-5S
mkdir rfam-5.8S

for file in ./data/*.fq;
do
	name=$(basename ${file} .fq)
	echo $name
	sortmerna --ref ./databases/rRNA_databases/silva-bac-16s-id90.fasta  --reads $file  --aligned ./aligned/bac-16S-mapped-$name\
	 --other ./aligned/unmapped-bac-16S-$name \
	 --fastx ./aligned/unmapped-bac-16S-$name --threads 7
		rm -rfv /home/oscarmwaura/sortmerna/run/kvdb
		mv ./aligned/bac-16S* bac-16S
done

for file in ./aligned/unmapped-bac-16S*;
do
	name=$(basename ${file} .fq)
	prefix="unmapped-bac-16S-"
        newname=${name#"$prefix"}
	echo "${newname}"
	sortmerna --ref ./databases/rRNA_databases/silva-bac-23s-id98.fasta --reads $file --aligned ./aligned/bac-23S-mapped-$newname\
	 --other ./aligned/unmapped-bac-23S-$newname --fastx ./aligned/unmapped-bac-23S-$newname.fq --threads 7
               rm -rfv /home/oscarmwaura/sortmerna/run/kvdb
		mv ./aligned/bac-23* bac-23S
done


for file in ./aligned/unmapped-bac-23*;
do
  	name=$(basename ${file} .fq)
        prefix="unmapped-bac-23S-"
        newname=${name#"$prefix"}
        echo "${newname}"
        sortmerna --ref ./databases/rRNA_databases/silva-arc-16s-id95.fasta --reads $file --aligned ./aligned/arc-16S-mapped-$newname\
         --other ./aligned/unmapped-arc-16S-$newname --fastx ./aligned/unmapped-arc-16S-$newname.fq --threads 7
               rm -rfv /home/oscarmwaura/sortmerna/run/kvdb
                mv ./aligned/arc-16S* arc-16S
done

for file in ./aligned/unmapped-arc-16S*;
do
  	name=$(basename ${file} .fq)
        prefix="unmapped-arc-16S-"
        newname=${name#"$prefix"}
        echo "${newname}"
        sortmerna --ref ./databases/rRNA_databases/silva-arc-23s-id98.fasta --reads $file --aligned ./aligned/arc-23S-mapped-$newname\
         --other ./aligned/unmapped-arc-23S-$newname --fastx ./aligned/unmapped-arc-23S-$newname.fq --threads 7
               rm -rfv /home/oscarmwaura/sortmerna/run/kvdb
                mv ./aligned/arc-23S* arc-23S
done



for file in ./aligned/unmapped-arc-23S*;
do
  	name=$(basename ${file} .fq)
        prefix="unmapped-arc-23S-"
        newname=${name#"$prefix"}
        echo "${newname}"
        sortmerna --ref ./databases/rRNA_databases/silva-euk-18s-id95.fasta --reads $file --aligned ./aligned/euk-18S-mapped-$newname\
         --other ./aligned/unmapped-euk-18S-$newname --fastx ./aligned/unmapped-euk-18S-$newname.fq --threads 7
               rm -rfv /home/oscarmwaura/sortmerna/run/kvdb
                mv ./aligned/euk-18S* euk-18S
done

for file in ./aligned/unmapped-euk-18S*;
do
  	name=$(basename ${file} .fq)
        prefix="unmapped-euk-18S-"
        newname=${name#"$prefix"}
        echo "${newname}"
        sortmerna --ref ./databases/rRNA_databases/silva-euk-28s-id98.fasta --reads $file --aligned ./aligned/euk-28S-mapped-$newname\
         --other ./aligned/unmapped-euk-28S-$newname --fastx ./aligned/unmapped-euk-28S-$newname.fq --threads 7
               rm -rfv /home/oscarmwaura/sortmerna/run/kvdb
                mv ./aligned/euk-28S* euk-28S
done

for file in ./aligned/unmapped-euk-28S*;
do
        name=$(basename ${file} .fq)
        prefix="unmapped-euk-28S-"
        newname=${name#"$prefix"}
        echo "${newname}"
        sortmerna --ref ./databases/rRNA_databases/rfam-5s-database-id98.fasta --reads $file --aligned ./aligned/rfam-5S-mapped-$newname\
         --other ./aligned/unmapped-rfam-5S-$newname --fastx ./aligned/unmapped-rfam-5S-$newname.fq --threads 7
               rm -rfv /home/oscarmwaura/sortmerna/run/kvdb
                mv ./aligned/rfam-5S* rfam-5S
done

for file in ./aligned/unmapped-rfam-5S*;
do
  	name=$(basename ${file} .fq)
        prefix="unmapped-rfam-5S-"
        newname=${name#"$prefix"}
        echo "${newname}"
        sortmerna --ref ./databases/rRNA_databases/rfam-5.8s-database-id98.fasta --reads $file --aligned ./aligned/rfam-5.8S-mapped-$newname\
         --other ./aligned/unmapped-rfam-5.8S-$newname --fastx ./aligned/unmapped-rfam-5.8S-$newname.fq --threads 7
               rm -rfv /home/oscarmwaura/sortmerna/run/kvdb
                mv ./aligned/rfam-5.8S* rfam-5.8S
done



