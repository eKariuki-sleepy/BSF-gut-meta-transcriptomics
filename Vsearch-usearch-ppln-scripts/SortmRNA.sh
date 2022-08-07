### 5. Filtering and sorting rRNA Sequences with SortMeRNA

##### 5.1 Installation
```
conda install sortmerna
```
##### 5.2 Preparation
Before we can run the sortmerna command, we must first download the eukaryotic, archeal and bacterial rRNA databases. The sortmerna_db/ folder will be the location that we will keep the files necessary to run SortMeRNA. These databases only need to be created once, so any future RNAseq experiements can use these files.

##### 5.3 Download the sortmerna package (~2min) into sortmerna_db folder
```
    wget -P sortmerna_db https://github.com/biocore/sortmerna/archive/2.1b.zip

    # Decompress folder 
    unzip sortmerna_db/2.1b.zip -d sortmerna_db

    # Move the database into the correct folder
    mv sortmerna_db/sortmerna-2.1b/rRNA_databases/ sortmerna_db/

    # Remove unnecessary folders
    rm sortmerna_db/2.1b.zip
    rm -r sortmerna_db/sortmerna-2.1b

```
##### 5.4 Sorting 
```
**id85- < 50k reads
#This script takes in our fastq files and sorts the rRNA based on ribosomal subunits.

mkdir bac-16S
mkdir aligned

#The loop generates mapped rRNA reads ,unmapped rRNA reads with respect to each sample and a log file.
#Removes the kvdb file to allows next run.(sortmerna guidlines)

for file in ./*.fq;
do

        name=$(basename ${file} .fq)
       # echo $name
        sortmerna --ref ../sortmerna/sortmerna_db/rRNA_databases/silva-bac-16s-database-id85.fasta   --reads $file \
        --aligned ./aligned/bac-16S-85-mapped$name  --other ./bac-16S-85-umapped$name \
        --fastx ./bac-16S-85-unmapped$name  --threads 9
             rm -rfv /home/oscarmwaura/sortmerna/run/kvdb
             mv ./aligned/bac-16S* bac-16S

done


**id90 > 50k reads.
#This script takes in our fastq files and sorts the rRNA based on ribosomal subunits.

mkdir bac-16S
mkdir aligned

#The loop generates mapped rRNA reads ,unmapped rRNA reads with respect to each sample and a log file.
#Removes the kvdb file to allows next run.(sortmerna guidlines)

for file in ./*.fq;
do

        name=$(basename ${file} .fq)
       # echo $name
        sortmerna --ref ../sortmerna/sortmerna_db/rRNA_databases/silva-bac-16s-database-id85.fasta   --reads $file \
        --aligned ./aligned/bac-16S-85-mapped$name  --other ./bac-16S-85-umapped$name \
        --fastx ./bac-16S-85-unmapped$name  --threads 9
             rm -rfv /home/oscarmwaura/sortmerna/run/kvdb
             mv ./aligned/bac-16S* bac-16S

done
```
