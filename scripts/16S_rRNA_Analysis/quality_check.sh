#This script loads fastqc from hpc.
#Through a loop it runs the fastqc trough the sample and redirects the output to a directory named as fastqc.
module load fastqc/0.11.9 

for file in ./data/*.fq;
do
	fastqc $file -o fastqc1
done
#Loads multiqc from hpc and opens fastqc directory
module load multiqc/1.4 
cd ./fastqc1/
#Runs multiqc 
multiqc ./
