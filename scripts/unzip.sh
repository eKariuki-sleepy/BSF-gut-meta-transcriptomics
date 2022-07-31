#!/bin/bash
#This script takes in a zipped file and unzips it and redirects the output to a sortmerna_db directory
for i in sortmerna_db/*.zip;
do
	unzip $i -d sortmerna_db
done

# This command moves the rRNA dtabases into the sortmerna directory

mv sortmerna_db/sortmerna-2.1b/rRNA_databases/ sortmerna_db

#this code removes the unnecessary files and folders

rm -r ./sortmerna_db/sortmerna-2.1b
rm ./sortmerna_db/2.1b.zip
