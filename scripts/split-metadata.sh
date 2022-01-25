#!/bin/env/python

# this script reads in the metadata file and removes the subject column 
# which is the last column and writes a new metadata file without it 

with open("../Data/BSF-Metadata-File.tsv","r") as myfile:
    with open("../Data/new_metadata.tsv","w") as outfile:
        lines = myfile.readlines()
        for line in lines:
            line_split = line.split()
            fields = (line_split[0:4])
            joined = '\t'.join(fields)
            outfile.write(joined + '\n')
            
# This script takes in the metadata
# Creates the tsv files from the sample IDs
# Writes the metadata of each sample in the tsv files located in the Data folder

data = "../Data/new_metadata.tsv"

with open(data,"r") as file:
    for line in file:
        if line.startswith('#'):
            header = line
        else:
            
            with open("../Data/" + line.split()[0]+'.tsv', 'w') as write_file:
                write_file.write(header)
                write_file.write(line)
           
