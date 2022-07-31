#!/bin/env/python3
import os
import glob
#this script takes in the Fastq files, runs blast on them and then writes the outputs in the txt files
for file in glob.glob('*.fq'):
    new = file.rsplit('.',1)[0] + '.txt'
    filename = '/opt/data/oscarmwaura/transcriptome/bac-16S/blast_dir/16S_ribosomal_RNA'

    blast = "blastn " + "-db " + filename  + " -query " + file + " -out " + new
    os.system (blast)
#test-18S.py (END)
