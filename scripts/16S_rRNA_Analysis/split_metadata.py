#!/bin/env/python

#This script splits the metadata file into five metadata files.

with open("../metadata/BSF-Metadata-File.tsv","r") as myfile:
	lines=myfile.readlines()
	print(lines)
	names=[]
	
