#!/bin/env/python3

import os
import glob

#This script takes in the blast output and extracts the accession numbers


for file in glob.glob ("*.txt"):
	name = file.split("_",2)[1] + "_accesion" + ".txt"
	with open(file,'r') as myfile:
		with open(name,"w") as outfile:
			outfile.write('These are the accession numbers \n')
			#outfile.write(field)	
			lines = myfile.readlines()
			for i in lines:
				if i.startswith('NR'):
					lines = i.split()
					field = lines[0]
					join = "".join(field)
					outfile.write(join + "\n")
		print(name)
