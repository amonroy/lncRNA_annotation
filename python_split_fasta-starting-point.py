#!/usr/bin/python

##A Monroy
## this is to loop through fasta files and split them into new files with names based on 
## their fasta name


import sys
import re

file_to_split = sys.argv[1]


with open(file_to_split, 'r') as f:
	while f:
		currentline = f.readline()
		if currentline.startswith('>'):
			data = currentline.strip()
			header = re.findall('>(\S+)\s', data)
			fasta_name = header[0]
			print "This is fasta_name: %s" %fasta_name
			out = open (fasta_name + ".fa", 'w')
			out.write("%s\n" %data)
			currentline = f.readline
		if not currentline.startswith('>'):
			data = currentline.strip()
			out.write("%s" %data)
			currentline = f.readline
		if not currentline:
			break
		
	#currentline = f.readline
	
	
			
			
out.close()




#####seems to do everything as I expect except will not print last fasta???????f


