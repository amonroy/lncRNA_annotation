###This is the start, rewriting (yet again, starting off of step_by_step_flanking)


### starting materials:
##sys.argv[1] is the gff file
##for practice I have been using..."1_prac_ncRNA_and_genes.txt"

import sys
#import datetime
import itertools

## starting variables to be used later:
#today = datetime.date.today()
window_length = 4 # assigns how many upstream genes and (the same number of) downstream genes
fly_types_list = ['Dsim', 'Dsec', 'Dyak', 'Dere', 'Dana', 'Dpse', 'Dper', 'Dwil', 'Dmoj', 'Dvir', 'Dgri']

#this works, now can add the actual gff file, not the modified one
def parse_gff_nc_genes():
	"""This function takes the Dmel gff file from argv[1] and parses it. 
	It returns lines with 'ncRNA' and 'genes' in them. And removes the sequences at the end of the file.
	It uses sys and datetime"""
	
	gff3 = sys.argv[1]

	#today = datetime.date.today()

	#ncRNA_gene_out = open('1_dmel_genes_ncRNA_%s_out.txt' %today, 'w')

	f= open(gff3, 'r')
	gff = [] # initiates gff list
	mel_ncRNA_coords = {}
	while True:
		line = f.readline()
		if not line.startswith('#'):
			data = line.strip().split("\t")
			if True == True:
	
				try:
		
					if data[2] == 'gene':
						#print "This is gene data", data
						gff.append(data,)
						#print "This is gff", gff
						#ncRNA_gene_out.write('%s' %line)
			
					elif data[2] == 'ncRNA':
						#ncRNA_gene_out.write('%s' %line)
						#print "This is ncRNA data", data
						gff.append(data,)
						post_data = data[8].split(';')[0]
						mel_ncRNA_coords[post_data] = data[0], data[3], data[4], data[6]
						#print data
						#print mel_ncRNA_coords_dict
						#quit()
				
					else:
						continue
				
				except IndexError:
					if line.startswith('>'):
						break
					elif line.startswith('A'):
						break
					elif line.startswith('T'):
						break
					elif line.startswith('G'):
						break
					elif line.startswith('C'):
						break
				
			if not line:
				break

	#ncRNA_gene_out.close()			
	f.close()
	#print "This is gff:", gff
	return gff, mel_ncRNA_coords			

gff_list, mel_coords_dict = parse_gff_nc_genes()
print gff_list
print mel_coords_dict


