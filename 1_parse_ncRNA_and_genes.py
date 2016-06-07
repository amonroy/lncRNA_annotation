import sys

gff3 = sys.argv[1]

ncRNA_gene_out = open('1_ncRNA_and_genes_out-2.txt', 'w')

f= open(gff3, 'r')

while True:
	line = f.readline()
	if not line.startswith('#'):
		data = line.strip().split("\t")
	
		if True == True:
	
			try:
		
				if data[2] == 'gene':
					ncRNA_gene_out.write('%s' %line)
			
				elif data[2] == 'ncRNA':
					ncRNA_gene_out.write('%s' %line)
				
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

ncRNA_gene_out.close()			
f.close()			
