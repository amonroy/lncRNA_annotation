
### this is my working document for my current homolog search project


import sys
import datetime
import itertools
import os
import re


today = datetime.date.today()
window_length = 4
fly = sys.argv[3]

#"python lncRNA_annotation_doc/flybase_flanking_genes.py out/1_dmel_protein_ncRNA_2016-08-22_out.txt flybase/gene_orthologs_fb_2016_03.tsv Dsec flybase/dmel-all-ncRNA-r6.11.fasta flybase/dsec-all-chromosome-r1.3.fasta"
#python lncRNA_annotation_doc/flybase_flanking_genes.py out/1_dmel_protein_ncRNA_2016-08-22_out.txt flybase/gene_orthologs_fb_2016_03.tsv Dsim flybase/dmel-all-ncRNA-r6.11.fasta flybase/dsim-all-chromosome-r2.02.fasta out/dmel-on-dsim-2.blstn out/tblstn-dmel-on-dsim.tblstn flybase/dmel-all-translation-r6.11.fasta
#sys.argv[1] = gff (modified)
#sys.argv[2] = flybase ortholg file, "gene_orthologs_fb_2016_03.tsv"
#sys.argv[3] = fly name (ie Dsec, Dyak, ...)
#sys.argv[4] =  ncRNA fasta file
#sys.argv[5] = scaffold file, species specific
#sys.argv[6] = ncRNA blast format 6
#sys.argv[7] = protein blast format 6
#sys.argv[8] = protein fasta

#protein = sys.argv[6] (?) protein fasta, now sys,argv[8]
#file = sys.argv[4] (ncRNA blast), now sys.argv[6]
#blast_file = sys.argv[5] (protein blast), now sys.argv[7]

def mel_gff_list():
	"""This function takes the modified gff3 file and creates a list"""
	mod_gff3 = sys.argv[1]
	with open(mod_gff3, 'r') as f:
		gff = [line.strip().split('\t') for line in f]
		f.close()
	return gff
	#gff_list ex/:
	#['2L', 'FlyBase', 'gene', '7529', '9484', '.', '+', '.', 'ID=FBgn0031208;Name=CG11023;Ontology_term=SO:0000010,SO:0000087,GO:0016929,GO:0016926;Dbxref=FlyBase:FBan0011023,FlyBase_Annotation_IDs:CG11023,GB_protein:ACZ94128,GB_protein:AAO41164,]
	# ['2L', 'FlyBase', 'ncRNA', '286383', '288292', '.', '+', '.', 'ID=FBtr0347595;Name=CR46263-RA;Parent=FBgn0267996;Dbxref=FlyBase_Annotation_IDs:CR46263-RA;score_text=Weakly Supported;score=0'],

def make_start_dictionary(list):
	"""This function makes a dictionary using the start location as a key and the name as the entry."""
	start_dict = {}
	for i in list:
		#if i[2] == 'protein':
		start = i[3]
		name = i[8].split(';')[0].split('=')[1]
		#print "This is start", start
		#print "this is name", name
		start_dict.setdefault(start,[]).append(name)
		#else:
			#continue
	#print start_dict
	#'14659561': ['FBpp0082426', 'FBpp0082427'], '7438647': ['FBpp0289592'], '13862097': ['FBpp0086712', 'FBpp0308975'], '25477102': ['FBpp0084211'], '17407882': ['FBpp0082878'], 
#'17799003': ['FBpp0082952', 'FBpp0312514'], '11167232': ['FBpp0075947', 'FBpp0075948'], 
#'334286': ['FBpp0297520'], '18554336': ['FBpp0289650', 'FBpp0309333'], '1108701': ['FBpp0077619'], '5139745': ['FBpp0076869'], 
#'7518900': ['FBpp0081076'], '6040702': ['FBpp0293888', 'FBpp0293889'], 
#'12379652': ['FBpp0087067'], '13974233': ['FBpp0080188', 'FBpp0301776', 'FBpp0310237'],
	#quit()
	return start_dict

def mel_ncRNA_strand(list): #2 function
	"""This function takes the gff_list and makes an ncRNA_list"""
	strand_dic = {} #initiates list
	for i in list:
		if i[2] == 'ncRNA':
			#print i
			#quit()
			idRNA = i[8].split(';')[0].split('=')[1]
			strand = i[6]
			start = i[3]
			strand_dic[idRNA] = [strand, start]
			
	#print strand_dic
	return strand_dic
	#'FBtr0344953': '-', 'FBtr0345411': '+', 'FBtr0345412': '-', 'FBtr0334097': '-',

def mel_ncRNA_chrom():
	"""This has to read in the ncRNA fasta file, and find the chromosome name associated with the ncRNA"""
	ncRNA_to_chrom_dict = {}
	ncRNA = sys.argv[4]
	with open(ncRNA, 'r') as n:
		for line in n:
			if line.startswith('>'):
				data = line.strip().split(';')
				#print data
				#quit()
				ncRNA = data[0].split(' ')[0][1:]
				chrom = [data[1].split(':')[0]]
				ncRNA_to_chrom_dict[ncRNA]= chrom
			else:
				continue
		return ncRNA_to_chrom_dict
		#'FBtr0340236': ' loc=2L', 'FBtr0340585': ' loc=2L', 'FBtr0340584': ' loc=2L'


def indexing_location(rna_strand_dic,location_dic, window_length):
	"""This function is going to make a dictionary with upstream genes and downstream genes for each lncRNA"""
	print "I am in indexing"
	#location_dic '14659561': ['FBpp0082426', 'FBpp0082427'], '7438647': ['FBpp0289592']
	fbgn_id_dict = {}
	sorted_starts = sorted(location_dic.keys()) #making an ordered list of the key in location_dic
	#print "length of sorted_starts", len(sorted_starts) # = 18651
	for k,v in rna_strand_dic.iteritems():
		#print "This is k, v[1]:", k, v[1]# FBtr0340236  loc=2L
		#get the index of the lncRNA
		rna_index = sorted_starts.index(v[1]) #finds the index in sorted_starts of rna currently on ... k is whatever, but v[1] is the start location of that k
		up_counter = 1
		down_counter = 1
		uplist, downlist = [], []
		#match = 
		#'FBpp"
		while len(uplist) < window_length: #this is counting how many genes we've found
			#up_counter = up_counter + 1 #this is counting the index
		#for i in range(rna_index - window_length, rna_index): #upstream
			#print location_dic.get(sorted_starts[i])
			try:
				#print location_dic.get(sorted_starts[i])[0] #type = list
				#print location_dic.get(sorted_starts[rna_index - up_counter])
				#['FBpp0312346']
				#quit()
			#v_up = len(location_dic.get(sorted_starts[rna_index - up_counter])) -1  #this is counting the length of the value list in dic
				#print v_up
			#while v_up > 0:
			#	print location_dic.get(sorted_starts[rna_index - up_counter][v_up])
					#['FBpp0312346']
				for i in location_dic.get(sorted_starts[rna_index - up_counter]):
					match = re.search("FBpp*", i)
					if match:
						#print "Match:", i
						uplist.append(i)
						up_counter += 1
						break
					elif not match:
						print "no match:", location_dic.get(i)
						up_counter += 1
						break	
					else:
						print "No, I am here... why????"
						#up_counter += 1
						break
			except IndexError:
				print "This number gave me a problem", rna_index- up_counter
				#gets stuck at 18649, length of sorted starts = 18651
				break
				#up_counter += 1
		#print "This is uplist", uplist
		while len(downlist) < window_length:
			try:
				#print location_dic.get(sorted_starts[rna_index + down_counter])
				for i in location_dic.get(sorted_starts[rna_index - down_counter]):
					match = re.search("FBpp*", i)
					if match:
						#print "Match:", i
						downlist.append(i)
						down_counter += 1
						break
					elif not match:
						#print "no match:", location_dic.get(i)
						down_counter += 1
						break	
					else:
						print "No, I am here... why????"
						#up_counter += 1
						break
			except IndexError:
				print "This number gave me a problem", rna_index- down_counter
				break
		fbgn_id_dict[k]=uplist, downlist
	#print len(fbgn_id_dict)
	#print fbgn_id_dict
	print len(fbgn_id_dict)
	for k,v in fbgn_id_dict.iteritems():
		if len(v[0]) < 4:
			print k , v
			print rna_strand_dic[k][1]
			search = rna_strand_dic[k][1]
			unknown =sorted_starts.index(search)
			print unknown
					
		if len(v[1]) < 4:
			print k, v
			print rna_strand_dic[k][1]
			search = rna_strand_dic[k][1]
			unknown =sorted_starts.index(search)
			print unknown
			#except IndexError:
			#	print "There was an indexerror"
			#	continue
		
		#quit()
		#for i in range( (rna_index + 1), rna_index + (window_length + 1)): #downstream
		#s	try:
				#print location_dic.get(sorted_starts[i])[0]
		#		downlist.append(location_dic.get(sorted_starts[i])[0])
		#	except IndexError:
		#		continue
		#fbgn_id_dict[k]=uplist, downlist
		#'FBtr0343766': (['FBpp0112245', 'FBtr0346835', 'FBpp0082602', 'FBpp0075299'], ['FBpp0075301', 'FBpp0075302', 'FBpp0080297', 'FBpp0271810']), 
		#'FBtr0343760': (['FBpp0083204', 'FBpp0085644', 'FBpp0083205', 'FBtr0345374'], ['FBpp0088212', 'FBpp0080829', 'FBpp0080849', 'FBpp0304540']), 
		#'FBtr0343761': (['FBpp0312349', 'FBpp0088246', 'FBpp0077268', 'FBpp0303136'], ['FBpp0112952', 'FBpp0070495', 'FBpp0070496', 'FBpp0072923']), 
		#'FBtr0343762': (['FBtr0347071', 'FBpp0083564', 'FBpp0083562', 'FBpp0085263'], ['FBpp0303394', 'FBpp0083559', 'FBpp0071659', 'FBpp0070372'])
	return fbgn_id_dict
	#quit()
	
def mel_ncRNA_up_down_dict(rna_chrom_dic, gff_list, window_length):
	"""This function takes our two lists, ncRNA and gff_list and makes a dictionary. fbgn_id_dict where the key is ncRNA and the values are upstream genes and downstream genes. """
	fbgn_id_dict = {}
	for k in rna_chrom_dic.iterkeys():
		#print "This is r in rna_list", r
		for i in gff_list:
			data = i[8].split(';')[0].split('=')[1]
			if data == k:
				#print r, i #
				#post_data = data.split('=')[1] # ex/ FBgn0031208
				#ncRNA_gff_dict[post_data] = i[0], i[3], i[4]
				#[['2L'], 'FlyBase', 'gene', ['7529'], ['9484'], '.', '+', '.', 'ID=FBgn0031208;Name=CG11023;...']
				index = gff_list.index(i)# indexing so more efficient to move backward
				upstream = 0
				counter = 0
				downstream = 0
				anticounter = 0
				ticker=0
				up,down = [],[] #initiating upstream and downstream gene stuff
				while upstream < window_length:
					#print index
					counter = counter + 1
					#this is how we can move backward
					#print gff_list[index-counter]
					if gff_list[index-counter][2] == 'protein':
						#print gff_list[index-counter][2]
						#print gff_list[index-counter]
						#quit()
						#print gff_list[index-counter]
						#print gff_list[index-counter][3]
						#print gff_list[index-counter][4]
						info = gff_list[index-counter][8].split(';')[0].split('=')[1]
						#info_list = [gff_list[index-counter][3], gff_list[index-counter][4]]
						#print "This is info", info
						#print "This is info_list", info_list
						#This is info FBpp0309212
						#prot_dict = {}
						#prot_dict[info] = info_list
						#print prot_dict
						#quit()
						upstream = upstream + 1
						
					else:
					   continue
					   
					up.append(info,)
					
				while downstream < window_length:
					#adding one for each iteration
					anticounter = anticounter + 1
					# the problem is that the last set of ncRNA downstream genes only goes to 3
					try:
						test = gff_list[index+anticounter][2]
						if gff_list[index+anticounter][2] == 'protein':
							info2 = gff_list[index+anticounter][8].split(';')[0].split('=')[1]
							
							downstream = downstream + 1
						else:
							continue
					except IndexError:
						downstream = downstream + 1
						continue
						
					down.append(info2,)
					
				idRNA = k
				#"front-back-gene-id-dict"
				fbgn_id_dict[idRNA] = [up, down]
				#print fbgn_id_dict
				#quit()
	#print fbgn_id_dict
	#'FBtr0343760': [['FBpp0311265', 'FBpp0074577', 'FBpp0074576', 'FBpp0311065'], ['FBpp0074565', 'FBpp0074566', 'FBpp0074575', 'FBpp0074573']], 
	#'FBtr0343761': [['FBpp0305506', 'FBpp0305505', 'FBpp0070493', 'FBpp0070491'], ['FBpp0070495', 'FBpp0070496', 'FBpp0070511', 'FBpp0070494']], 
	return fbgn_id_dict
	#fbgn_id_dict
	#'FBtr0345733': [['FBgn0266879', 'FBgn0266878', 'FBgn0267987', 'FBgn0051973'], ['FBgn0067779', 'FBgn0266322', 'FBgn0031213', 'FBgn0031214']]
	#{'FBtr0336987': (['FBgn0265149', 'FBgn0262252', 'FBgn0031235', 'FBgn0263465'], ['FBgn0022246', 'FBgn0031238', 'FBgn0031239', 'FBgn0265150']), 'FBtr0309810': (['FBgn0263584', 'FBgn0031209', 'FBgn0002121', 'FBgn0031208'], ['FBgn0051973', 'FBgn0267987', 'FBgn0266878', 'FBgn0266879']),}

def prot_to_gn(nc_pp_dic, pp_gn_dic):
	new_dic = dict()
	for k,v in nc_pp_dic.iteritems():
		#print "This is v:", v
		#print "This is v[0]:", v[0]
		#print "This is v[1]:", v[1]
		up_list = []
		down_list = []
		for i in v[0]:
			new_up= pp_gn_dic[i]
			up_list.append(new_up)
		for j in v[1]:
			new_down = pp_gn_dic[j]
			down_list.append(new_down)
		new_dic[k]= [up_list, down_list]
	print new_dic
	#'FBtr0347262': [['FBgn0011766', 'FBgn0011766', 'FBgn0011766', 'FBgn0011766'], ['FBgn0038893', 'FBgn0038894', 'FBgn0259113', 'FBgn0051176']], 
	return new_dic
	

def mel_gene_set(dict): # this uses the flanking genes, specifically
	"""This function finds unique mel genes, and puts them in a set (what is returned), so we don't get the same coords twice. It takes fbgn_id_dict. This is so we have the mel genes that we need coordinates for in the non-mel species', ie we're using this to find the orthologs that we care about"""
	mel_gene_set = set()
	for k, v in dict.iteritems():
		#v[0] is up, v[1] is down
		#print "this is v:", v
		for mg in v[0]:
			mel_gene_set.add(mg)
		for mg in v[1]:
			mel_gene_set.add(mg)
	return mel_gene_set

def map_mel_gene_to_ortho_gene(set):
	"""This function maps other another species' orthologs to dmel genes"""
	mapping = dict()
	with open(sys.argv[2], 'r') as orthos: #this is finding the  ortho_gene_coords
		for line in orthos:
			if not line.startswith('#') and not line.startswith('\n'):
				data = line.split ('\t')
				###switch this to "fly" sys.argv
				if fly in data[6]:
					if data[0] in set:
						coord = data[8].split("..")
						try:
							if 'nonp' in mapping[data[0]]:
								mapping[data[0]][data[5]] = [data[5], data[7], coord[0], coord[1]]
						except KeyError:
							mapping[data[0]]={}
							mapping[data[0]]['nonp']= [data[5], data[7], coord[0], coord[1]]
				
	return mapping
	# mapping ex/
	#'FBgn0085212': {'nonp': ['FBgn0167159', 'scaffold_3', '5320158', '5336608']}, 
	#'FBgn0266101': {'FBgn0180013': ['FBgn0180013', 'scaffold_0', '1697815', '1706253'], 'nonp': ['FBgn0180014', 'scaffold_0', '1707282', '1708727']}	

def ortho_up_down_dict(fbgn_id_dict, map):
	"""This function maps the mel lncRNA to its upstream and downstream orthologs in a different species"""
	rna_ortho = dict()
	for k, v in fbgn_id_dict.iteritems():
		before, after = [], []
		code = []
		for i, gene1 in enumerate(v[0]):
			ortho_gene = map.get(gene1, None) # get returns a value for a given key
			if ortho_gene is None:
				code.append('0')
				continue
			if i <= window_length:
				code.append('1')
				before.append(ortho_gene)
				
		for i, gene2 in enumerate(v[1]):
			ortho_gene = map.get(gene2, None)
			if ortho_gene is None:
				code.append('0')
				continue
			if i <= window_length:
				code.append('1')
				after.append(ortho_gene)
		score =''.join(code)
		rna_ortho[k] = [before,after, score] #from tuple to list
	return rna_ortho
#'FBtr0345411': [[], [{'nonp': ['FBgn0179621', 'scaffold_0', '3565044', '3567184']}], '00000100'], 'FBtr0345412': [[], [{'nonp': ['FBgn0179621', 'scaffold_0', '3565044', '3567184']}, {'nonp': ['FBgn0180128', 'scaffold_0', '3596896', '3597756']}], '00001001'], 'FBtr0334097': [[{'nonp': ['FBgn0179631', 'scaffold_0', '3361223', '3367367']}, {'nonp': ['FBgn0180115', 'scaffold_0', '3376393', '3378580']}],

def ortho_final_coord(ortho_dict):#rna_ortho_dict,
	"""This function finds the end of the front gene ortholog and the front of the back gene ortholog, to give a dictionary with a putative start and putative stop"""
	###### dictionary of dictionaries? when there's two scaffolds that match, does it just overwrite? I think so. Double check this
	final_coord_dict = dict()
	for k, v in ortho_dict.iteritems():
		upstream = v[0]
		downstream = v[1]
		uscafs = set()
		dscafs = set()
		for d in upstream:
			for key, value in d.iteritems():
				uscafs.add(d[key][1]) #scaffold
		for e in downstream:
			for key1, value1 in e.iteritems():
				dscafs.add(e[key1][1])
		common_scaf = uscafs.intersection(dscafs)

		for x in common_scaf:
			upos = []
			for d in upstream:
				for key2, value2 in d.iteritems():
					if d[key2][1] == x:
						upos.append(d[key2][2])
						upos.append(d[key2][3])
			dpos = []
			for e in downstream:
				for key3, value3 in e.iteritems():
					if e[key3][1] == x:
						dpos.append(e[key3][2])
						dpos.append(e[key3][3])
			upos_num = [int(n) for n in upos]
			dpos_num = [int(n) for n in dpos]
			merged_pos = upos_num + dpos_num
			try:
				if x in final_coord_dict[k]:
					print "It looks like this scaffold showed up twice", x
				#final_coord_dict[k] = [x, min(merged_pos), max(merged_pos)]
			except KeyError:
				final_coord_dict[k] = {}
				final_coord_dict[k][x] = [x, min(merged_pos), max(merged_pos)]
		else:
			continue
	#print final_coord_dict
	return final_coord_dict
	#{'FBtr0344114': {'scaffold_39': ['scaffold_39', 318497, 337630]}, 'FBtr0345732': {'scaffold_14': ['scaffold_14', 31600, 75558]}, , 'FBtr0344645': {'scaffold_0': ['scaffold_0', 3565044, 3628549]}, 
	
def ortho_read_scaffolds():
	"""This reads in non-mel-species chromosome (err scaffold) file to make a dictionary of the key = name of scaffold and values of the scaffold sequence"""
	#user_input = raw_input("Enter the path of the scaffold file: ")
	#assert os.path.exists(user_input), "I did not find the file at, "+str(user_input)
	#print ("Hooray we found your file!")
	#print "I am in karl_read"
	non_mel_scaff = sys.argv[5]
	with open(non_mel_scaff, 'r+') as species:
		key = ''
		sequence = []
		scaffold = dict()
		for line in species:
			if line.startswith('>'):
				if key and sequence:
					scaffold[key] = ''.join(sequence)
				key = line.split(' ')[0][1:]
				sequence = []
			else:
				sequence.append(line.rstrip())
		scaffold[key] = ''.join(sequence)
	return scaffold

def mel_ncrna_seq_dict():
	"""This reads in flybase mel ncRNA file and makes a dictionary of ncRNA id (key) and sequence (value)"""
	RNA_file = sys.argv[4]
	with open(RNA_file, 'r+') as file:
		key = ''
		sequence = []
		ncRNA = dict()
		for line in file:
			if line.startswith('>'):
				if key and sequence:
					ncRNA[key] = ''.join(sequence)
				#print "This is line"
				#print line
				#>FBtr0076899 type=ncRNA; loc=3L:join(6976322..6977006,6977169..6977479,6977647..6979222); ID=FBtr0076899; name=CR32385-RA; dbxref=FlyBase:FBtr0076899,FlyBase_Annotation_IDs:CR32385-RA,REFSEQ:NR_001944,REFSEQ:NR_001944; score=9; score_text=Strongly Supported; MD5=25c515079ad37d428daae3062ba15e0c; length=2572; parent=FBgn0047205; release=r6.11; species=Dmel; 
				key = line.split(' ')[0][1:]
				sequence = []
			else:
				sequence.append(line.rstrip())
		ncRNA[key] = ''.join(sequence)
	return ncRNA

def reverse_complement(strng):
	"""This function returns the reverse complement of a string of nucleotides"""
	#print "This is your starting string:"
	#print type(strng)
	reverse = strng[::-1]
	#print "This is the reverse of your string"
	#print reverse
	complement = []
	for i in reverse:
		if i == 'a' or i == 'A':
			complement.append('t')
		if i == 't' or i == 'T':
			complement.append('a')
		if i == 'c' or i == 'C' :
			complement.append('g')
		if i == 'g' or i == 'G':
			complement.append('c')
	#print "This is the reverse complement"
	#print complement
	rev_com_str = ''.join(complement)
	#print "This is the string rev com"
	#print rev_com_str
	return rev_com_str

def merge_seq_and_strand(mel_seq_dic, strand_dic):
	"""This function adds strandedness to important mel_seq_info"""
	mel_complete = {}
	for k,v in mel_seq_dic.iteritems():
		try: ## remove this try/ except!!!!
			strand = strand_dic[k]
			mel_complete[k] = v, strand
			#print mel_complete[k]
			#print mel_complete[k][1]
		except KeyError:
			continue
	#quit()
	return mel_complete
	#'FBtr0344953': ('ACAAATATACTTTGAAATCTTTTGGTCGCCAAATC...', '-')
	
def mel_ortho_output_sequence(scaf_dict, coord_dict, mel_seqs, mel_ud_dict):
	"""This print the sequence of putative lncRNA in non-mel species, using start and stop and scaffold"""
	#make out file to be written to
	out_file = open('out_flybase_%s_ncRNA_mel_and_ortho_seq_%s.txt' %(fly, today), 'w')
	for k, v in coord_dict.iteritems():
		#This is k: FBtr0344113
		#This is v: {'scaffold_39': ['scaffold_39', 290730, 337630]}
		#print "This is v:", v
		#print "This is v[1]", v[1]
		for l,w in v.iteritems():
			#print "This is l", l
			#print "This is w", w
			if w[1] < w[2]:
				start = int(w[1])+1
				end = int(w[2])+1
			if w[2] < w[1]:
				start = int(w[2])+1
				end = int(w[1])+1
		if mel_seqs[k][1] == '+':
			out_file.write(">%s;mel;%s\n%s\n" %(k, '+', mel_seqs[k][0]))
		if mel_seqs[k][1] == '-':
			r_string = mel_seqs[k][0]
			rev_com = reverse_complement(r_string)
			#print "This is rev_com", rev_com
			out_file.write(">%s;dmel;%s\n%s\n" %(k, '-', rev_com))
		out_file.write(">%s;%s;%s;%s\n%s\n" %(k, l, fly, mel_ud_dict[k][2], scaf_dict[l][start:end]))
	out_file.close()
	
def blast_coord(coord_dict):
	"""This function will find all of the hits of dmel ncRNA against non-mel genome and return the best scaffold and the start and stop"""
	file = sys.argv[6]
	blast_coord_dict = dict()
	final_blast_dict = dict()
	with open(file, 'r') as f:
		for line in f:
			data= line.strip().split('\t')
			if data[0] in blast_coord_dict:
				blast_coord_dict[data[0]].append(data[1:])
			else:
				blast_coord_dict[data[0]]= [data[1:]]
	
	for k,v in blast_coord_dict.iteritems():
		if k in coord_dict:
			sstart = []
			send = []
			for list in v:
				uscafs =set()
				uscafs.add(list[0])
				sstart.append(int(list[7]))
				send.append(int(list[8]))
			#here, since we are finding it compared to the lncRNA, we want more buffer, so the lowest and highest between up and down
			scaf= [i for i in uscafs]
			final_blast_dict[k] = scaf[0], min(sstart), max(send)
	return blast_coord_dict, final_blast_dict


def find_best_hit():
	"""This reads in your protein blast output (in format 6) and gives you the data for the best hit per protein. Makes a dictionary."""
	#print "In find_best_hit()"
	blast_file = sys.argv[7]
	protein_dict = {}
	with open(blast_file, 'rU') as f:
		for line in f:
			data = line.strip().split('\t')
			key = data[0]
			if key in protein_dict:
				if float(protein_dict[key][3]) < float(data[10]):
					continue
				if float(protein_dict[key][3])> float(data[10]):
					protein_dict[key]= [data[1], data[8], data[9], data[10]]
				if float(protein_dict[key][3]) == float(data[10]):
					continue
			else:
				protein_dict[key]= [data[1], data[8], data[9], data[10]]
	return protein_dict
	#'FBpp0297216': ['scaffold_1', '11744116', '11743913', '1.05e-16'], 
	#'FBpp0297217': ['scaffold_9', '2334999', '2335109', '5.05e-26'], 
	#'FBpp0297219': ['scaffold_25', '588672', '588544', '3.43e-25']
def match_pp_to_gn():
	"""This has to read in the stupid protein files, and find the gn name associated with the protein"""
	pp_to_gn_dict = {}
	protein = sys.argv[8]
	with open(protein, 'r') as p:
		for line in p:
			if line.startswith('>'):
				data = line.strip().split(';')
				pp = data[2].split('=')[1]
				gn = data[4].split(',')[0].split('=')[1]
				pp_to_gn_dict[pp]=gn
			
			else:
				continue
		return pp_to_gn_dict
		#'FBpp0088107': 'FBgn0033139', 
		#'FBpp0088104': 'FBgn0016032', 
		#'FBpp0088105': 'FBgn0033135', 
		#'FBpp0298338': 'FBgn0263113', 

def dmelgn_to_ortho(p2g_dict, best_hit_dict): # 11_fun
	"""This function takes pp_to_gn_dict and protein_dict and makes dmelgn_to_ortho_dict. It links the dmel gn to the best hit."""
	dmelgn_to_ortho = {}
	for k,v in best_hit_dict.iteritems():
		#print "this is k", k
		new_key = p2g_dict[k]
		#print "This is new_key", new_key
		dmelgn_to_ortho[new_key] = {}
		dmelgn_to_ortho[new_key]['best_hit']= new_key, best_hit_dict[k][0],best_hit_dict[k][1], best_hit_dict[k][2]

	return dmelgn_to_ortho# needs to be a dict of dict
	#FBgn0032006': {'best_hit': ('FBgn0032006', 'scaffold_3', '3723783', '3722836')}

def trouble_shoot(mel_ncRNA_dic, mel_ud_gn_dic, ortho_ud_gn_dic, ortho_final_coord, mel_prot_blast, mel_rna_blast, mel_ud_pp_dic):
#(mel_ncRNA_chrom_obj,mel_ud_gn_dict_obj, ortho_ud_gn_dict, ortho_final_coord_obj, mel_to_ortho_prot_map, final_bls_dict)
#prot_blast, rna_blast, obj1, obj2, 
	"""This outputs a table that gives relevant troubleshooting info for each key (ncRNA gene)"""
	trouble_out = open('out_troubleshoot_%s_ncRNA_flanking_genes_%s.txt' %(fly, today), 'w')
	trouble_out.write("#ncRNA_id\t#ortho_found\t#mel_chrom\t#up_proteins\t#down_proteins\t#up_genes\t#down_genes\t#ortho_gene_score\t#flybase_ortho_scaffold\t#prot_ortho_scaffold\t#rna_ortho_scaffold\n")
                            #k,         seqs, mel_chrom,    pp_up,     pp_down, p_up, p_down, score,        fb_ortho_scaffs,       prot_ortho_scaff,       rna_ortho_scaff)
##ncRNA_id	##ortho_found	#mel_chrom	#up_proteins	#down_proteins	#up_genes	#down_genes	#ortho_gene_score	#flybase_ortho_scaffold	#prot_ortho_scaffold	#rna_ortho_scaffold
	for k,v in mel_ncRNA_dic.iteritems():

		#out_file.write("%s\tncRNAbls:%s\tortho:%s\tprotbls:%s\n" %(k, final_bls_dict[k][0], ortho_final_coord_obj[k][0], new_ortho_final_coord_obj[k][0]))
		mel_chrom = mel_ncRNA_dic[k]
		if k in ortho_final_coord:
			seqs= 'yes'
		else:
			seqs= 'no'
		try:
			pup = mel_ud_pp_dic[k][0]
			pp_up = ';'.join(pup)
		except KeyError:
			pp_up = "KeyError"
		#print "pp_up:", pp_up
		try:
			pdown = mel_ud_pp_dic[k][1]
			pp_down = ';'.join(pdown)
		except KeyError:
			pp_down = "KeyError"
		try:
			up = mel_ud_gn_dic[k][0]
			p_up = ';'.join(up)
		except KeyError:
			p_up = "KeyError"
		try:
			down = mel_ud_gn_dic[k][1]
			p_down = ';'.join(down)
		except KeyError:
			p_down = "KeyError"
		try:
			score = ortho_ud_gn_dic[k][2]
		except KeyError:
			score = "KeyError"
		fb_ortho_scaffs = []
		try:
			for v1 in ortho_final_coord[k]:
				#print "This is v1", v1
				fb_ortho_scaffs.append(v1)
		except KeyError:
			fb_ortho_scaffs.append("KeyError")
		p_fb_ortho_scaffs = ';'.join(fb_ortho_scaffs)
		try:		
			
			if mel_prot_blast[k]:
				key_list=[]
				for kay in mel_prot_blast[k].iterkeys():
					#print kay
					key_list.append(kay)
					
				prot_ortho_scaff = ';'.join(key_list)
		except KeyError:
			#'FBtr0343761': {'scaffold_10': ['scaffold_10', 2846888, 2975383]}
			#print mel_prot_blast[k]
			prot_ortho_scaff= "KeyError"
		try:
			rna_ortho_scaff = mel_rna_blast[k][0]
		except KeyError:
			rna_ortho_scaff= "KeyError"
		trouble_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(k, seqs, mel_chrom, pp_up, pp_down, p_up, p_down, score, fb_ortho_scaffs, prot_ortho_scaff, rna_ortho_scaff))
		#print k, mel_chrom, p_up, p_down, score, fb_ortho_scaffs, prot_ortho_scaff, rna_ortho_scaff
		#       %s  %s  %s    %s  %s    %s     %s                  %s              %s
		
	trouble_out.close()
	
		
		
			
		
		
		



### making sure the functions are executed in the right order
mel_gff_obj = mel_gff_list()

mel_loc_dictionary = make_start_dictionary(mel_gff_obj)

mel_strand =mel_ncRNA_strand(mel_gff_obj)

mel_ncRNA_chrom_obj = mel_ncRNA_chrom()

mel_ud_pp_dict_obj = indexing_location(mel_strand,mel_loc_dictionary, window_length)

#mel_ud_pp_dict_obj = mel_ncRNA_up_down_dict(mel_ncRNA_chrom_obj, mel_gff_obj, window_length)

prot_to_gene = match_pp_to_gn()

mel_ud_gn_dict_obj = prot_to_gn(mel_ud_pp_dict_obj, prot_to_gene)

##need an intermediate here where I replace the protein with a gn

mel_genes_obj = mel_gene_set(mel_ud_gn_dict_obj)

ortho_map = map_mel_gene_to_ortho_gene(mel_genes_obj)

ortho_ud_gn_dict = ortho_up_down_dict(mel_ud_gn_dict_obj, ortho_map)

ortho_final_coord_obj = ortho_final_coord(ortho_ud_gn_dict)

ortho_scaff_dict = ortho_read_scaffolds()

mel_rna_seqs_only = mel_ncrna_seq_dict()

mel_rna_seq_complete = merge_seq_and_strand(mel_rna_seqs_only, mel_strand)

mel_ortho_output_sequence(ortho_scaff_dict, ortho_final_coord_obj, mel_rna_seq_complete, ortho_ud_gn_dict )


bls_dict, final_bls_dict =blast_coord(ortho_final_coord_obj) ##p1##  #8_fun(#7_obj)
#8a_obj, #8b_obj
pp_dict = find_best_hit() #9_fun
#9_obj
#print pp_dict
 #10_fun
#10_obj
#'FBpp0290552': 'FBgn0032006', 'FBpp0290553': 'FBgn0032006', 'FBpp0290550': 'FBgn0039932',
mel_to_ortho_prot_map =dmelgn_to_ortho(prot_to_gene, pp_dict) #11_fun(#10_obj, #9_obj)

new_ortho_ud_gn_dict = ortho_up_down_dict(mel_ud_gn_dict_obj, mel_to_ortho_prot_map)

new_ortho_final_coord_obj = ortho_final_coord(new_ortho_ud_gn_dict) #7_fun(#6.2_obj)

trouble_shoot(mel_ncRNA_chrom_obj,mel_ud_gn_dict_obj, ortho_ud_gn_dict, ortho_final_coord_obj, new_ortho_final_coord_obj, final_bls_dict, mel_ud_pp_dict_obj)
