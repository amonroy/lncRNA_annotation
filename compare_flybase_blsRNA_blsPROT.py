import sys
import datetime
import itertools
import os

today = datetime.date.today()
window_length = 4
fly = sys.argv[3]
#sys.argv[1] mel_gff_list
#sys.argv[2] ortho_file
#sys.argv[4] ncRNA blast file
#sys.argv[5] protein blast file
#sys.argv[6] protein fasta file
#sys.argv[7] ncRNA fasta file

def mel_gff_list(): #1 function
	"""This function takes the modified gff3 file and creates a list"""
	mod_gff3 = sys.argv[1]
	with open(mod_gff3, 'r') as f:
		gff = [line.strip().split('\t') for line in f]
		f.close()
	return gff
	#gff_list ex/:
	#['2L', 'FlyBase', 'gene', '7529', '9484', '.', '+', '.', 'ID=FBgn0031208;Name=CG11023;Ontology_term=SO:0000010,SO:0000087,GO:0016929,GO:0016926;Dbxref=FlyBase:FBan0011023...'],
	# ['2L', 'FlyBase', 'ncRNA', '286383', '288292', '.', '+', '.', 'ID=FBtr0347595;Name=CR46263-RA;Parent=FBgn0267996;Dbxref=FlyBase_Annotation_IDs:CR46263-RA;score_text=Weakly Supported;score=0'], 

def mel_ncRNA_list(list): #2 function
	"""This function takes the gff_list and makes an ncRNA_list"""
	ncRNA = [] #initiates list
	for i in list:
		if i[2] == 'ncRNA':
			preidRNA = i[8].split(';')[0]
			ncRNA.append(preidRNA)
	return ncRNA
	#['ID=FBtr0309810', 'ID=FBtr0347585', 'ID=FBtr0345732', 'ID=FBtr0345733', 'ID=FBtr0344052']

def mel_ncRNA_up_down_dict(rna_list, gff_list, window_length): #3 function
	"""This function takes our two lists, ncRNA and gff_list and makes a dictionary. mel_ncRNA_up_down_dict where the key is ncRNA and the values are upstream genes and downstream genes. """
	up_down_dict = {}
	for r in rna_list:
		for i in gff_list:
			data = i[8].split(';')[0]
			if data == r:
				index = gff_list.index(i)# indexing so more efficient to move backward
				upstream = 0
				counter = 0
				downstream = 0
				anticounter = 0
				ticker=0
				up,down = [],[] #initiating upstream and downstream gene stuff
				while upstream < window_length:
					counter = counter + 1
					#this is how we can move backward
					if gff_list[index-counter][2] == 'gene':
						pre_info = gff_list[index-counter][8].split(';')[0]
						info = pre_info.split('=')[1]
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
						if gff_list[index+anticounter][2] == 'gene':
							pre_info2 = gff_list[index+anticounter][8].split(';')[0]
							info2 = pre_info2.split('=')[1]
							downstream = downstream + 1
						else:
							continue
					except IndexError:
						downstream = downstream + 1
						continue
						
					down.append(info2,)
					
				idRNA = r.split('=')[1]
				up_down_dict[idRNA] = [up, down] #changed from tuple to list?
	return up_down_dict
	#up_down_dict
	#'FBtr0345055': [['FBgn0266688', 'FBgn0267475', 'FBgn0264001', 'FBgn0036377'], ['FBgn0266595', 'FBgn0262415', 'FBgn0262813', 'FBgn0036380']],
	
def mel_gene_set(dict): # this uses the flanking genes, specifically
	"""This function finds unique mel genes, and puts them in a set (what is returned), so we don't get the same coords twice. It takes up_down_dict. This is so we have the mel genes that we need coordinates for in the non-mel species', ie we're using this to find the orthologs that we care about"""
	mel_gene_set = set()
	for k, v in dict.iteritems():
		for mg in v[0]:
			mel_gene_set.add(mg)
		for mg in v[1]:
			mel_gene_set.add(mg)
	return mel_gene_set

def map_mel_gene_to_ortho_gene(set): #5_fun
	"""This function maps other another species' orthologs to dmel genes"""
	mapping = dict()
	with open(sys.argv[2], 'r') as orthos: #this is finding the  ortho_gene_coords
		not_count = 0
		yes_count = 0
		for line in orthos:
			if not line.startswith('#') and not line.startswith('\n'):
				data = line.split ('\t')
				if fly in data[6]:
					if data[0] in set:
						coord = data[8].split("..")
						try:
							if 'nonp' in mapping[data[0]]:
								mapping[data[0]][data[5]] = [data[5], data[7], coord[0], coord[1]]
						except KeyError:
							mapping[data[0]] = {}
							mapping[data[0]]['nonp']=[data[5], data[7], coord[0], coord[1]]
						yes_count = yes_count + 1
					if data[0] not in set:
						not_count = not_count + 1  
	print "This is yes_count", yes_count
	print "This is not_count", not_count
	return mapping #5_obj
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

def ortho_final_coord(ortho_dict):#rna_ortho_dict,
	"""This function finds the end of the front gene ortholog and the front of the back gene ortholog, to give a dictionary with a putative start and putative stop"""
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
			final_coord_dict[k] = [x, min(merged_pos), max(merged_pos)]
		else:
			continue
			
	return final_coord_dict

def blast_coord(coord_dict):
	"""This function will find all of the hits of dmel ncRNA against non-mel genome and return the best scaffold and the start and stop"""
	file = sys.argv[4]
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

#this works so far
def find_best_hit():
	"""This reads in your protein blast output (in format 6) and gives you the data for the best hit per protein. Makes a dictionary."""
	print "In find_best_hit()"
	blast_file = sys.argv[5]
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
	protein = sys.argv[6]
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
	
def mel_ncRNA_chrom():
	"""This has to read in the ncRNA fasta file, and find the chromosome name associated with the ncRNA"""
	ncRNA_to_chrom_dict = {}
	ncRNA = sys.argv[7]
	with open(ncRNA, 'r') as n:
		for line in n:
			if line.startswith('>'):
				data = line.strip().split(';')
				ncRNA = data[0].split(' ')[0][1:]
				chrom = data[1].split(':')[0]
				ncRNA_to_chrom_dict[ncRNA]= chrom
			else:
				continue
		return ncRNA_to_chrom_dict


mel_gff_obj = mel_gff_list() #1_fun
#1_obj
mel_ncRNA_obj = mel_ncRNA_list(mel_gff_obj) #2_fun(#1_obj)
#2_obj
mel_ud_gn_dict_obj =mel_ncRNA_up_down_dict(mel_ncRNA_obj, mel_gff_obj, window_length)#3_fun(#2_obj, #1_obj, window_length)
#3_obj
mel_genes_obj = mel_gene_set(mel_ud_gn_dict_obj)#4_fun(#3_obj)
#4_obj
#set(['FBgn0027567', 'FBgn0027564', 'FBgn0027565', 'FBgn0027569', 'FBgn0266044', 'FBgn0266045', 'FBgn0266046', ...])
ortho_map = map_mel_gene_to_ortho_gene(mel_genes_obj) #5_fun(#4_obj)
#5_obj
#'FBgn0085212': {'nonp': ['FBgn0167159', 'scaffold_3', '5320158', '5336608']}, 
#'FBgn0266101': {'FBgn0180013': ['FBgn0180013', 'scaffold_0', '1697815', '1706253'], 'nonp': ['FBgn0180014', 'scaffold_0', '1707282', '1708727']}, 
ortho_ud_gn_dict = ortho_up_down_dict(mel_ud_gn_dict_obj, ortho_map) #6_fun(#3_obj,#5_obj)
#6.1_obj
ortho_final_coord_obj = ortho_final_coord(ortho_ud_gn_dict)#7_fun(#6_obj)
#7_obj
bls_dict, final_bls_dict =blast_coord(ortho_final_coord_obj) ##p1##  #8_fun(#7_obj)
#8a_obj, #8b_obj
pp_dict = find_best_hit() #9_fun
#9_obj
#print pp_dict
prot_to_gene = match_pp_to_gn() #10_fun
#10_obj
#'FBpp0290552': 'FBgn0032006', 'FBpp0290553': 'FBgn0032006', 'FBpp0290550': 'FBgn0039932',
mel_to_ortho_map =dmelgn_to_ortho(prot_to_gene, pp_dict) #11_fun(#10_obj, #9_obj)
#11_obj
#print "This is mel_to_ortho_map:"
#print mel_to_ortho_map
#quit()
new_ortho_ud_gn_dict = ortho_up_down_dict(mel_ud_gn_dict_obj, mel_to_ortho_map) #6_fun(#3_obj, #11_obj)
#6.2_obj
#print "new_ortho_ud_gn_dict"
#print new_ortho_ud_gn_dict
#'FBtr0340588': [[{'best_hit': ('FBgn0031514', 'scaffold_5', '1419855', '1418905')}, {'best_hit': ('FBgn0031513', 'scaffold_5', '1412344', '1413807')}], [{'best_hit': ('FBgn0031515', 'scaffold_5', '1435284', '1434919')}, {'best_hit': ('FBgn0031516', 'scaffold_5', '1438869', '1437721')}, {'best_hit': ('FBgn0031518', 'scaffold_5', '1458413', '1459063')}, {'best_hit': ('FBgn0031517', 'scaffold_5', '1453682', '1454188')}], '01101111'], 'FBtr0340585': [[{'best_hit': ('FBgn0085200', 'scaffold_7', '2005863', '2006237')}], [{'best_hit': ('FBgn0261597', 'scaffold_7', '2064366', '2064076')}, {'best_hit': ('FBgn0086707', 'scaffold_7', '2068026', '2067118')}, {'best_hit': ('FBgn0032679', 'scaffold_7', '2075296', '2073395')}, {'best_hit': ('FBgn0032680', 'scaffold_7', '2076161', '2076550')}], '00011111'], 'FBtr0340584': [[{'best_hit': ('FBgn0085200', 'scaffold_7', '2005863', '2006237')}, {'best_hit': ('FBgn0011559', 'scaffold_7', '1985427', '1984948')}], [{'best_hit': ('FBgn0261597', 'scaffold_7', '2064366', '2064076')}, {'best_hit': ('FBgn0086707', 'scaffold_7', '2068026', '2067118')}, {'best_hit': ('FBgn0032679', 'scaffold_7', '2075296', '2073395')}], '00110111'], 'FBtr0340586': [[], [{'best_hit': ('FBgn0031513', 'scaffold_5', '1412344', '1413807')}, {'best_hit': ('FBgn0031514', 'scaffold_5', '1419855', '1418905')}, {'best_hit': ('FBgn0031515', 'scaffold_5', '1435284', '1434919')}], '00001101']
###observation: #11_obj and #5_obj have to have the same data structure ###

new_ortho_final_coord_obj = ortho_final_coord(new_ortho_ud_gn_dict) #7_fun(#6.2_obj)
ncRNA_chrom =mel_ncRNA_chrom() #12_fun
#12_obj




## 8a and 8b are used mostly here
out_file = open('out_compare_%s_ncRNAbls_flybase_proteinbls_%s.txt' %(fly, today), 'w')
ortho_count = 0
bls_ortho_match = 0
keys_tried = 0
key_errors = 0
bad_key_list = []
for k,v in ortho_final_coord_obj.iteritems():#change this to lncRNA genes
	try:
		keys_tried = keys_tried + 1
		out_file.write("%s\tncRNAbls:%s\tortho:%s\tprotbls:%s\n" %(k, final_bls_dict[k][0], ortho_final_coord_obj[k][0], new_ortho_final_coord_obj[k][0]))

		if ortho_final_coord_obj[k][0] == new_ortho_final_coord_obj[k][0]:
			ortho_count = ortho_count + 1
		if ortho_final_coord_obj[k][0] == new_ortho_final_coord_obj[k][0] and ortho_final_coord_obj[k][0] == final_bls_dict[k][0]:
			bls_ortho_match = bls_ortho_match + 1
	except KeyError:
		key_errors = key_errors + 1
		bad_key_list.append(k)
		out_file.write("%s\tKEY_ERROR\n" %k)
out_file.close()

bad_key_rna = 0
bad_key_ortho = 0

for x in bad_key_list:
	try:
		final_bls_dict[x]
	except KeyError:
		bad_key_rna = bad_key_rna+1
	try:
		new_ortho_final_coord_obj[x]
	except KeyError:
		bad_key_ortho = bad_key_ortho +1	

print fly
print "keys_tried", keys_tried
print "bls_ortho_match", bls_ortho_match
print "ortho_count", ortho_count
print "key_errors", key_errors
print "bad_key_ortho", bad_key_ortho
print "bad_key_rna", bad_key_rna

for j in mel_ncRNA_obj:
	bad = False
	i =j.split('=')[1]
	if i not in ortho_final_coord_obj:
		bad = True
		print "%s did not have a common scaffold" %i
	if i not in bls_dict:
		bad = True
		print "%s did not have a ncRNA blast hit" %i
	if i not in final_bls_dict:
		bad = True
		print "%s did not have a ncRNA blast common scaffold" %i
	if i not in new_ortho_final_coord_obj:
		bad = True
		print "%s did not have a protein blast ortholog" %i
	if bad:
		print "ortholog_score:", ortho_ud_gn_dict[i][2]
		print "location:", ncRNA_chrom[i]
