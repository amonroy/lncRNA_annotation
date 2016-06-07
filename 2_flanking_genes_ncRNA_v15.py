### TODO today (todays notes to myself)
#
#work on quality control

### TODO (notes to myself)
# this is not finished! stuck as to what to do with ncRNA and genes after; what's my output; have a basic idea
# need to do something to keep track of "quality"
# need to make sure all information we need is accessible

### version: "2_flanking_genes_ncRNA_v14.py"

### PURPOSE:
# The purpose of this script is to take known lncRNA "regions" in Dmel and find homologs of Dmel lncRNAs in non-mel flies.
# The way I am doing this is by finding [window length] upstream protein coding "genes" and [window length] downstream coding "genes".
# I am doing this, because the genes' homologs have been found across non dmel species, and all compare back to dmel, including gene coordinates.
# This is the best way I can think of narrowing down non-annotated, non-coding DNA in non Dmel species

# The goal of this script is to have at the end:
#      Dmel lcnRNA name : ncgene region in nonmel species (start, stop)
#      *the start on stop are based on the last coordinate of the upstream gene, and the first coordinate of the downstream gene

### PARAMETERS:
# You can input different window_length, ie how many upstream and downstream genes you want
# you can changes sys.argv[3] to Dsim, Dsec, Dyak, Dere, Dana, etc. ; only limited by what is in the ortholog file


### STARTING MATERIALS:
#parsed down gff3 file from script 1, called " ", found " "

#need sys to read argv (ie. to read in modified gff3 file from command line)
#

### ENDING MATERIALS:

# no ending materials as of now, but anticipating a quality score thing, and a "work from here" file

### SYS ARGV:
# [1] input file, modified gff3
# [2] gene orthologs file from flybase
# [3] shortened name of nonmel fly (the way that flybase calls is, Dsec = Drosophila sechellia


import sys

mod_gff3 = sys.argv[1]
flyname = sys.argv[3]

#creating object "out" that I will write, 'w' too
ncRNA_flank = open ('2_out_%s_flanking_genes_ncRNA.txt' %flyname, 'w')


# initiating dictionaries:

# dmel ncRNA from gff ('ncRNA_gff_dict')

# front back gene ids ('fbgn_id_dict') 

#
ncRNA_gff_dict = {}

fbgn_id_dict = {}

#initiating ncRNA list
ncRNA = []

#this window length (how many upstream and downstream genes) can change
window_length = 4 #normally two, changing it to 3 for now (many genes don't have up and down)
# 4 doesn't work, something going on at line 88

#opening the file sys[1] as f
with open (mod_gff3, 'r') as f:
	#making gff_list of each line, which is also a list, so I end with a list of lists, this is most efficient way that I know how to do this.
	gff_list = [line.strip().split('\t') for line in f]
	f.close()
	
	#gff_list ex/:
	#[['2L', 'FlyBase', 'gene', '7529', '9484', '.', '+', '.', 'ID=FBgn0031208;Name=CG11023;Ontology_term=SO:0000010,SO:0000087,GO:0016929,GO:0016926;Dbxref=FlyBase:FBan0011023,FlyBase_Annotation_IDs:CG11023,GB_protein:ACZ94128,GB_protein:AAO41164,GB:AI944728,GB:AJ564667,GB_protein:CAD92822,GB:BF495604,UniProt/TrEMBL:Q86BM6,INTERPRO:IPR003653,GB_protein:AGB92323,UniProt/TrEMBL:M9PAY1,OrthoDB7_Drosophila:EOG796K1P,OrthoDB7_Diptera:EOG7X1604,EntrezGene:33155,UniProt/TrEMBL:E1JHP8,UniProt/TrEMBL:Q6KEV3,OrthoDB7_Insecta:EOG7Q8QM7,OrthoDB7_Arthropoda:EOG7R5K68,OrthoDB7_Metazoa:EOG7D59MP,InterologFinder:33155,BIOGRID:59420,FlyAtlas:CG11023-RA,GenomeRNAi:33155;gbunit=AE014134;derived_computed_cyto=21A5-21A5'], 
	# ['2L', 'FlyBase', 'gene', '9839', '21376', '.', '-', '.', 'ID=FBgn0002121;Name=l(2)gl;fullname=lethal (2) giant larvae;Alias=Lgl,lgl,lethal giant larvae,lethal giant larve,lethal giant larva,lethal(2)giant larvae,Complementation group 2.1,Lethal Giant Larvae,dlgl,p127l(2)gl,LGL,l(2) giant larva,CG2671,L(2)GL,p127,l(2)giant larvae,D-LGL,l(2),gl,l[[2]]gl,l-gl,lethal-giant-larvae,Lethal giant larvae,Lethal (2) giant larvae,L(2)gl,Lethal (2) giant larva,Lethal-giant-larvae,MENE (2L)-B,lethal(2) giant larvae,p127[l(2)gl],lethal(2)-giant larvae,lethal-2-giant larvae,l(2) giant larvae,lethal- giant-larvae,Lethal(2)giant larvae,Lethal-2-giant larvae;Ontology_term=SO:0000010,SO:0000087,GO:0005578,GO:0005886,GO:0007269,GO:0016082,GO:0008021,GO:0008283,GO:0016334,GO:0016336,GO:0016333,GO:0016335,GO:0016327,GO:0005829,GO:0045175,GO:0016332,GO:0045184,GO:0007399,GO:0005938,GO:0005737,GO:0007179,GO:0045197,GO:0045196,GO:0002009,GO:0005918,GO:0008105,GO:0045167,GO:0008104,GO:0045746,GO:0007423,GO:0008285,GO:0001738,GO:0016323,GO:0007391,GO:0005856,GO:0030154,GO:0042127,GO:0005614,GO:0045159,GO:0035072,GO:0007559,GO:0045200,GO:0008360,GO:0019991,GO:0007406,GO:0051726,GO:0051668,GO:0007314,GO:0016325,GO:0030036,GO:0030863,GO:0035070,GO:0055059,GO:0035212,GO:0035293,GO:0090163,GO:0048730,GO:0000132,GO:0098725,GO:0060429,GO:0007293,GO:0045176,GO:0072697,GO:0000149,SO:0000548,GO:0005920,GO:0017022,GO:0004860,GO:0006469;Dbxref=FlyBase:FBan0002671,FlyBase_Annotation_IDs:CG2671,INTERPRO:IPR015943,GB_protein:AAN10503,GB_protein:AAG22256,GB_protein:AAN10502,GB_protein:AAN10501,GB_protein:AAF51570,GB_protein:AAG22255,INTERPRO:IPR017986,GB:AA246243,GB:AW942062,GB:AY051654,GB_protein:AAK93078,GB:BH809482,GB:CZ471313,GB:CZ482024,GB:CZ484691,GB:M17022,GB_protein:AAA28671,GB_protein:AAA28672,GB:X05426,GB_protein:CAA29007,UniProt/Swiss-Prot:P08111,INTERPRO:IPR000664,INTERPRO:IPR001680,INTERPRO:IPR013577,GB_protein:AGB92324,UniProt/TrEMBL:M9NCX1,UniProt/TrEMBL:M9PBJ2,OrthoDB7_Drosophila:EOG7CW2GT,OrthoDB7_Diptera:EOG7DRVK2,GB_protein:AFH03479,GB_protein:AFH03478,GB_protein:AFH03481,GB_protein:AFH03480,EntrezGene:33156,INTERPRO:IPR013905,BDGP_clone:PC00404,OrthoDB7_Insecta:EOG7SRGKH,OrthoDB7_Arthropoda:EOG7ZDD82,OrthoDB7_Metazoa:EOG79W94C,InterologFinder:33156,FlyAtlas:CG2671-RB,BIOGRID:59421,Fly-FISH:CG2671,GenomeRNAi:33156,INTERACTIVEFLY:/cytoskel/lethl2g1.htm;gbunit=AE014134;derived_computed_cyto=21A5-21A5'], 
	# ['2L', 'FlyBase', 'gene', '21823', '25155', '.', '-', '.', 'ID=FBgn0031209;Name=Ir21a;fullname=Ionotropic receptor 21a;Alias=CT8983,IR21a,CG2657,DmelIR21a,ionotropic receptor 21a,ir21a;Ontology_term=SO:0000010,SO:0000087,GO:0015276,GO:0016021,GO:0050907,GO:0004970,GO:0016020;Dbxref=FlyBase:FBan0002657,FlyBase_Annotation_IDs:CG2657,GB_protein:AAF51569,GB:CZ468165,UniProt/TrEMBL:Q9VPI2,INTERPRO:IPR001320,OrthoDB7_Drosophila:EOG77MN6S,OrthoDB7_Diptera:EOG725R0H,EntrezGene:33157,OrthoDB7_Insecta:EOG7GZ19M,FlyAtlas:Stencil:2L:25151:23928:GENSCAN%3BCG2657-RA,GenomeRNAi:33157;gbunit=AE014134;derived_computed_cyto=21A5-21B1'], ['2L', 'FlyBase', 'gene', '21952', '24237', '.', '+', '.', 'ID=FBgn0263584;Name=CR43609;Ontology_term=SO:0000011,SO:0000087,SO:0000077;Dbxref=FlyBase_Annotation_IDs:CR43609,EntrezGene:12797867,GenomeRNAi:12797867;derived_computed_cyto=21A5-21B1'], ['2L', 'FlyBase', 'ncRNA', '21952', '24237', '.', '+', '.', 'ID=FBtr0309810;Name=CR43609-RA;Parent=FBgn0263584;Dbxref=FlyBase_Annotation_IDs:CR43609-RA,REFSEQ:NR_047865;score_text=Strongly Supported;score=9'], ['2L', 'FlyBase', 'gene', '25402', '65404', '.', '-', '.', 'ID=FBgn0051973;Name=Cda5;fullname=Chitin deacetylase-like 5;Alias=FBgn0031210,FBgn0031211,FBgn0063656,CG2761,CG2776,BcDNA:RH43162,DmCDA5,CG31973;Ontology_term=SO:0000010,SO:0000087,GO:0005975,GO:0006030,GO:0016810,GO:0005576,GO:0008061,SO:0000548;Dbxref=FlyBase:FBan0031973,FlyBase_Annotation_IDs:CG31973,GB_protein:AAF51567,GB_protein:ABI31281,GB_protein:AAF51568,GB:AY129461,GB_protein:AAM76203,GB:BI613902,GB:CZ468962,GB:CZ486696,GB:CZ486697,UniProt/TrEMBL:Q9VPI3,UniProt/TrEMBL:Q9VPI4,INTERPRO:IPR002509,INTERPRO:IPR002557,INTERPRO:IPR011330,UniProt/TrEMBL:M9NEJ4,UniProt/TrEMBL:M9NCL3,UniProt/TrEMBL:M9NE36,UniProt/TrEMBL:M9NCX5,OrthoDB7_Drosophila:EOG72GCS8,OrthoDB7_Diptera:EOG7WHTSK,GB_protein:AFH03484,GB_protein:AFH03485,GB_protein:AFH03483,GB_protein:ABV53594,GB_protein:AFH03482,EntrezGene:33158,UniProt/TrEMBL:Q0E8V4,UniProt/TrEMBL:A8DYS5,BDGP_clone:RH43162,UniProt/TrEMBL:Q8MQI4,OrthoDB7_Insecta:EOG7F877T,OrthoDB7_Arthropoda:EOG713694,OrthoDB7_Metazoa:EOG7JHM4H,InterologFinder:33158,BIOGRID:59423,FlyAtlas:CG31973-RA,GenomeRNAi:33158;gbunit=AE014134;derived_compute
	
	#making an additional list of ncRNAs, only, actually is id ex/ "ID=FBtr0309810"
	for i in gff_list: #this is just finding all the ncRNAS, a list of lists situation
		if i[2] == 'ncRNA':
			# i = ['2L'[0], 'FlyBase', 'ncRNA', '21952'[3], '24237'[4], '.', '+', '.', 'ID=FBtr0309810;Name=CR43609-RA;Parent=FBgn0263584;Dbxref=FlyBase_Annotation_IDs:CR43609-RA,REFSEQ:NR_047865;score_text=Strongly Supported;score=9']
			#print i
			
			preidRNA = i[8].split(';')[0] #this is to pull out relevant info in last column
			ncRNA.append(preidRNA)
	#make an ncRNA_gff_dict so I can get relevant chromosome_arm, start, and stop info
	for r in ncRNA:
		for i in gff_list:
			data = i[8].split(';')[0] #honing in on relevant information
			if data == r:
				post_data = data.split('=')[1] #want just id (ie, FBgn0031208), not id=blah (ie, ID=FBgn0031208)
				ncRNA_gff_dict[post_data]= i[0], i[3], i[4]
				index = gff_list.index(i)# indexing so more efficient to move backward
				upstream = 0
				counter = 0
				downstream = 0
				anticounter = 0
				ticker=0
				#initiating upstream down stream list
				up,down = [],[]				
	###getting window length # upstream and downstream gene stuff
				while upstream < window_length:
					#print index
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
					#print index
					#just adding 1 for each iteration
					anticounter = anticounter +1
#					print index+anticounter
#					print gff_list[index+anticounter][2]
					#the problem is that last set of ncRNA downstream genes only goes to 3
					try:
						test= gff_list[index+anticounter][2]
						#print "this is test", test
						if gff_list[index+anticounter][2] == 'gene':
							pre_info2 = gff_list[index+anticounter][8].split(';')[0]
							info2 =pre_info2.split('=')[1]
					
							downstream = downstream + 1
						else:
							continue
					except IndexError:
						#print "there was an index error"
						downstream = downstream + 1
						continue
						
					down.append(info2,)
				
				idRNA = r.split('=')[1]
				#"front-back-gene-id-dict"
				fbgn_id_dict[idRNA]= (up,down) # adding information to the dictionary
				#del up_down[:] #... then removing it from the list

#print fbgn_id_dict
#exit()


#ncRNA_gff_dict ex/ #mel-ncRNA
#'FBtr0336987': ('2L', '248909', '249812'), 
#'FBtr0309810': ('2L', '21952', '24237'), 
#'FBtr0345733': ('2L', '66318', '66524'), 

#fbgn_id_dict ex/
#'FBtr0336987': (['FBgn0265149', 'FBgn0262252', 'FBgn0031235', 'FBgn0263465'], ['FBgn0022246', 'FBgn0031238', 'FBgn0031239', 'FBgn0265150']), 
#'FBtr0309810': (['FBgn0263584', 'FBgn0031209', 'FBgn0002121', 'FBgn0031208'], ['FBgn0051973', 'FBgn0267987', 'FBgn0266878', 'FBgn0266879']),
			
	#for r in ncRNA:	

# this is to remove ... actually... not sure... 
# I think this is to find the coordinates for mel genes, which are in a set 
#(so we don't find the coordinates twice for the same gene)		
mel_gene_list = set()
for k,v in fbgn_id_dict.iteritems():
	for mg in v[0]:
		mel_gene_list.add(mg)
	for mg in v[1]:
		mel_gene_list.add(mg)
#print mel_gene_list


#putting the mel gene homologs? with their relevant coordinates in the appropriate species
mapping = dict() #mapping is a dict
with open(sys.argv[2], 'r') as orthos:
	for line in orthos:
		if not line.startswith('#') and not line.startswith('\n'):
			data = line.split('\t')
			#print data
			#print "This is data[5]:", data[5]
			#print "This is data[6]:", data[6]
			
			#sys.argb[3] has been set as 'Dsec', 'Dsim', or whatevs
			if sys.argv[3] in data[6]:
				#just connect dmel gene with d(other) coordinates
				#narrowing it down to mel genes that are upstream and downstream of the lncRNA
				if data[0] in mel_gene_list:
					#info2 =pre_info2.split('=')[1]
					coord = data[8].split ("..")
					mapping[data[0]] = data[5], data[7], coord[0], coord[1]
orthos.close()	
#print mapping
# mapping ex/
#FBgn0001142': ('FBgn0171629', 'scaffold_14', '131505', '132704'), 
#'FBgn0002121': ('FBgn0171620', 'scaffold_14', '19278', '25107'), 
#'FBgn0022246': ('FBgn0171602', 'scaffold_14', '248707', '249826'), 
#'FBgn0025683': ('FBgn0171613', 'scaffold_14', '122003', '129535'), 
#'FBgn0025686': ('FBgn0168201', 'scaffold_841', '3036', '6791'), 

rna_ortho = dict()
for k, v in fbgn_id_dict.iteritems():
	before, after = [], []
	for i,gene1 in enumerate(v[0]):
		ortho_gene = mapping.get(gene1, None)
		if ortho_gene is None:
			continue
		if i <= window_length:
			#print "This is ortho_gene:", ortho_gene
			before.append(ortho_gene)	
			
	for i, gene2 in enumerate(v[1]):
		ortho_gene = mapping.get(gene2, None)
		if ortho_gene is None:
			continue
		if i <= window_length:
			after.append(ortho_gene)
		
		# not sure if this two should be hard coded
		#maybe this should be the same as 'window_length' used earlier? How would that work?
		
		
			#print "This is after ortho_gene:", ortho_gene
	rna_ortho[k] = (before, after)
	# values are tuples with before and after lists
#print rna_ortho
#exit()

###THIS DOESNT WORK RIGHT HERE
#bad_keys = []
#for k, v in rna_ortho.iteritems():
#	#print "This is v:", v
#	print "This is v[0]:", v[0]
#	print "This is v[1]:", v[1]
	#print "\n"
#	for x in v[0]:
#		#print x[1]
#		for y in v[1]: 
#			#print y[1]
#			if x[1] != y[1]:
#				print " %s not %s" %(x[1],y[1])
#				bad_keys.append(k)

# if they are on different scaffolds
#this is important information for quality control	
#del_keys = set(bad_keys)	
#print "del_keys:", del_keys

### stuck here!	
#print rna_ortho

#for i in del_keys:
#	del rna_ortho[i]
	
#rna_ortho ex/
#FBtr0344032 ([('FBgn0171613', 'scaffold_14', '122003..129535')], [('FBgn0171614', 'scaffold_14', '115786..118485'), ('FBgn0171629', 'scaffold_14', '131505..132704')])
#FBtr0309810 ([('FBgn0171619', 'scaffold_14', '27012..28335')], [('FBgn0171620', 'scaffold_14', '19278..25107'), ('FBgn0171618', 'scaffold_14', '31600..58198')])
#FBtr0345733 ([], [('FBgn0171624', 'scaffold_14', '67886..70036'), ('FBgn0171625', 'scaffold_14', '71736..75558')])
#FBtr0345732 ([], [('FBgn0171618', 'scaffold_14', '31600..58198'), ('FBgn0171624', 'scaffold_14', '67886..70036')])
#FBtr0344053 ([('FBgn0171624', 'scaffold_14', '67886..70036')], [('FBgn0171625', 'scaffold_14', '71736..75558'), ('FBgn0171617', 'scaffold_14', '82723..86571')])

#next up, compare Dmel ncRNA with non-mel --- maybe a set?

###haven't used below, maybe i need to. Do I need a set?
#diff = set(ncRNA_gff_dict.keys())- set(rna_ortho.keys())
#print diff


nc_coord_dict = dict()
for key in rna_ortho.keys():
	#print key
	if key in ncRNA_gff_dict:
		#print "*****" * 8
		#print key, ncRNA_gff_dict[key], rna_ortho[key] #prints	dmel and d other coordinates
		#print key, rna_ortho[key]
		# This is rna_ortho[key]: ([], [('FBgn0171639', 'scaffold_14', '231564', '231935'), ('FBgn0171602', 'scaffold_14', '248707', '249826'), ('FBgn0171601', 'scaffold_14', '250192', '251903'), ('FBgn0171599', 'scaffold_14', '261281', '264250')])
		up, down = [],[]
		#print "####" * 10
		#print "? upstream for key:", key
		#print "This is rna_ortho[%s]" %(key), rna_ortho[key]
		#This is rna_ortho[FBtr0336988] ([('FBgn0171599', 'scaffold_14', '261281', '264250')], [('FBgn0171601', 'scaffold_14', '250192', '251903'), ('FBgn0171640', 'scaffold_14', '266985', '269107'), ('FBgn0171598', 'scaffold_14', '269747', '271438'), ('FBgn0171641', 'scaffold_14', '271807', '272472')])
		#print "This is rna_ortho[key][0]", rna_ortho[key][0]
		#print "This is rna_ortho[key][1]", rna_ortho[key][1]
		#print type(rna_ortho[key][0])
		#this doesn't work right for one upstream situations
		#something i need to address
		#print "This is rna_ortho[key][0][1]", rna_ortho[key][0][1]
		if len(rna_ortho[key][0]) != 0:
			up.append(rna_ortho[key][0][0][1])
			up.append(rna_ortho[key][0][0][3])
			#keep this for later
			#try:
				#print "Trying"
				#print rna_ortho[key][0][1]
				#print rna_ortho[key][0][1][1]
				#print 'This is [0][1][1]', rna_ortho[key][0][1][1]
				#up.append(rna_ortho[key][0][1][3])
				#print "This is up", up
				
				
			#except IndexError:
				#continue
				#print "no more upstreams"
				
		
		else:
			print "skipped it all, no upstream!", rna_ortho[key][0]
		
		
		#print "%%%%" * 10
		#print "? downstream for key:", key
		#print "This is rna_ortho[key][1]", rna_ortho[key][1]
		
		if len(rna_ortho[key][1]) != 0:
			
			#print rna_ortho[key][1][0][1]
			#print rna_ortho[key][1][0]
			down.append(rna_ortho[key][1][0][1])
			down.append(rna_ortho[key][1][0][2])
			#print "This is the most updated down", down
			
			#keep this for later
			#try:
				#print "trying the next one..."
				#print "This is rna_ortho[%s][1][1]:" %(key), rna_ortho[key][1][1]
			#except IndexError:
				#continue
				#print "no 2nd downstream!"
		
		else:
			print "skipped it all, no downstream!", rna_ortho[key][1]
		
		nc_coord_dict[key] = (up, down)
		
#print nc_coord_dict
#'FBtr0336987': (['scaffold_14', '231935'], ['scaffold_14', '248707']),
#'FBtr0309810': (['scaffold_14', '28335'], ['scaffold_14', '31600']), 
#'FBtr0345733': (['scaffold_14', '58198'], ['scaffold_14', '67886']), 
#'FBtr0345732': (['scaffold_14', '58198'], ['scaffold_14', '67886']), 
#'FBtr0344053': (['scaffold_14', '70036'], ['scaffold_14', '71736']), 
#'FBtr0344052': (['scaffold_14', '70036'], ['scaffold_14', '71736']), 
#'FBtr0347585': (['scaffold_14', '58198'], ['scaffold_14', '67886']), 
#'FBtr0347595': (['scaffold_14', '284282'], ['scaffold_841', '3036']),
#'FBtr0347594': (['scaffold_14', '284282'], ['scaffold_841', '3036']),
#'FBtr0336988': (['scaffold_14', '264250'], ['scaffold_14', '266985']), 
#'FBtr0336984': (['scaffold_14', '240248'], ['scaffold_14', '231564']), 
#'FBtr0336985': (['scaffold_14', '240248'], ['scaffold_14', '231564']), 
#'FBtr0336986': (['scaffold_14', '240248'], ['scaffold_14', '231564']), 
#'FBtr0344032': (['scaffold_14', '129535'], ['scaffold_14', '131505']), 
#'FBtr0336836': (['scaffold_14', '132704'], ['scaffold_14', '138261']), 
#'FBtr0336837': (['scaffold_14', '132704'], ['scaffold_14', '138261']

for k, v in nc_coord_dict.iteritems():
	print "This is k:", k
	print "This is v:", v
	#print "This is v[0]:", v[0]
	#print "This is v[1]:", v[1]
	if len(v[0]) !=0 and len(v[1]) != 0:
		if v[0][0] == v[1][0]:
			#print "%s = %s" %(v[0][0], v[1][0])
			ncRNA_flank.write("%s\t%s\t%s\t%s\n" %(k, v[0][0], v[0][1], v[1][1]))
			#ncRNA_gene_out.write('%s' %line)
			
		else:
			print "%s != %s" %(v[0][0], v[1][0])
			#this might be important for quality control
			#think about how to address inversions?
		
	#print "\n"
#	for x in v[0]:
#		#print x[1]
#		for y in v[1]: 
#			#print y[1]
#			if x[1] != y[1]:
#				print " %s not %s" %(x[1],y[1])
#				bad_keys.append(k)





#downstream from here I want 1) to do something to measure how good a hit is (something like that)
#2) to move forward to the next step.  Looks like 
#<dmel lncRNA><non-mel chr, start, stop>
#<FBtr0344032><scaffold_14, 129535, 115786 (inversion?),
ncRNA_flank.close()
exit()
#then get dmel ncRNA coordinates? Do I need this
# What do I want this output to look like?


