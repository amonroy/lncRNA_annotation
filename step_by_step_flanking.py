## this is me going through and looking closely at 2_flanking_genes_ncRNA_2_v1.py
## right now I am trying to reduced excess computing by just running through all the fly types at once.
## need to make the fly input interactive
## when do I make the script stop?
## should I only let it input those 12?

# to run this code in Desktop/Projects/lncRNA directory on my laptop
#python lncRNA_annotation_doc/step_by_step_flanking.py out/1_prac_ncRNA_and_genes.txt flybase/gene_orthologs_fb_2016_03.tsv Dsec

#also, probably a good idea to turn list of tuples into a list of list...

# the parsed down gff only has genes and ncRNA.


#python lncRNA_annotation_doc/step_by_step_flanking.py out/1_2_prac_ncRNA_and_genes.txt flybase/gene_orthologs_fb_2016_03.tsv Dsec out/dmel-on-dsec-2.blstn
#              sys[0]                                             sys[1]                      sys[2]                        sys[3]    sys[4]


import sys
import datetime
import itertools
import os
today = datetime.date.today()
window_length = 4 #assigns how many upstream genes and downstream genes are desired

#fly_types_list = ['Dsim', 'Dsec', 'Dyak', 'Dere', 'Dana', 'Dpse', 'Dper', 'Dwil', 'Dmoj', 'Dvir', 'Dgri']
#what fly?
# ok or sorry that does not exist
#try again
# how to get it to quit?

#I think I should change this to be interactive
fly = sys.argv[3]
#this works
def gff_list():
	"""This function takes the modified gff3 file and creates a list"""
	mod_gff3 = sys.argv[1]
	with open(mod_gff3, 'r') as f:
		gff = [line.strip().split('\t') for line in f]
		f.close()
	return gff
	  #gff_list ex/:
	#['2L', 'FlyBase', 'gene', '7529', '9484', '.', '+', '.', 'ID=FBgn0031208;Name=CG11023;Ontology_term=SO:0000010,SO:0000087,GO:0016929,GO:0016926;Dbxref=FlyBase:FBan0011023,FlyBase_Annotation_IDs:CG11023,GB_protein:ACZ94128,GB_protein:AAO41164,GB:AI944728,GB:AJ564667,GB_protein:CAD92822,GB:BF495604,UniProt/TrEMBL:Q86BM6,INTERPRO:IPR003653,GB_protein:AGB92323,UniProt/TrEMBL:M9PAY1,OrthoDB7_Drosophila:EOG796K1P,OrthoDB7_Diptera:EOG7X1604,EntrezGene:33155,UniProt/TrEMBL:E1JHP8,UniProt/TrEMBL:Q6KEV3,OrthoDB7_Insecta:EOG7Q8QM7,OrthoDB7_Arthropoda:EOG7R5K68,OrthoDB7_Metazoa:EOG7D59MP,InterologFinder:33155,BIOGRID:59420,FlyAtlas:CG11023-RA,GenomeRNAi:33155;gbunit=AE014134;derived_computed_cyto=21A5-21A5'], ['2L', 'FlyBase', 'gene', '9839', '21376', '.', '-', '.', 'ID=FBgn0002121;Name=l(2)gl;fullname=lethal (2) giant larvae;Alias=Lgl,lgl,lethal giant larvae,lethal giant larve,lethal giant larva,lethal(2)giant larvae,Complementation group 2.1,Lethal Giant Larvae,dlgl,p127l(2)gl,LGL,l(2) giant larva,CG2671,L(2)GL,p127,l(2)giant larvae,D-LGL,l(2),gl,l[[2]]gl,l-gl,lethal-giant-larvae,Lethal giant larvae,Lethal (2) giant larvae,L(2)gl,Lethal (2) giant larva,Lethal-giant-larvae,MENE (2L)-B,lethal(2) giant larvae,p127[l(2)gl],lethal(2)-giant larvae,lethal-2-giant larvae,l(2) giant larvae,lethal- giant-larvae,Lethal(2)giant larvae,Lethal-2-giant larvae;Ontology_term=SO:0000010,SO:0000087,GO:0005578,GO:0005886,GO:0007269,GO:0016082,GO:0008021,GO:0008283,GO:0016334,GO:0016336,GO:0016333,GO:0016335,GO:0016327,GO:0005829,GO:0045175,GO:0016332,GO:0045184,GO:0007399,GO:0005938,GO:0005737,GO:0007179,GO:0045197,GO:0045196,GO:0002009,GO:0005918,GO:0008105,GO:0045167,GO:0008104,GO:0045746,GO:0007423,GO:0008285,GO:0001738,GO:0016323,GO:0007391,GO:0005856,GO:0030154,GO:0042127,GO:0005614,GO:0045159,GO:0035072,GO:0007559,GO:0045200,GO:0008360,GO:0019991,GO:0007406,GO:0051726,GO:0051668,GO:0007314,GO:0016325,GO:0030036,GO:0030863,GO:0035070,GO:0055059,GO:0035212,GO:0035293,GO:0090163,GO:0048730,GO:0000132,GO:0098725,GO:0060429,GO:0007293,GO:0045176,GO:0072697,GO:0000149,SO:0000548,GO:0005920,GO:0017022,GO:0004860,GO:0006469;Dbxref=FlyBase:FBan0002671,FlyBase_Annotation_IDs:CG2671,INTERPRO:IPR015943,GB_protein:AAN10503,GB_protein:AAG22256,GB_protein:AAN10502,GB_protein:AAN10501,GB_protein:AAF51570,GB_protein:AAG22255,INTERPRO:IPR017986,GB:AA246243,GB:AW942062,GB:AY051654,GB_protein:AAK93078,GB:BH809482,GB:CZ471313,GB:CZ482024,GB:CZ484691,GB:M17022,GB_protein:AAA28671,GB_protein:AAA28672,GB:X05426,GB_protein:CAA29007,UniProt/Swiss-Prot:P08111,INTERPRO:IPR000664,INTERPRO:IPR001680,INTERPRO:IPR013577,GB_protein:AGB92324,UniProt/TrEMBL:M9NCX1,UniProt/TrEMBL:M9PBJ2,OrthoDB7_Drosophila:EOG7CW2GT,OrthoDB7_Diptera:EOG7DRVK2,GB_protein:AFH03479,GB_protein:AFH03478,GB_protein:AFH03481,GB_protein:AFH03480,EntrezGene:33156,INTERPRO:IPR013905,BDGP_clone:PC00404,OrthoDB7_Insecta:EOG7SRGKH,OrthoDB7_Arthropoda:EOG7ZDD82,OrthoDB7_Metazoa:EOG79W94C,InterologFinder:33156,FlyAtlas:CG2671-RB,BIOGRID:59421,Fly-FISH:CG2671,GenomeRNAi:33156,INTERACTIVEFLY:/cytoskel/lethl2g1.htm;gbunit=AE014134;derived_computed_cyto=21A5-21A5'],
	# ['2L', 'FlyBase', 'ncRNA', '286383', '288292', '.', '+', '.', 'ID=FBtr0347595;Name=CR46263-RA;Parent=FBgn0267996;Dbxref=FlyBase_Annotation_IDs:CR46263-RA;score_text=Weakly Supported;score=0'], ['2L', 'FlyBase', 'gene', '287252', '289144', '.', '-', '.', 'ID=FBgn0025686;Name=Amnionless;fullname=Amnionless ortholog;Alias=FBgn0031246,CG11592,CK02467,BEST:CK02467,dAMN,Amnionless;Ontology_term=SO:0000010,SO:0000087,GO:0046331,GO:0097206,GO:0016021,GO:0097017;Dbxref=FlyBase:FBan0011592,FlyBase_Annotation_IDs:CG11592,GB_protein:AAF51514,GB:AA141784,GB:CZ468687,UniProt/TrEMBL:Q9VPN2,GB_protein:AGB92350,OrthoDB7_Drosophila:EOG7CGKJK,EntrezGene:33199,BDGP_clone:IP03221,OrthoDB7_Diptera:EOG774804,INTERPRO:IPR026112,OrthoDB7_Insecta:EOG7G266G,OrthoDB7_Arthropoda:EOG7P65FW,OrthoDB7_Metazoa:EOG7ZGX2W,InterologFinder:33199,FlyAtlas:CG11592-RA,GenomeRNAi:33199;gbunit=AE014134;derived_computed_cyto=21B7-21B7'], ['2L', 'FlyBase', 'gene', '292419', '293222', '.', '+', '.', 'ID=FBgn0031247;Name=CG11562;Alias=FBgn0063011,BcDNA:RE44650;Ontology_term=SO:0000010,SO:0000087,GO:0005739,GO:0003674,GO:0008150;Dbxref=FlyBase:FBan0011562,FlyBase_Annotation_IDs:CG11562,GB_protein:AAF51513,GB:AI520524,GB:AI945841,GB:AY119645,GB_protein:AAM50299,GB:BE662187,GB:BI358003,UniProt/TrEMBL:Q9VPN3,OrthoDB7_Drosophila:EOG7HTW3H,OrthoDB7_Diptera:EOG7200K9,EntrezGene:33200,BDGP_clone:RE44650,OrthoDB7_Insecta:EOG7B9454,OrthoDB7_Arthropoda:EOG7RK278,OrthoDB7_Metazoa:EOG78H3X3,FlyAtlas:CG11562-RA,INTERPRO:IPR031568,Fly-FISH:CG11562,GenomeRNAi:33200;gbunit=AE014134;derived_computed_cyto=21B7-21B7'], ['2L', 'FlyBase', 'gene', '292959', '294681', '.', '-', '.', 'ID=FBgn0017457;Name=U2af38;fullname=U2 small nuclear riboprotein auxiliary factor 38;Alias=FBgn0010626,U2AF38,U2AF,dU2AF38,DU2AF38,CG3582,dU2AF[38],l(2)06751,u2af38,U2AF 38;Ontology_term=GO:0089701,SO:0000010,SO:0000087,GO:0000398,GO:0008187,GO:0005681,GO:0005686,GO:0000381,GO:0005634,GO:0003729,GO:0007052,GO:0071011,GO:0008380,GO:0000166,GO:0046872;Dbxref=FlyBase:FBan0003582,FlyBase_Annotation_IDs:CG3582,GB_protein:AAF51512,GB:AA264081,GB:AA820431,GB:AC004115,GB:AC008371,GB:AI061776,GB:AI455418,GB:AI944553,GB:AQ026079,GB:AY058537,GB_protein:AAL13766,GB:U67066,GB_protein:AAB17271,UniProt/Swiss-Prot:Q94535,INTERPRO:IPR000504,INTERPRO:IPR000571,INTERPRO:IPR009145,INTERPRO:IPR012677,GB_protein:AGB92351,UniProt/TrEMBL:M9PBM1,OrthoDB7_Drosophila:EOG7FRM2M,OrthoDB7_Diptera:EOG700KS6,EntrezGene:33201,BDGP_clone:LD24048,OrthoDB7_Insecta:EOG76QSHP,OrthoDB7_Arthropoda:EOG7KMJ7T,OrthoDB7_Metazoa:EOG70089G,apodroso:10448-U2af38[k14504],InterologFinder:33201,FlyAtlas:CG3582-RA,BIOGRID:59457,Fly-FISH:CG3582,GenomeRNAi:33201;gbunit=AE014134;derived_computed_cyto=21B7-21B8']]

#this works  
def ncRNA_list(list):
	"""This function takes the gff_list and makes an ncRNA_list"""
	ncRNA = [] #initiates list
	for i in list:
		if i[2] == 'ncRNA':
			preidRNA = i[8].split(';')[0]
			#[ID=FBgn0031208];Name=CG11023;Ontology_term=SO:0000010,SO:0000087,GO:0016929,GO:0016926;Dbxref=FlyBase:FBan0011023,FlyBase_Annotation_IDs:CG11023,GB_protein:ACZ94128,GB_protein:AAO41164,GB:AI944728,GB:AJ564667,GB_protein:CAD92822,GB:BF495604,UniProt/TrEMBL:Q86BM6,INTERPRO:IPR003653,GB_protein:AGB92323,UniProt/TrEMBL:M9PAY1,OrthoDB7_Drosophila:EOG796K1P,OrthoDB7_Diptera:EOG7X1604,EntrezGene:33155,UniProt/TrEMBL:E1JHP8,UniProt/TrEMBL:Q6KEV3,OrthoDB7_Insecta:EOG7Q8QM7,OrthoDB7_Arthropoda:EOG7R5K68,OrthoDB7_Metazoa:EOG7D59MP,InterologFinder:33155,BIOGRID:59420,FlyAtlas:CG11023-RA,GenomeRNAi:33155;gbunit=AE014134;derived_computed_cyto=21A5-21A5'
			ncRNA.append(preidRNA)
	return ncRNA
	#['ID=FBtr0309810', 'ID=FBtr0347585', 'ID=FBtr0345732', 'ID=FBtr0345733', 'ID=FBtr0344052', 'ID=FBtr0344053', 'ID=FBtr0344032', 'ID=FBtr0336836', 'ID=FBtr0336837', 'ID=FBtr0336984', 'ID=FBtr0336985', 'ID=FBtr0336986', 'ID=FBtr0336987', 'ID=FBtr0336988', 'ID=FBtr0347594', 'ID=FBtr0347595']
			
####This part is by far the most confusing! Go back through and comment it up			
#this works too! I am great at everything
def ncRNA_gff_dict(rna_list, gff_list, window_length):
	"""This function takes our two lists, ncRNA and gff_list and makes two great dictionaries. fbgn_id_dict where the key is ncRNA and the values are upstream genes and downstream genes. ncRNA_gff_dict where the key is ncRNA gene and the values are chromosome, start, stop, for dmel"""
	#ncRNA_gff_dict = {}
	fbgn_id_dict = {}
	for r in rna_list:
		#print "This is r in rna_list", r
		for i in gff_list:
			data = i[8].split(';')[0]
			if data == r:
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
				#"front-back-gene-id-dict"
				fbgn_id_dict[idRNA] = (up, down)
	return fbgn_id_dict
	#fbgn_id_dict
	#{'FBtr0336987': (['FBgn0265149', 'FBgn0262252', 'FBgn0031235', 'FBgn0263465'], ['FBgn0022246', 'FBgn0031238', 'FBgn0031239', 'FBgn0265150']), 'FBtr0309810': (['FBgn0263584', 'FBgn0031209', 'FBgn0002121', 'FBgn0031208'], ['FBgn0051973', 'FBgn0267987', 'FBgn0266878', 'FBgn0266879']), 'FBtr0345733': (['FBgn0266879', 'FBgn0266878', 'FBgn0267987', 'FBgn0051973'], ['FBgn0067779', 'FBgn0266322', 'FBgn0031213', 'FBgn0031214']), 'FBtr0345732': (['FBgn0266878', 'FBgn0267987', 'FBgn0051973', 'FBgn0263584'], ['FBgn0266879', 'FBgn0067779', 'FBgn0266322', 'FBgn0031213']), 'FBtr0344053': (['FBgn0266322', 'FBgn0067779', 'FBgn0266879', 'FBgn0266878'], ['FBgn0031213', 'FBgn0031214', 'FBgn0002931', 'FBgn0031216']), 'FBtr0344052': (['FBgn0266322', 'FBgn0067779', 'FBgn0266879', 'FBgn0266878'], ['FBgn0031213', 'FBgn0031214', 'FBgn0002931', 'FBgn0031216']), 'FBtr0347585': (['FBgn0267987', 'FBgn0051973', 'FBgn0263584', 'FBgn0031209'], ['FBgn0266878', 'FBgn0266879', 'FBgn0067779', 'FBgn0266322']), 'FBtr0347595': (['FBgn0267996', 'FBgn0267995', 'FBgn0031245', 'FBgn0031244'], ['FBgn0025686', 'FBgn0031247', 'FBgn0017457']), 'FBtr0347594': (['FBgn0267995', 'FBgn0031245', 'FBgn0031244', 'FBgn0003444'], ['FBgn0267996', 'FBgn0025686', 'FBgn0031247', 'FBgn0017457']), 'FBtr0336988': (['FBgn0265150', 'FBgn0031239', 'FBgn0031238', 'FBgn0022246'], ['FBgn0031240', 'FBgn0086912', 'FBgn0086856', 'FBgn0086855']), 'FBtr0336984': (['FBgn0265151', 'FBgn0266557', 'FBgn0031233', 'FBgn0031232'], ['FBgn0265153', 'FBgn0265152', 'FBgn0263465', 'FBgn0031235']), 'FBtr0336985': (['FBgn0265153', 'FBgn0265151', 'FBgn0266557', 'FBgn0031233'], ['FBgn0265152', 'FBgn0263465', 'FBgn0031235', 'FBgn0262252']), 'FBtr0336986': (['FBgn0265152', 'FBgn0265153', 'FBgn0265151', 'FBgn0266557'], ['FBgn0263465', 'FBgn0031235', 'FBgn0262252', 'FBgn0265149']), 'FBtr0344032': (['FBgn0266304', 'FBgn0025683', 'FBgn0031220', 'FBgn0031219'], ['FBgn0001142', 'FBgn0265074', 'FBgn0265075', 'FBgn0051975']), 'FBtr0336836': (['FBgn0265074', 'FBgn0001142', 'FBgn0266304', 'FBgn0025683'], ['FBgn0265075', 'FBgn0051975', 'FBgn0051976', 'FBgn0051974']), 'FBtr0336837': (['FBgn0265075', 'FBgn0265074', 'FBgn0001142', 'FBgn0266304'], ['FBgn0051975', 'FBgn0051976', 'FBgn0051974', 'FBgn0031224'])}

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
#set(['FBgn0001142', 'FBgn0266322', 'FBgn0265075', 'FBgn0002121', 'FBgn0265153', 'FBgn0265152', 'FBgn0022246', 'FBgn0265150', 'FBgn0025683', 'FBgn0025686', 'FBgn0031233', 'FBgn0031232', 'FBgn0031235', 'FBgn0031219', 'FBgn0031239', 'FBgn0031238', 'FBgn0031214', 'FBgn0031213', 'FBgn0031224', 'FBgn0267995', 'FBgn0265074', 'FBgn0267996', 'FBgn0002931', 'FBgn0067779', 'FBgn0003444', 'FBgn0086912', 'FBgn0266557', 'FBgn0266304', 'FBgn0031247', 'FBgn0086856', 'FBgn0086855', 'FBgn0265151', 'FBgn0265149', 'FBgn0263584', 'FBgn0051974', 'FBgn0051975', 'FBgn0051976', 'FBgn0031208', 'FBgn0031209', 'FBgn0031220', 'FBgn0262252', 'FBgn0266879', 'FBgn0266878', 'FBgn0051973', 'FBgn0263465', 'FBgn0267987', 'FBgn0031240', 'FBgn0031244', 'FBgn0031245', 'FBgn0017457', 'FBgn0031216'])														 
def ortho_mapping(set):
	"""This function maps other another species' orthologs to dmel genes"""
	mapping = dict()
	with open(sys.argv[2], 'r') as orthos: #this is finding the  ortho_gene_coords
		for line in orthos:
			if not line.startswith('#') and not line.startswith('\n'):
				data = line.split ('\t')
				###switch this to "fly" sys.argv
				if sys.argv[3] in data[6]:
					if data[0] in mel_genes_obj:
						#print data
						coord = data[8].split("..")
						mapping[data[0]] = data[5], data[7], coord[0], coord[1]
						#exit()
	return mapping
	# mapping ex/
#FBgn0001142': ('FBgn0171629', 'scaffold_14', '131505', '132704'), 
#'FBgn0002121': ('FBgn0171620', 'scaffold_14', '19278', '25107'), 
#'FBgn0022246': ('FBgn0171602', 'scaffold_14', '248707', '249826'), 
#'FBgn0025683': ('FBgn0171613', 'scaffold_14', '122003', '129535'), 
#'FBgn0025686': ('FBgn0168201', 'scaffold_841', '3036', '6791'),			  
				 
def mel_rna_ortho(fbgn_id_dict, map):
	"""This function maps the mel lncRNA to its upstream and downstream orthologs in a different species"""
	rna_ortho = dict()
	for k, v in fbgn_id_dict.iteritems():
		before, after = [], []
		for i, gene1 in enumerate(v[0]):
			ortho_gene = map.get(gene1, None) # get returns a value for a given key
			if ortho_gene is None:
				continue
			if i <= window_length:
				before.append(ortho_gene)
				
		for i, gene2 in enumerate(v[1]):
			ortho_gene = map.get(gene2, None)
			if ortho_gene is None:
				continue
			if i <= window_length:
				after.append(ortho_gene)
		rna_ortho[k] = (before,after)
	return rna_ortho
	#'FBtr0344032': ([('FBgn0171613', 'scaffold_14', '122003', '129535'), ('FBgn0171614', 'scaffold_14', '115786', '118485')], [('FBgn0171629', 'scaffold_14', '131505', '132704'), ('FBgn0171610', 'scaffold_14', '138261', '139727')]), 
	#'FBtr0309810': ([('FBgn0171619', 'scaffold_14', '27012', '28335'), ('FBgn0171620', 'scaffold_14', '19278', '25107'), ('FBgn0177629', 'scaffold_8', '1392069', '1396697')], [('FBgn0171618', 'scaffold_14', '31600', '58198')]), 
	#'FBtr0345733': ([('FBgn0171618', 'scaffold_14', '31600', '58198')], [('FBgn0171624', 'scaffold_14', '67886', '70036'), ('FBgn0171625', 'scaffold_14', '71736', '75558')]),

def final_coord(ortho_dict):#rna_ortho_dict,
	"""This function finds the end of the front gene ortholog and the front of the back gene ortholog, to give a dictionary with a putative start and putative stop"""
	final_coord_dict = dict()
	for k, v in ortho_dict.iteritems():
		upstream = v[0]
		downstream = v[1]
		uscafs = set()
		dscafs = set()
		for gene in upstream:
			uscafs.add(gene[1])
		for gene in downstream:
			dscafs.add(gene[1])
		#print "up", uscafs
		#print "down", dscafs
		common_scaf = uscafs.intersection(dscafs)
		#different_scaf = uscafs.difference(dscafs)
		#print "common", common_scaf
		#if len(different_scaf) >= 1:
		#	print "different", different_scaf, k
		for x in common_scaf:
			#print "This is x", x
			upos = []
			for gene in upstream:
				#print "This is gene", gene
				if gene[1] == x:
					#print "This is gene[1]", gene[1]
					#print "This is common_scaf", common_scaf
					#upos.append(gene[0])
					upos.append(gene[2])
					upos.append(gene[3])
					#print upos
			dpos = []
			for gene in downstream:
				if gene[1] == x:
					#dpos.append(gene[0])
					dpos.append(gene[2])
					dpos.append(gene[3])
		#print "This is upos", upos
		#print "This is upos[1:]", upos[1:]
		#print "This is dpos", dpos
		#ex upos : ['3815439', '3822866', '3808823', '3809996']
		#ex dbos : ['3823313', '3826021', '3826740', '3828621', '3829156', '3829994', '3831313', '3855168']
		#quit()
		#final_coord_dict[k] = [x, max(upos), min(dpos)]
		#print type(upos[1]), its a string, want int!
		upos_num = [int(n) for n in upos]
		#print "This is upos", upos
		#print "this is upos_num", upos_num
		dpos_num = [int(n) for n in dpos]
		merged_pos = upos_num + dpos_num
		print "merged", merged_pos
		print min(merged_pos)
		final_coord_dict[k] = [x, min(merged_pos), max(merged_pos)]
		print k, final_coord_dict[k]
		#avg_upos = sum(upos_num)/len(upos_num)
		#avg_dpos = sum(dpos_num)/len(dpos_num)
		#print k
		#print upos, dpos
		#print max(upos_num), min(dpos_num), min(upos_num), max(upos_num)
		#if max(upos_num) < min(dpos_num):
		#	final_coord_dict[k] = [x, min(upos_num), max(dpos_num)]
		#elif min(upos_num) > max(upos_num):
		#	final_coord_dict[k] = [x, min(dpos_num), max(upos_num)]
		#print k, x, max(upos), min(dpos)
	return final_coord_dict
		#'FBtr0342867': ['scaffold_0', '3442611', '3447776'], 'FBtr0342862': ['scaffold_0', '3442611', '3447776']
	
		#FBtr0344032 Scf_2L 114960 117805
		#FBtr0309810 Scf_2L 4304 17164
		#FBtr0345733 Scf_2L 47089 52876

	#quit()
		
		#v[0], v[1] are lists of tuples
		#FBtr0344032 ([('FBgn0171613', 'scaffold_14', '122003', '129535'), ('FBgn0171614', 'scaffold_14', '115786', '118485')], [('FBgn0171629', 'scaffold_14', '131505', '132704'), ('FBgn0171610', 'scaffold_14', '138261', '139727')])

		#quit()
###PSEUDO CODE:::


### this is potentially not necessary
def get_all_u_scaffs(coord_dict):
	all_unique_scaffs= set()
	for k,v in coord_dict.iteritems():
		#print "This is v, and v[0]", v, v[0]
		all_unique_scaffs.add(v[0])
	return all_unique_scaffs

def karl_read_scaffolds(filepath):
	"""Need to check this works"""
	user_input = raw_input("Enter the path of the scaffold file: ")
	assert os.path.exists(user_input), "I did not find the file at, "+str(user_input)
	print ("Hooray we found your file!")
	#print "I am in karl_read"
	with open(user_input, 'r+') as species:
		key = ''
		sequence = []
		fly = dict()
		for line in species:
			if line.startswith('>'):
				if key and sequence:
					fly[key] = ''.join(sequence)
				key = line.split(' ')[0][1:]
				sequence = []
			else:
				sequence.append(line.rstrip())
		fly[key] = ''.join(sequence)
	return fly
	

def output_sequence(scaf_dict, coord_dict, mel_seqs):
	"""This print the sequence of putative lncRNA in non-mel species, using start and stop and scaffold"""
	for k, v in coord_dict.iteritems():
		#print "This is k", k
		#print "This is v", v
		#print "This is type(v)", type(v)
		#print "This is v[0]", v[0]
		if v[1] < v[2]:
			start = int(v[1])+1
			end = int(v[2])+1
		if v[2] < v[1]:
			start = int(v[2])+1
			end = int(v[1])+1
		print '>'+k+';'+'dmel'+'\n',mel_seqs[k] 
		print '>'+k+';'+fly+'\n',scaf_dict[v[0]][start:end]
		
def ncrna_seq_dict():
	"""This reads in flybase file and makes a dictionary of ncRNA id (key) and sequence (value)"""
	with open('flybase/dmel-all-ncRNA-r6.11.fasta', 'r+') as file:
		key = ''
		sequence = []
		ncRNA = dict()
		for line in file:
			if line.startswith('>'):
				if key and sequence:
					ncRNA[key] = ''.join(sequence)
				key = line.split(' ')[0][1:]
				sequence = []
			else:
				sequence.append(line.rstrip())
		ncRNA[key] = ''.join(sequence)
	return ncRNA
	#'FBtr0345520': 'TCAGTTCAATTTAACTACTCATCGGCAACGGTTTGAATTGAAAATGGCATCTTCGATGGACTTACTCTGGCTTTTTGGCATTTTGGTTATCAGCAGCGTGCTACTGGAAACTTCGACTCTGAACAGCGAAGTTGTGGAGTCCGTAACCGTTGCTTCGAAGTCAAGCAACTTCAGTGAAAATGTGAAAATAAAACTGGAAAAAACCTGGAATCATAATGGAACTATGCAGGCGCAACCGTACAATATTCTGGAGACCAACTGTCCTCATGGCTATGTCCTGGCAAATAAACATTGCCACAAGCGTGCCTAAACGTATTGAAATAAATGTAAATGTTACGTAAAGATAATTTCTGTTATAGCACATTTTGAATTGAAGATCATATGTTCGCATATAATTACACAATTATAAAATATATATGATTTAAAAGCATCTCTGTAAGACAACTTTAAATAGCTCTTGGCATTCTCCTAACTTGGCACTCAATAAAAATTGCACCTTATGGCATTG'

def blast_coord(coord_dict):
	"""This function will find all of the hits of dmel ncRNA against non-mel genome and return the best scaffold and the start and stop"""
	file = sys.argv[4]
	blast_coord_dict = dict()
	final_blast_dict = dict()
	with open(file, 'r') as f:
		for line in f:
			data= line.strip().split('\t')
			#print "This is data:", data
	   # key = data[0]
		#try: 
			if data[0] in blast_coord_dict:
				blast_coord_dict[data[0]].append(data[1:])
			else:
				blast_coord_dict[data[0]]= [data[1:]]
		#except:
			#blast_coord_dict[data[0]]\t
			#data_list.append(data[1:]) 
			#blast_coord_dict[data[0]]= data_list
		#print "This is blast_coord_dict", blast_coord_dict()
	
	
	for k,v in blast_coord_dict.iteritems():
		if k in coord_dict:
		#print "This is k", k
		#print "this is v", v
			sstart = []
			send = []
			for list in v:
				uscafs =set()
				uscafs.add(list[0])
				sstart.append(int(list[7]))
				send.append(int(list[8]))
			#here, since we are finding it compared to the lncRNA, we want more buffer, so the lowest and highest between up and down
			final_blast_dict[k] = uscafs, min(sstart), max(send)
			#want small, large
			#print type(list[7])
			#print type(list[8])
			#if max(list[7]) > min(list[8]):
			#	print '*'+k, max(list[8]), min(list[7])
			#if max(list[7]) < min(list[8]):
			#	print k, max(list[7]), min(list[8])

	return blast_coord_dict, final_blast_dict

#print sys.argv[4]
gff_obj =gff_list()
#print gff_obj
#quit()
ncRNA_obj = ncRNA_list(gff_obj)
#print ncRNA_obj
#quit()
#print ncRNA_obj
fbgn_dict_obj, ncRNA_gff_dict_obj =ncRNA_gff_dict(ncRNA_obj, gff_obj, window_length)
#print fbgn_dict_obj
#print ncRNA_gff_dict_obj
mel_genes_obj = mel_gene_set(fbgn_dict_obj)
#print mel_genes_obj
ortho_map = ortho_mapping(mel_genes_obj)
#print ortho_map
rna_ortho_dict = mel_rna_ortho(fbgn_dict_obj, ortho_map)
#print rna_ortho_dict
#quit()
#python lncRNA_annotation_doc/step_by_step_flanking.py out/1_prac_ncRNA_and_genes.txt 
final_coord_obj = final_coord(rna_ortho_dict)
#print final_coord_obj
rna_seqs_obj = ncrna_seq_dict()
#print rna_seqs_obj
#quit()

#scaffs_set = get_all_u_scaffs(final_coord_obj)
#print "This is scaffs_set", scaffs_set
#input = read_in_file_from_user()
#fly_dict = karl_read_scaffolds(input)
#scaff_dict = karl_read_scaffolds(input)
#for k,v in scaff_dict.iteritems():
#	print k, v
#bls_dict, final_bls_dict =blast_coord(final_coord_obj)

#for k,v in final_coord_obj.iteritems():
#	try:
#		print k, v, final_bls_dict[k]
#	except KeyError:
#		print "%s was causing an error" %(k)


#flybase/dsec-all-chromosome-r1.3.fasta

#This is just to quickly compare the output from the lncRNA keys created from the ncRNA file, versus the ncRNA found in the GFF3 file
#They are all the same length == 2914
#print len(ncRNA_gff_dict_obj)
#print len(rna_seqs_obj)
#print len(ncRNA_obj)


#print len(set(scaff_dict.values()))
##14730
#output_sequence(scaff_dict, final_coord_obj,rna_seqs_obj)
#cd 
#what all-chromosome file for flies (under fasta sequence)
#
#print scaffs_set

##DGRP? initial pass to compare against mel.
	## how well do you do?
	## might not be any variation, but intersting if it is
	## each individual raleigh line is its own species
	
##classify by structure
#PCA multidimensional scaling? what can the computer use to handle it?
#naively look for structure against all 10 species


#scaffold_1 type=golden_path_region; loc=scaffold_1:1..14215200; ID=scaffold_1; dbxref=GB:CH480816,REFSEQ:NW_001999690; MD5=390802af22260f1bc1b969b928cafa29; length=14215200; release=r1.3; species=Dsec;
#>scaffold_0 type=golden_path_region; loc=scaffold_0:1..21120651; ID=scaffold_0; dbxref=GB:CH480815,REFSEQ:NW_001999689; MD5=6d2fa4d8ab670b07e9100cce3fdb1c22; length=21120651; release=r1.3; species=Dsec;