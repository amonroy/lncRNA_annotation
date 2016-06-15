## this is me going through and looking closely at 2_flanking_genes_ncRNA_2_v1.py

import sys
window_length = 4 #assigns how many upstream genes and downstream genes are desired

#this works
def gff_list():
    """This function takes the modified gff3 file and creates a list"""
    mod_gff3 = sys.argv[1]
    with open(mod_gff3, 'r') as f:
        gff = [line.strip().split('\t') for line in f]
        f.close()
    return gff
        
        
        #gff_list ex/:
	    #[['ID=FBtr0309810', 'ID=FBtr0347585', 'ID=FBtr0345732', 'ID=FBtr0345733', 'ID=FBtr0344052', 'ID=FBtr0344053', 'ID=FBtr0344032', 'ID=FBtr0336836', 'ID=FBtr0336837', 'ID=FBtr0336984', 'ID=FBtr0336985', 'ID=FBtr0336986', 'ID=FBtr0336987', 'ID=FBtr0336988', 'ID=FBtr0347594', 'ID=FBtr0347595']
	  
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
            
####This part is by far the most confusing! Go back through and comment it up            
#this works too! I am great at everything
def ncRNA_gff_dict(rna_list, gff_list, window_length):
    """This function takes our two lists, ncRNA and gff_list and makes a great dictionary"""
    ncRNA_gff_dict = {}
    fbgn_id_dict = {}
    for r in rna_list:
        #print "This is r in rna_list", r
        for i in gff_list:
            data = i[8].split(';')[0]
            if data == r:
                #print r, i #
                post_data = data.split('=')[1] # ex/ FBgn0031208
                ncRNA_gff_dict[post_data] = i[0], i[3], i[4]
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
    return fbgn_id_dict, ncRNA_gff_dict
    #fbgn_id_dict
    #{'FBtr0336987': (['FBgn0265149', 'FBgn0262252', 'FBgn0031235', 'FBgn0263465'], ['FBgn0022246', 'FBgn0031238', 'FBgn0031239', 'FBgn0265150']), 'FBtr0309810': (['FBgn0263584', 'FBgn0031209', 'FBgn0002121', 'FBgn0031208'], ['FBgn0051973', 'FBgn0267987', 'FBgn0266878', 'FBgn0266879']), 'FBtr0345733': (['FBgn0266879', 'FBgn0266878', 'FBgn0267987', 'FBgn0051973'], ['FBgn0067779', 'FBgn0266322', 'FBgn0031213', 'FBgn0031214']), 'FBtr0345732': (['FBgn0266878', 'FBgn0267987', 'FBgn0051973', 'FBgn0263584'], ['FBgn0266879', 'FBgn0067779', 'FBgn0266322', 'FBgn0031213']), 'FBtr0344053': (['FBgn0266322', 'FBgn0067779', 'FBgn0266879', 'FBgn0266878'], ['FBgn0031213', 'FBgn0031214', 'FBgn0002931', 'FBgn0031216']), 'FBtr0344052': (['FBgn0266322', 'FBgn0067779', 'FBgn0266879', 'FBgn0266878'], ['FBgn0031213', 'FBgn0031214', 'FBgn0002931', 'FBgn0031216']), 'FBtr0347585': (['FBgn0267987', 'FBgn0051973', 'FBgn0263584', 'FBgn0031209'], ['FBgn0266878', 'FBgn0266879', 'FBgn0067779', 'FBgn0266322']), 'FBtr0347595': (['FBgn0267996', 'FBgn0267995', 'FBgn0031245', 'FBgn0031244'], ['FBgn0025686', 'FBgn0031247', 'FBgn0017457']), 'FBtr0347594': (['FBgn0267995', 'FBgn0031245', 'FBgn0031244', 'FBgn0003444'], ['FBgn0267996', 'FBgn0025686', 'FBgn0031247', 'FBgn0017457']), 'FBtr0336988': (['FBgn0265150', 'FBgn0031239', 'FBgn0031238', 'FBgn0022246'], ['FBgn0031240', 'FBgn0086912', 'FBgn0086856', 'FBgn0086855']), 'FBtr0336984': (['FBgn0265151', 'FBgn0266557', 'FBgn0031233', 'FBgn0031232'], ['FBgn0265153', 'FBgn0265152', 'FBgn0263465', 'FBgn0031235']), 'FBtr0336985': (['FBgn0265153', 'FBgn0265151', 'FBgn0266557', 'FBgn0031233'], ['FBgn0265152', 'FBgn0263465', 'FBgn0031235', 'FBgn0262252']), 'FBtr0336986': (['FBgn0265152', 'FBgn0265153', 'FBgn0265151', 'FBgn0266557'], ['FBgn0263465', 'FBgn0031235', 'FBgn0262252', 'FBgn0265149']), 'FBtr0344032': (['FBgn0266304', 'FBgn0025683', 'FBgn0031220', 'FBgn0031219'], ['FBgn0001142', 'FBgn0265074', 'FBgn0265075', 'FBgn0051975']), 'FBtr0336836': (['FBgn0265074', 'FBgn0001142', 'FBgn0266304', 'FBgn0025683'], ['FBgn0265075', 'FBgn0051975', 'FBgn0051976', 'FBgn0051974']), 'FBtr0336837': (['FBgn0265075', 'FBgn0265074', 'FBgn0001142', 'FBgn0266304'], ['FBgn0051975', 'FBgn0051976', 'FBgn0051974', 'FBgn0031224'])}
    #ncRNA_gff_dict
    #{'FBtr0336987': ('2L', '248909', '249812'), 'FBtr0309810': ('2L', '21952', '24237'), 'FBtr0345733': ('2L', '66318', '66524'), 'FBtr0345732': ('2L', '65999', '66242'), 'FBtr0344053': ('2L', '71039', '73836'), 'FBtr0344052': ('2L', '71039', '73836'), 'FBtr0347585': ('2L', '54817', '55767'), 'FBtr0347595': ('2L', '286383', '288292'), 'FBtr0347594': ('2L', '283807', '284255'), 'FBtr0336988': ('2L', '262203', '263003'), 'FBtr0336984': ('2L', '219651', '220050'), 'FBtr0336985': ('2L', '223519', '224708'), 'FBtr0336986': ('2L', '228037', '228588'), 'FBtr0344032': ('2L', '122837', '124031'), 'FBtr0336836': ('2L', '135402', '136128'), 'FBtr0336837': ('2L', '136144', '136701')}
#this works!
def mel_gene_set(dict):
    """This function finds unique mel genes, and puts them in a set (what is returned), so we don't get the same coords twice. It takes fbgn_id_dict"""
    mel_gene_set = set()
    for k, v in dict.iteritems():
        for mg in v[0]:
            mel_gene_set.add(mg)
        for mg in v[1]:
            mel_gene_set.add(mg)
    return mel_gene_set
                                                          
def ortho_mapping(set):
    """This function maps other another species' orthologs to dmel genes"""
    mapping = dict()
    with open(sys.argv[2], 'r') as orthos:
        for line in orthos:
            if not line.startswith('#') and not line.startswith('\n'):
                data = line.split ('\t')
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

gff_obj =gff_list()
ncRNA_obj = ncRNA_list(gff_obj)
#print ncRNA_obj
fbgn_dict_obj, ncRNA_gff_dict_obj =ncRNA_gff_dict(ncRNA_obj, gff_obj, window_length)
#print fbgn_dict_obj
#print ncRNA_gff_dict_obj
mel_genes_obj = mel_gene_set(fbgn_dict_obj)
#print mel_genes_obj
ortho_map = ortho_mapping(mel_genes_obj)
#print ortho_map
rna_ortho_dict = mel_rna_ortho(fbgn_dict_obj, ortho_map)
print rna_ortho_dict
#python lncRNA_annotation_doc/step_by_step_flanking.py out/1_prac_ncRNA_and_genes.txt 

