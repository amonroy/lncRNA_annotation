## this is me going through and looking closely at 2_flanking_genes_ncRNA_2_v1.py

import sys
window_length = 4 #how many upstream genes and downstream genes

def gff_list():
    """This function takes the modified gff3 file and creates a list"""
    mod_gff3 = sys.argv[1]
    with open(mod_gff3, 'r') as f:
        gff = [line.strip().split('\t') for line in f]
        f.close()
    return gff
        
        
        #gff_list ex/:
	    #[['2L', 'FlyBase', 'gene', '7529', '9484', '.', '+', '.', 'ID=FBgn0031208;Name=CG11023;Ontology_term=SO:0000010,SO:0000087,GO:0016929,GO:0016926;Dbxref=FlyBase:FBan0011023,FlyBase_Annotation_IDs:CG11023,GB_protein:ACZ94128,GB_protein:AAO41164,GB:AI944728,GB:AJ564667,GB_protein:CAD92822,GB:BF495604,UniProt/TrEMBL:Q86BM6,INTERPRO:IPR003653,GB_protein:AGB92323,UniProt/TrEMBL:M9PAY1,OrthoDB7_Drosophila:EOG796K1P,OrthoDB7_Diptera:EOG7X1604,EntrezGene:33155,UniProt/TrEMBL:E1JHP8,UniProt/TrEMBL:Q6KEV3,OrthoDB7_Insecta:EOG7Q8QM7,OrthoDB7_Arthropoda:EOG7R5K68,OrthoDB7_Metazoa:EOG7D59MP,InterologFinder:33155,BIOGRID:59420,FlyAtlas:CG11023-RA,GenomeRNAi:33155;gbunit=AE014134;derived_computed_cyto=21A5-21A5'], 
	    # ['2L', 'FlyBase', 'gene', '9839', '21376', '.', '-', '.', 'ID=FBgn0002121;Name=l(2)gl;fullname=lethal 

def ncRNA_list(list):
    """This function takes the gff_list and makes an ncRNA_list"""
    ncRNA = [] #initiates list
    for i in gff_list:
        if i[2] == 'ncRNA':
            preidRNA = i[8].split(';')[0]
            #[ID=FBgn0031208];Name=CG11023;Ontology_term=SO:0000010,SO:0000087,GO:0016929,GO:0016926;Dbxref=FlyBase:FBan0011023,FlyBase_Annotation_IDs:CG11023,GB_protein:ACZ94128,GB_protein:AAO41164,GB:AI944728,GB:AJ564667,GB_protein:CAD92822,GB:BF495604,UniProt/TrEMBL:Q86BM6,INTERPRO:IPR003653,GB_protein:AGB92323,UniProt/TrEMBL:M9PAY1,OrthoDB7_Drosophila:EOG796K1P,OrthoDB7_Diptera:EOG7X1604,EntrezGene:33155,UniProt/TrEMBL:E1JHP8,UniProt/TrEMBL:Q6KEV3,OrthoDB7_Insecta:EOG7Q8QM7,OrthoDB7_Arthropoda:EOG7R5K68,OrthoDB7_Metazoa:EOG7D59MP,InterologFinder:33155,BIOGRID:59420,FlyAtlas:CG11023-RA,GenomeRNAi:33155;gbunit=AE014134;derived_computed_cyto=21A5-21A5'
            ncRNA.append(preidRNA)
    return ncRNA
            
####This part is by far the most confusing! Go back through and comment it up            
def ncRNA_gff_dict(rna_list, gff_list, window_length):
    """This function takes our two lists, ncRNA and gff_list and makes a great dictionary"""
    ncRNA_gff_dict = {}
    for r in rna_list:
        for i in gff_list:
            data = i[8].split(';')[0]
            if data == r:
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
        return fbgn_id_dict
        return ncRNA_gff_dict
                                                          
gff_list()
ncRNA_list(gff_list)
ncRNA_gff_dict(ncRNA, gff, window_length)