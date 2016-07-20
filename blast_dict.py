####The purpose of this script is to double check the other script (flanking, yadda yadda) is working.
#### I am just looking out the output of blastn non-mel fly against the annotated dmel ncRNA file. Those files are in the "out" folder

import sys

def blast_coord():
	"""This function will find all of the hits of dmel ncRNA against non-mel genome and return the best scaffold and the start and stop"""
	file = sys.argv[1]
	blast_coord_dict = dict()
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
		print k, uscafs, min(sstart), max(send)
			#want small, large
			#print type(list[7])
			#print type(list[8])
			#if max(list[7]) > min(list[8]):
			#	print '*'+k, max(list[8]), min(list[7])
			#if max(list[7]) < min(list[8]):
			#	print k, max(list[7]), min(list[8])

	return blast_coord_dict
	
bls_dict =blast_coord()
#print bls_dict
