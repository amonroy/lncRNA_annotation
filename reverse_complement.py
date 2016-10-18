


def reverse_complement(string):
	"""This function returns the reverse complement of a string of nucleotides"""
	print "This is your starting string:"
	print string
	reverse = string[::-1]
	print "This is the reverse of your string"
	print reverse
	complement = []
	for i in reverse:
		if i == 'a':
			complement.append('t')
		if i == 't':
			complement.append('a')
		if i == 'c':
			complement.append('g')
		if i == 'g':
			complement.append('c')
	print "This is the reverse complement"
	print complement
	rev_com_str = ''.join(complement)
	print "This is the string rev com"
	print rev_com_str