###practice on HSR-omega, which is located on 3R
### this script is designed to pull out specific starts and stop on chromosome and the 
### sequence in between

import sys

#"~/Desktop/Projects/RNA_annotation/FlyBase/Dyak/dyak_3R_prac.fa"

in_file = sys.argv[1]

end_start = 1
end_stop = 50
start = end_start-1
stop = end_stop+1

with open (in_file, 'r') as i: #open ('3R_fafrag_prac.fa', 'w') as out:
	header = i.readline().rstrip()
	print header
	for line in i:
		#header = line
		#print header
		sequence = line.strip('\n')
		#print sequence
		print sequence[end_start:end_stop]
			
		#print sequence[start:stop]
	

	