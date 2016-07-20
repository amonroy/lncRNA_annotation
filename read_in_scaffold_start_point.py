#!/usr/bin/python

import sys
#import re
sequence = ""
found = False

file_to_split = sys.argv[1]
scaffold = raw_input("What's the scaffold?")
f = open(file_to_split, 'r')

for line in f:
    if found:
        #print line
        if line.startswith('>'):
            break
        sequence += line.strip()
        
    else:
        if line.strip().startswith('>'+scaffold):
            found = True
            print "This is the line I wanted!!!!!!!!!!!!\n", line
print "\n"
print "*"*60

print "This is sequence:\n", sequence
f.close()