# FPSAC: Fast Phylogenetic Scaffolding of Ancient Contigs.
# May 2013.
# Cedric Chauve: cedric.chauve@sfu.ca

# Generate a file with the extremities of a set of CARs

# argument 1: Input  Scaffold order doubled
# argument 2: Output CARs extremities, one per line

import sys

# command_line="python "
# for a in sys.argv:
#     command_line=command_line+' '+a
# print "------> COMMAND: "+command_line

scaffold_order_file=open(sys.argv[1],"r").readlines()
output_file=open(sys.argv[2],"w")
for s in scaffold_order_file:
    if s[0]=="#":
        car_name=s.rstrip('\n')
    elif s[0] !=">":
        s1=s.split(' ')
        ext1=s1[1]
        ext2=s1[-2]
        output_file.write(ext1+"\t"+car_name+"\n")
        output_file.write(ext2+"\t"+car_name+"\n")
