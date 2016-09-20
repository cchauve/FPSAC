# FPSAC: Fast Phylogenetic Scaffolding of Ancient Contigs.
# May 2013.
# Cedric Chauve: cedric.chauve@sfu.ca

# Correcting the fact that markers order is lost in repeat spanning intervals during weighting
# Minor issue to address in the next release.
  
# argument 1 = Input  weighted repeat spanning intervals
# argument 2 = Input  repeat spanning intervals
# argument 3 = Output file

import sys

# command_line="python "
# for a in sys.argv:
#     command_line=command_line+' '+a
# print "------> COMMAND: "+command_line

intervalsfile1 = open(sys.argv[1],"r").readlines()
intervalsfile2 = open(sys.argv[2],"r").readlines()
outputfile=open(sys.argv[3],"w")

weightslist = []
i=0
while i<len(intervalsfile1):
    weight=intervalsfile1[i].split('|')[1].split(';')[0]
    weightslist.append(weight)
    i=i+1
i=0
while i<len(intervalsfile2):
    id=intervalsfile2[i].split('|')[0]
    interval=intervalsfile2[i].split(';')[1]
    outputfile.write(id+"|"+weightslist[i]+";"+interval)
    i=i+1

