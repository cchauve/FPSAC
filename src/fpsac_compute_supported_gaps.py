# FPSAC: Fast Phylogenetic Scaffolding of Ancient Contigs.
# May 2013.
# Cedric Chauve: cedric.chauve@sfu.ca

# Selecting gaps with supported length, following a Dollo criterion

# argument 1 = Input  gaps file
# argument 2 = Input  species pairs
# argument 3 = Input  0/1 (0=unsupported gaps not present, 1=all gaps added for unsupported gaps)
# argument 4 = Input  ratio for extending/shrinking gaps length defining support (float, between 0 and 1)
# argument 5 = Output selected gaps

import sys
from operator import itemgetter

# command_line="python "
# for a in sys.argv:
#     command_line=command_line+' '+a
# print "------> COMMAND: "+command_line

# ----- Reading species pairs ---------------------------------------------------------
speciespairslist=[]
speciespairsfile=open(sys.argv[2],"r").readlines()
for l in speciespairsfile:
    if l[0] != "#":
        species1=l.split(' ')[0]
        species2=l.split(' ')[1].rstrip('\n')
        speciespairslist.append((species1,species2))
        speciespairslist.append((species2,species1))

# ----- Reading parameters ------------------------------------------------------------
par_complement=int(sys.argv[3])
par_stretching=float(sys.argv[4])

# ----- Reading gaps ------------------------------------------------------------------
gapslist={}
gapsfile=open(sys.argv[1],"r").readlines()
outputfile=open(sys.argv[5],"w")
for l in gapsfile:
    if l[0]==">":
        gap_header=l.rstrip('\n')
        gap_min=-1
        gap_max=-1
        gap_supported_min=-1
        gap_supported_max=-1
        gap_index=0
    elif len(l)>1:
        gapslist[gap_index]=l
        gap_index+=1
    else:
        for i in range(gap_index):
            g1=gapslist[i]
            species1=g1.split('.')[0]
            start1=int(g1.split(' ')[0].split(':')[1].split('-')[0])
            end1=int(g1.split(' ')[0].split(':')[1].split('-')[1])
            lg1=float(end1-start1+1)
            if lg1>=0 and (lg1<gap_min or gap_min==-1):
                gap_min=int(lg1)
            if lg1>gap_max:
                gap_max=int(lg1)
            for j in range(i+1,gap_index):
                g2=gapslist[j]
                species2=g2.split('.')[0]
                start2=int(g2.split(' ')[0].split(':')[1].split('-')[0])
                end2=int(g2.split(' ')[0].split(':')[1].split('-')[1])
                lg2=float(end2-start2+1)
                if (species1,species2) in speciespairslist:
                    if (lg2>=lg1*(1-par_stretching) and lg2<=lg1*(1+par_stretching)) or (lg1>=lg2*(1-par_stretching) and lg1<=lg2*(1+par_stretching)):
                        if min(lg1,lg2)>=0 and (min(lg1,lg2)<gap_supported_min or gap_supported_min==-1):
                            gap_supported_min=int(min(lg1,lg2))
                        if max(lg1,lg2)>gap_supported_max:
                            gap_supported_max=int(max(lg1,lg2))
        if gap_supported_min>-1:
            outputfile.write(gap_header+" supported_length_interval "+str(gap_supported_min)+"-"+str(gap_supported_max)+"\n")
        elif par_complement==1:
            outputfile.write(gap_header+" unsupported_length_interval "+str(gap_min)+"-"+str(gap_max)+"\n")
        for i in range(gap_index):
            g=gapslist[i]
            species=g.split('.')[0]
            start=int(g.split(' ')[0].split(':')[1].split('-')[0])
            end=int(g.split(' ')[0].split(':')[1].split('-')[1])
            lg=end-start+1
            if (lg>=gap_supported_min and lg<=gap_supported_max) or (par_complement==1):
                outputfile.write(g)
        outputfile.write("\n")



