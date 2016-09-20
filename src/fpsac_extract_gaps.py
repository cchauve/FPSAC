# FPSAC: Fast Phylogenetic Scaffolding of Ancient Contigs.
# May 2013.
# Cedric Chauve: cedric.chauve@sfu.ca

# Extracting gaps coordinates from kept adjacencies and repeat spanning intervals

# argument 1 = Input  kept adjacencies file
# argument 2 = Input  kept repeat spanning intervals file
# argument 3 = Input  ancestral content file
# argument 4 = Output gaps

import sys
from operator import itemgetter

# command_line="python "
# for a in sys.argv:
#     command_line=command_line+' '+a
# print "------> COMMAND: "+command_line

# ------ Auxiliary function --------------------------------------------------------
# Tests if two doubled markers originate from the same undoubled one
# Assumption: marker 1 < marker2
def test_marker_adj(marker1, marker2):
    if marker1%2!=0 and marker2==marker1+1:
        return True
    else:
        return False            

# ----- Reading profiles ----------------------------------------------------------
copynumbers={}
profilesfile = open(sys.argv[3],"r").readlines()
for l in profilesfile:
    undoubled_marker=int(l.split('>')[1].split(' ')[0])
    copynumber=int(l.split('>')[1].split(' ')[1])
    copynumbers[2*undoubled_marker]=copynumber
    copynumbers[2*undoubled_marker-1]=copynumber

# ----- Reading adjacencies ----------------------------------------------------------
gap_id=1
output=open(sys.argv[4],"w")
adjacenciesfile = open(sys.argv[1],"r").readlines()
for l in adjacenciesfile:
    a1=int(l.split(':')[1].split(' ')[0])
    a2=int(l.split(':')[1].split(' ')[1]) # a2 assumed to be lower than a1
    gap=l.rstrip('\n').split('#')[1:]
    c1=copynumbers[a1]
    c2=copynumbers[a2]
#    if c1==1 and c2==1 and not test_marker_adj(min(a1,a2),max(a1,a2)): # no repeat
    if c1==1 or c2==1 and not test_marker_adj(min(a1,a2),max(a1,a2)): # at least one non-repeat
# we chose solution 2 because if the non-repeat repeat adjacency is present sometimes outside of a repeat spanning interval
# we will miss some of its extant occurrences
# problem: if we have a configuration b a b where b is a repeat and a a non-repeat, we see two extant occurrences
        output.write(">gap_"+str(gap_id)+" adjacency "+str(a1)+"-"+str(a2)+"\n")
        for g in gap:
            orientation=g.split(' ')[1].rstrip( ) # orientation of the adjacenc
            # + means it was seen as a1 a2 in the genome
            # - means it was seen as a2 a1 in the genome
            sgn="+" # indicates on which strand to get the gap sequence
            # rule: gap is taken on on the same strand than a1, the smaller marker
            if (orientation=="+" and a1%2==1) or (orientation=="-" and a1%2==0):
                sgn="-"
            output.write(g.split(' ')[0]+" "+sgn+"\n")
        output.write("\n")
        gap_id+=1

# ----- Reading intervals ----------------------------------------------------------
intervalsfile = open(sys.argv[2],"r").readlines()
for l in intervalsfile:
    id=int(l.split('|')[0])
    weight=float(l.split('|')[1].split(';')[0])
    species=l.split('|')[1].split(';')[1].split(':')[0]
    interval=l.split(':')[1].split('#')[0].rstrip(' ').split(' ')
    gap1=l.rstrip('\n').split('#')[1:]
    species={}
    gap2={}
    orientations={}
    for i in range(len(gap1)):
        species[i]=gap1[i].split(':')[0]
        orientations[i]=gap1[i].split(' ')[1]
        gap2[i]=gap1[i].split(':')[1].rstrip( ).split(' ')[0].split(',')
        if orientations[i]=="-":
            gap2[i].reverse()
    str_interval=' '.join(interval)
    for i in range(len(interval)-1):
        a1=int(interval[i])
        a2=int(interval[i+1])
        c1=copynumbers[a1]
        c2=copynumbers[a2]
        amin=min(a1,a2)
        amax=max(a1,a2)
#        if i%2==0: # avoiding gaps between the two extremities of a given marker
        if i%2==0 and c1>1 and c2>1: # avoiding gaps between the two extremities of a given marker and involving a non-repeat
            output.write(">gap_"+str(gap_id)+" adjacency "+str(amin)+"-"+str(amax)+" pos "+str(i+1)+" in repeat_spanning_interval "+str_interval+"\n")
            # adjacency written with smaller marker first, to be consistent with non-interval adjacencies
            for g in range(len(gap1)):
                sgn="+" # same rule than for adjacencies: gap sequence taken on the strand of the smaller marker
                if (orientations[g]=="+"): # interval written as seen in the genome
                    if (amin==a1 and amin%2==1) or (amin==a2 and amin%2==0): 
                        sgn="-"
                else: # interval seen reversed in the genome
                    if (amin==a2 and amin%2==1) or (amin==a1 and amin%2==0): 
                        sgn="-"
                output.write(species[g]+":"+gap2[g][i]+" "+sgn+"\n")
            output.write("\n")
            gap_id+=1
