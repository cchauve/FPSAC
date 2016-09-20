# FPSAC: Fast Phylogenetic Scaffolding of Ancient Contigs.
# May 2013.
# Cedric Chauve: cedric.chauve@sfu.ca

# Greedy filtering of repeat spanning intervals using copy numbers and
# kept adjacencies

# argument 1 = Input  repeat spanning intervals file
# argument 2 = Input  kept adjacencies file
# argument 3 = Input  ancestral content file
# argument 4 = Output kept intervals file
# argument 5 = Output discarded intervals file

import sys
from operator import itemgetter

# command_line="python "
# for a in sys.argv:
#     command_line=command_line+' '+a
# print "------> COMMAND: "+command_line

# ----- Reading adjacencies ----------------------------------------------------------
adjacencieslist = []
adjacenciesfile = open(sys.argv[2],"r").readlines()
i=0
while i<len(adjacenciesfile):
    adjacency1=adjacenciesfile[i].split(':')[1].split(' ')[0]
    adjacency2=adjacenciesfile[i].split(':')[1].split(' ')[1]
    adjacencieslist.append((adjacency1,adjacency2))
    adjacencieslist.append((adjacency2,adjacency1))
    i=i+1

# ----- Reading profiles ----------------------------------------------------------
copynumbers={}
profilesfile = open(sys.argv[3],"r").readlines()
i=0
while i<len(profilesfile):
    undoubled_marker=int(profilesfile[i].split('>')[1].split(' ')[0])
    copynumber=int(profilesfile[i].split('>')[1].split(' ')[1])
    copynumbers[2*undoubled_marker]=copynumber
    copynumbers[2*undoubled_marker-1]=copynumber
    i=i+1

# ----- Reading intervals ----------------------------------------------------------
intervalslist = []
intervalsfile = open(sys.argv[1],"r").readlines()
i=0
while i<len(intervalsfile):
    id=int(intervalsfile[i].split('|')[0])
    weight=float(intervalsfile[i].split('|')[1].split(';')[0])
    species=intervalsfile[i].split('|')[1].split(';')[1].split(':')[0]
    interval=intervalsfile[i].split(':')[1].split('#')[0].rstrip(' ')
    gaps1=intervalsfile[i].rstrip( ).split(' #')
    gaps='#'+' #'.join(gaps1[1:len(gaps1)])
    intervalslist.append((id,weight,species,interval,gaps))
    i=i+1

intervalslist.sort(key=itemgetter(1))
intervalslist.reverse()

# ------ Auxiliary function --------------------------------------------------------
# Tests if two doubled markers originate from the same undoubled one
# Assumption: marker 1 < marker2
def test_marker_adj(marker1, marker2):
    if marker1%2!=0 and marker2==marker1+1:
        return True
    else:
        return False

# ----- Filtering intervals ----------------------------------------------------------
kept_intervals=[]
discarded_intervals=[]
trackcopynumbers={}  # records used copies of each marker
for i in copynumbers:
    trackcopynumbers[i]=0
    
for interval in intervalslist:
    markers=interval[3].rstrip('\n').split(' ')
    keep=1
    for i in range(len(markers)):
        if i<len(markers)-1:
            if int(markers[i]) < int(markers[i+1]):
                marker1=int(markers[i])
                marker2=int(markers[i+1])
            else:
                marker1=int(markers[i+1])
                marker2=int(markers[i])
                
            if (not (markers[i],markers[i+1]) in adjacencieslist) and (not test_marker_adj(marker1,marker2)):
                keep=0
        trackcopynumbers[int(markers[i])]=trackcopynumbers[int(markers[i])]+1
        if trackcopynumbers[int(markers[i])] > copynumbers[int(markers[i])]:
            keep=0

    if keep==1:
        kept_intervals.append(interval)
    else:
        discarded_intervals.append(interval)
        for i in range(len(markers)):
            trackcopynumbers[int(markers[i].rstrip('\n'))]=trackcopynumbers[int(markers[i].rstrip('\n'))]-1

kept_intervals.sort(key=itemgetter(0))
discarded_intervals.sort(key=itemgetter(0))

# ----- Writing output ----------------------------------------------------------
outputkeptfile=open(sys.argv[4],"w")
for interval in kept_intervals:
    id=int(interval[0])
    weight=interval[1]
    species=interval[2]
    sequence=interval[3]
    gaps=interval[4]
    outputkeptfile.write(str(id)+"|"+str(weight)+";"+species+":"+sequence+' '+gaps+'\n')
outputkeptfile.close()

outputdiscardedfile=open(sys.argv[5],"w")
for interval in discarded_intervals:
    id=interval[0]
    weight=interval[1]
    species=interval[2]
    sequence=interval[3]
    gaps=interval[4]
    outputdiscardedfile.write(str(id)+"|"+str(weight)+";"+species+":"+sequence+' '+gaps+'\n')
outputdiscardedfile.close()
