# FPSAC: Fast Phylogenetic Scaffolding of Ancient Contigs.
# May 2013.
# Cedric Chauve: cedric.chauve@sfu.ca

# Computing a scaffold from adjacencies and repeat spanning intervals

import sys

# argument 1 = repeat spanning intervals
# argument 2 = adjacencies
# argument 3 = copy numbers
# argument 4 = scaffolds output file (doubled markers)
# argument 5 = unassigned adjacencies output file (doubled markers)
# argument 6 = ancestor name

# command_line="python "
# for a in sys.argv:
#     command_line=command_line+' '+a
# print "------> COMMAND: "+command_line

#------ READING INPUT ------------------------------------------------------------
# ----- Globl variables ----------------------------------------------------------
copynumbers={}          # Copy numbers of doubled markers
markers={}              # Markers adjacencies markers[2x]=2x-1, markers[2x-1]=2x
undoubled_markers = {}  # undoubled_marker[x]=y <=> x=2y or 2y-1
marker_weight = "10000" # Weight of marker adjacencies

available_linkages={} # Indexed by markers
# Each non-repeat marker is extremity of at most one linkage
# at any time. Repeats can be extremities of several linkages, i.e. a list
# A non-repeat that is extremity of an interval is also extremity of an adjacency
# that is not recorded. By construction, if a marker is extremity of an interval, 
# it appears first in the list.
# Also important, extremities of scaffold are half markers, i.e. not marker adjacencies.

# ----- Reading markers ----------------------------------------------------------
profilesfile = open(sys.argv[3],"r").readlines()
i=0
while i<len(profilesfile):
    undoubled_marker=int(profilesfile[i].split('>')[1].split(' ')[0])    
    markers[2*undoubled_marker]=2*undoubled_marker-1
    markers[2*undoubled_marker-1]=2*undoubled_marker
    undoubled_markers[2*undoubled_marker]=undoubled_marker
    undoubled_markers[2*undoubled_marker-1]=undoubled_marker
    # max id of a marker
    max_marker=2*undoubled_marker
    i=i+1

# ----- Initialisations --------------------------------------------------------
for marker in range(max_marker+1):
    available_linkages[marker]=[]
    copynumbers[marker]=0

# ----- Reading copy numbers --------------------------------------------------------
i=0
while i<len(profilesfile):
    undoubled_marker=int(profilesfile[i].split('>')[1].split(' ')[0])
    copynumber=int(profilesfile[i].split('>')[1].split(' ')[1])
    copynumbers[2*undoubled_marker]=copynumber
    copynumbers[2*undoubled_marker-1]=copynumber
    i=i+1

# ----- Reading intervals ----------------------------------------------------------
intervalsfile = open(sys.argv[1],"r").readlines()
used_adjacencies=[]
i=0
while i<len(intervalsfile):
    id=intervalsfile[i].split('|')[0]
    interval1=intervalsfile[i].split(':')[1]
    interval_gaps=interval1[interval1.find(' #'):].rstrip()
    interval2=interval1[:interval1.find(' #')].split(' ')
    interval3=[]
    # Set interval3=[a r1 r2 ... rk b]
    for j in range(0,len(interval2)):
        interval3.append(int(interval2[j].rstrip('\n')))
    interval4=[]
    # Set interval4=[a r1 r2 ... rk b].reverse()
    for j in range(0,len(interval2)):
        interval4.append(int(interval2[j].rstrip('\n')))
    interval4.reverse()
    # Adding the intervals in the linkages list of both extremities
    marker1=int(interval2[0])
    marker2=int(interval2[len(interval2)-1].rstrip('\n'))
    available_linkages[marker1].append(interval3)
    available_linkages[marker2].append(interval4)
    # Recording adjacencies in intervals
    for j in range(1,len(interval3)):
        used_adjacencies.append([int(interval2[j-1]),int(interval2[j])])
        used_adjacencies.append([int(interval2[j]),int(interval2[j-1])])
    i=i+1

# ----- Reading adjacencies ----------------------------------------------------------
adjacenciesfile = open(sys.argv[2],"r").readlines()
i=0
while i<len(adjacenciesfile):
    weight=adjacenciesfile[i].split('|')[1].split(';')[0]
    marker1=int(adjacenciesfile[i].split(':')[1].split(' ')[0])
    marker2=int(adjacenciesfile[i].split(':')[1].split(' ')[1])
    if weight != marker_weight: # Not a marker adjacency
        adjacency1=[marker1,marker2]
        adjacency2=[marker2,marker1]
        if not adjacency1 in used_adjacencies:
            available_linkages[marker1].append(adjacency1)
        if not adjacency2 in used_adjacencies:
            available_linkages[marker2].append(adjacency2)
    i=i+1

# ----- COMPUTING SCAFFOLDS ----------------------------------------------------------
# ----- Algorithm --------------------------------------------------------------------
# Assumption: a scaffold is framed by two markers with left occurrences extactly 1
# Initialization: create an array of scaffolds defined by intervals and non-repeat 
# adjacencies
# For each marker
#   Check if it creates a circular scaffold in which case nothing is done
#   Check if it fuses two scaffolds in which case do it
# For each scaffold
#   Complete its extremities
#   

left_occurrences={}      # Left occurrences of markers that can be used for each marker
for marker in range(max_marker+1):
    left_occurrences[marker]=len(available_linkages[marker])
    # There is no point at using a marker more than the number of times it is extremity 
    # of an adj/interval

scaffold_array={}       # Array of scaffolds indexed by id
scaffold_id=0           # Max scaffold index
scaffold_extremities={} # Link extremities, scaffold id

def print_scaffolds(text):
    i=0
    for marker in range(max_marker+1):
        if scaffold_extremities[marker]!=-1:
            id=scaffold_extremities[marker]
            scaffold=scaffold_array[id]
            ext1=scaffold[0]
            ext2=scaffold[len(scaffold)-1]
            if ext1==marker:
                print str(i)+"\t"+text+"\t"+str(scaffold)
                i+=1

# Initialization of the list of scaffolds with the ones framed by non-repeats
# i.e. markers having a single left occurrence
for marker in range(max_marker+1):
    scaffold_extremities[marker]=-1
for marker in range(max_marker+1):
    linkage=available_linkages[marker]
    if len(linkage)==1:
        ext1=linkage[0][0]
        ext2=linkage[0][len(linkage[0])-1]
        if left_occurrences[ext1]==1 and left_occurrences[ext2]==1:
            scaffold_array[scaffold_id]=linkage[0]
            scaffold_extremities[ext1]=scaffold_id
            scaffold_extremities[ext2]=scaffold_id
            scaffold_id+=1
            available_linkages[marker]=[]

# Computing the type of a scaffold
CIRC=1
LIN=2
def check_scaffold_type(marker):
    scaffold=scaffold_array[scaffold_extremities[marker]]
    ext1=scaffold[0]
    ext2=scaffold[len(scaffold)-1]
    if markers[ext1]==ext2:
        return CIRC
    else:
        return LIN

# Extension of scaffolds: stage 1, fusion of scaffolds framed by non-repeats
for marker in range(max_marker+1):
    if left_occurrences[marker]==1 and scaffold_extremities[marker]!=-1:
        if check_scaffold_type(marker)!=CIRC:
            ext1=marker
            ext2=markers[marker]
            id1=scaffold_extremities[ext1]
            id2=scaffold_extremities[ext2]
            if id1!=-1 and id2!=-1:
                # Preparing the scaffolds
                if ext1==scaffold_array[id1][0]:
                    scaffold_array[id1].reverse()
                if ext2!=scaffold_array[id2][0]:
                    scaffold_array[id2].reverse()
                # Getting extremities
                ext1a=scaffold_array[id1][0]
                ext1b=scaffold_array[id1][len(scaffold_array[id1])-1]
                ext2a=scaffold_array[id2][0]
                ext2b=scaffold_array[id2][len(scaffold_array[id2])-1]
                # Concatenating scaffolds
                scaffold_array[id1]=scaffold_array[id1]+scaffold_array[id2]
                scaffold_array[id2]=[]
                # Updating extremities
                scaffold_extremities[ext1a]=id1
                scaffold_extremities[ext2b]=id1
                scaffold_extremities[ext1b]=-1
                scaffold_extremities[ext2a]=-1
                # Updating left occurrences
                left_occurrences[ext1b]-=1
                left_occurrences[ext2a]-=1


# Extension of scaffolds: stage 2, extending scaffold extremities
scaffold_extended={}   # records the number of times each scaffold has been extended
for marker in range(max_marker+1):
     if scaffold_extremities[marker]!=-1:
         if check_scaffold_type(marker)!=CIRC:
             id=scaffold_extremities[marker]
             scaffold_extended[id]=False

for marker in range(max_marker+1):
     if scaffold_extremities[marker]!=-1:
         if check_scaffold_type(marker)!=CIRC:
             id=scaffold_extremities[marker]
             if scaffold_extended[id]==False:
                 scaffold=scaffold_array[id]
                 ext2=scaffold[len(scaffold)-1]
                 next=markers[ext2]
                 scaffold.append(next)
                 scaffold_extremities[ext2]=-1
                 scaffold_extremities[next]=id
                 if len(available_linkages[next])==1:
                     adjacency=available_linkages[next][0]
                     scaffold.append(adjacency[1])
                     scaffold.append(markers[adjacency[1]])
                     scaffold_extremities[markers[adjacency[1]]]=id
                     scaffold_extremities[next]=-1
                     available_linkages[next]=[]
                 scaffold.reverse()
                 ext2=scaffold[len(scaffold)-1]
                 next=markers[ext2]
                 scaffold.append(next)
                 scaffold_extremities[ext2]=-1
                 scaffold_extremities[next]=id
                 if len(available_linkages[next])==1:
                     adjacency=available_linkages[next][0]
                     scaffold.append(adjacency[1])
                     scaffold.append(markers[adjacency[1]])
                     scaffold_extremities[markers[adjacency[1]]]=id
                     scaffold_extremities[next]=-1
                     available_linkages[next]=[]
                 scaffold_extended[id]=True

# ----- Writing SCAFFOLDS ------------------------------------------------------------
scaffoldfile = open(sys.argv[4],"w")
scaffoldfile.write(">"+sys.argv[6]+"\n")
def open_car(opening,id):
    scaffoldfile.write("#CAR"+str(id)+"\n"+opening)
def close_car(closing):
    scaffoldfile.write(closing+"\n")
i=1
for marker in range(max_marker+1):
    if scaffold_extremities[marker]!=-1:
        id=scaffold_extremities[marker]
        scaffold=scaffold_array[id]
        ext1=scaffold[0]
        ext2=scaffold[len(scaffold)-1]
        if ext1==marker:
            if check_scaffold_type(marker)==CIRC:
                open_car("_C ",i)
                # Correcting the circular scaffold
                scaffold_extremities[scaffold[len(scaffold)-1]]=-1
                scaffold.remove(ext1)           
                scaffold.append(ext1)
                scaffold_extremities[scaffold[0]]=id
                for m in scaffold:
                    scaffoldfile.write(str(m)+" ")
                close_car("C_")
            else:
                open_car("_Q ",i)
                for m in scaffold:
                    scaffoldfile.write(str(m)+" ")
                close_car("Q_")
            i+=1

# ----- Writing SCAFFOLDS ------------------------------------------------------------
unassignedfile = open(sys.argv[5],"w")
id=0
unassigned_weight=2*float(marker_weight)
for marker in range(max_marker+1):
    if available_linkages[marker]!=[]:
        for adjacency in available_linkages[marker]:
            unassignedfile.write(str(id)+"|"+str(unassigned_weight)+";unassigned_edge:"+str(adjacency[0])+" "+str(adjacency[1])+" #\n")
            id+=1
            
