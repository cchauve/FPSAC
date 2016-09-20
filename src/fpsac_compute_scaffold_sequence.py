# FPSAC: Fast Phylogenetic Scaffolding of Ancient Contigs.
# May 2013.
# Cedric Chauve: cedric.chauve@sfu.ca

# Computing a scaffold sequence, including gap sequences for gaps

# argument 1 = Input  families FASTA file
# argument 2 = Input  scaffold (undoubled)
# argument 3 = Input  gaps coordinates and lengths
# argument 4 = Input  gaps FASTA files prefixes
# argument 5 = Output scaffold FASTA file

import sys
import os.path
from operator import itemgetter

# command_line="python "
# for a in sys.argv:
#     command_line=command_line+' '+a
# print "------> COMMAND: "+command_line

MISSING_GAP_SEQ="NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"

# Auxiliary functions
def halve_marker(m):
    if m%2==0:
        return m/2 
    else:
        return (m+1)/2

base_complement = {"a":"t","A":"T","t":"a","T":"A","c":"g","C":"G","g":"c","G":"C","R":"","M":"","N":"N","n":""}
def complement(sequence):
    result = []
    for n in range(len(sequence)):
        result.append(base_complement[sequence[-(n+1)]])
    return ''.join(result)

# reading markers sequences
families_sequences_file=open(sys.argv[1],"r").readlines()
sequences={}
for l in families_sequences_file:
    if l[0]==">":
        fam=l.split('>')[1].rstrip('\n')
    else:
        sequences[int(fam)]=l.rstrip('\n')

# reading the scaffold
scaffold_file=open(sys.argv[2],"r").readlines()
order={}     # markers recorded in sequential order of reading in the scaffold
order_inv={} # inverse map from markers to position in order (valid only for non-repeat)
pos=1        # initial position
gaps_scaffold={} # list of gaps: ((adj1,adj2),gap_sequence) where adj1 and adj2 are undoubled markers (integers) and gap_sequence is a DNA sequence
for l in scaffold_file:
    if l[0]==">":
        genome_name=l.rstrip( )
    elif l[0]=="#":
        car_name=l[1:len(l)]
        car_type="L"  # by default, linear CAR
        car_start=pos # position of the first marker of the current scaffold
    elif l[0]=="_":
        car=l.rstrip('\n').rstrip(' ').split(' ') # list of markers in the current CAR
        for c in car:
            if c=="_C":
                car_type="C"  # circular CAR
            elif c!="_Q" and c!="Q_" and c!="C_":
                if c[0]=="-":
                    sign=1        # sign of the marker in the scaffold: 1=reversed, 0=unreversed
                    m=c[1:len(c)] # number of the marker
                else:
                    sign=0
                    m=c
                order[pos]=(m,car_name,sign,car_type,car_start)
                gaps_scaffold[pos]=((int(m),-1),MISSING_GAP_SEQ) # dummy gap, to be filled when reading the gaps
                order_inv[int(m)]=pos # reminer: not valid if m is a repeat but for its last occurrence 
                pos+=1

copy_number={}   # number of occurrences of each marker
for l in scaffold_file:
    if l[0]!=">" and l[0]!="#":
        car=l.rstrip('\n').rstrip(' ').split(' ') # list of markers in the current CAR
        for c in car:
            if c!="_C" and c!="_Q" and c!="Q_" and c!="C_":
                if c[0]=="-":
                    m=c[1:len(c)] # number of the marker
                else:
                    m=c
                copy_number[int(m)]=0
for l in scaffold_file:
    if l[0]!=">" and l[0]!="#":
        car=l.rstrip('\n').rstrip(' ').split(' ') # list of markers in the current CAR
        for c in car:
            if c!="_C" and c!="_Q" and c!="Q_" and c!="C_":
                if c[0]=="-":
                    m=c[1:len(c)] # number of the marker
                else:
                    m=c
                copy_number[int(m)]+=1

# reading the gaps to update gaps_scaffold
# entry ((adj1,adj2),gap_sequence) means there is a gap between markers adj1 and adj2, 
# adj1 appearing before adj2 in the scaffold order and gap_sequence being associated
# to the smaller marker
# entries are ordered in order of appearance in the scaffold order
gaps_file=open(sys.argv[3],"r").readlines()
for l in gaps_file:
    if l[0]==">":
        gap_line=l.rstrip( ).split(' ')
        gap=gap_line[0].split('>')[1]
        adj1=halve_marker(int(gap_line[2].split('-')[0]))
        adj2=halve_marker(int(gap_line[2].split('-')[1]))# assumption: adj1<adj2
        interval=l.split(' ')[len(l.split(' '))-1].rstrip('\n')
        if interval=="-1--1":
            lg1=0
            lg2=0
        else:
            lg1=int(gap_line[len(gap_line)-1].split('-')[0]) # lower bound of the length interval of the gap
            lg2=int(gap_line[len(gap_line)-1].split('-')[1]) # upper bound of the length interval of the gap

        if lg1!=0 and lg2!=0 and not os.path.isfile(sys.argv[4]+gap+".fasta_ancestral"):
            lg1=-1
            lg2=-1

        if lg1>0 or lg2>0: # the gap is not empty, so we need to find its sequence
            gap_sequence_file=open(sys.argv[4]+gap+".fasta_ancestral","r").readlines()
            for ll in gap_sequence_file:
                if ll[0]!=">":
                    gap_sequence=ll #.lower()
        elif lg1==0 and lg2==0: # the gap is empty
            gap_sequence=""
        elif lg1==-1 and lg2==-1: # error
            gap_sequence=MISSING_GAP_SEQ
            print "MISSING GAP SEQUENCE: "+l

        if len(gap_line)>5: # gap from a repeat spanning interval including two repeats
            ext1=halve_marker(int(gap_line[7])) # first non-repeat framing the interval
            ext2=halve_marker(int(gap_line[len(gap_line)-3])) # second non-repeat framing the interval
            pos=int(gap_line[4])/2 # position in the interval
            oext1=order_inv[ext1]  # position in the scaffold sequential order
            oext2=order_inv[ext2]  # same as above
            # assumption: a scaffold does not start by a repeat
            if oext1<oext2:        # the interval was read in the scaffold as written in the interval file
                a1=int(order[oext1+pos][0])   # leftmost marker of the adjacency
                a2=int(order[oext1+pos+1][0]) # rightmost marker of the adjacency
                gaps_scaffold[oext1+pos]=((a1,a2),gap_sequence)
            else: # the interval was read reversed
                a1=int(order[oext1-pos-1][0]) # same as above
                a2=int(order[oext1-pos][0])   # same as above
                gaps_scaffold[oext1-pos-1]=((a1,a2),gap_sequence)

        elif copy_number[adj1]==1 and copy_number[adj2]==1: # unique markers in the gap
            oadj1=order_inv[adj1] # position of smaller marker of adjacency in scaffold
            oadj2=order_inv[adj2] # position of larger marker of adjacency in scaffold
            if oadj1==oadj2-1:
                gaps_scaffold[oadj1]=((adj1,adj2),gap_sequence)
            elif oadj2==oadj1-1:
                gaps_scaffold[oadj2]=((adj2,adj1),gap_sequence)
            elif oadj1<oadj2: #circular chromosome finishing at oadj2
                gaps_scaffold[oadj2]=((adj2,adj1),gap_sequence)
            elif  oadj1>oadj2: #circular chromosome finishing at oadj1
                gaps_scaffold[oadj1]=((adj1,adj2),gap_sequence)

        else: # Mix of one repeat and one non-repeat
            if copy_number[adj1]==1:
                non_repeat=adj1                
                repeat=adj2
            else:
                non_repeat=adj2
                repeat=adj1
            onon_repeat=order_inv[non_repeat]
            if int(order[onon_repeat-1][0])==repeat and int(order[onon_repeat+1][0])!=repeat:
                gaps_scaffold[onon_repeat-1]=((repeat,non_repeat),gap_sequence)
            elif int(order[onon_repeat+1][0])==repeat and int(order[onon_repeat-1][0])!=repeat:
                gaps_scaffold[onon_repeat]=((non_repeat,repeat),gap_sequence)
            elif int(order[onon_repeat-1][0])==repeat and int(order[onon_repeat+1][0])==repeat:
                print "NON-REPEAT FRAMED BY REPEATS OF THE SAME MARKER "+str(repeat)+" "+str(non_repeat)

sequence=""    # ancestral genome sequence
current_car="" # current CAR whose sequence is being constructed
for g in gaps_scaffold:  # reading gaps by order of appearance in the scaffold
    gap=gaps_scaffold[g] # current gap
    m1=gap[0][0]         # current marker to write
    m1_info=order[g]     # info of current first marker of the gap
    car=m1_info[1]       # CAR of the current marker
    m2=gap[0][1]         # next marker, -1 if not in a gap seen in the gaps file

    if car!=current_car: # starting a new CAR
        if sequence!="": # the new CAR is not the first CAR
            sequence=sequence+"\n"
        sequence=sequence+genome_name+"."+car 
    current_car=car
    start_current_car=g # position of the first marker of the current CAR

    # adding sequence of the current marker
    sign=m1_info[2]  # sign, in the scaffold order, of the current marker
    if sign==0: # current marker not reversed
        sequence=sequence+sequences[m1]
    else:
        sequence=sequence+complement(sequences[m1])

    # adding sequence of the following gap if not at the end of a linear CAR
    if m2!=-1:
        if m1<m2: # the current marker is smaller than the next marker: the gap is associated to m1
            if sign==0: # smaller marker not reversed
                sequence=sequence+gap[1]
            else: # smaller marker reversed
                sequence=sequence+complement(gap[1])
        else: # the current marker is not smaller than the next marker, that is linked to the gap sequence
            # we need to find the orientation of m2
            # ennoying case: m2 is the beginning of the current circular scaffold
            if g<len(order) and int(order[g+1][0])==m2: # linear CAR or middle of a CAR
                sign=order[g+1][2]
            else: # m2 is not the following marker: end of circular CAR
                sign=order[start_current_car][2]
            if sign==0: # smaller marker not reversed
                sequence=sequence+gap[1]
            else: # smaller marker reversed
                sequence=sequence+complement(gap[1])
    else:
        sequence=sequence+MISSING_GAP_SEQ

output=open(sys.argv[5],"w")
output.write(sequence+'\n')
