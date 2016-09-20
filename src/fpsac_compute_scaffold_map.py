# FPSAC: Fast Phylogenetic Scaffolding of Ancient Contigs.
# May 2013.
# Cedric Chauve: cedric.chauve@sfu.ca

# Computing a scaffold map

# argument 1 = Input  families file, with contig names
# argument 2 = Input  families FASTA file
# argument 3 = Input  scaffold (undoubled)
# argument 4 = Input  gaps coordinates and lengths
# argument 5 = Input  gaps FASTA files prefixes
# argument 6 = Output map file

import sys
import os.path
from operator import itemgetter

# command_line="python "
# for a in sys.argv:
#     command_line=command_line+' '+a
# print "------> COMMAND: "+command_line

MISSING_GAP_SEQ_LG=len("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN")

# Auxiliary functions
def halve_marker(m):
    if m%2==0:
        return m/2 
    else:
        return (m+1)/2

# reading markers occurrences
families_file=open(sys.argv[1],"r").readlines()
families_extant={}
families_contigs={}
for l in families_file:
    if l[0]==">":
        fam=int(l.split('>')[1].rstrip('\n'))
        families_extant[fam]=[]
        families_contigs[fam]=[]
    elif len(l)>1:
        l1=l.rstrip('\n').split(' ')        
        families_extant[fam].append(l1[0]+" "+l1[1])
        contig_occurrences=l1[2].split(',')
        for contig_occurrence in contig_occurrences:
            if not contig_occurrence in families_contigs[fam]:
                families_contigs[fam].append(contig_occurrence)

# reading markers sequences
families_sequences_file=open(sys.argv[2],"r").readlines()
families_sequences_length={} # lengths of mqrkers sequences
for l in families_sequences_file:
    if l[0]==">":
        fam=int(l.split('>')[1].rstrip('\n'))
    else:
        families_sequences_length[fam]=len(l)-1

# reading the scaffold
scaffold_file=open(sys.argv[3],"r").readlines()
order={}     # markers recorded in sequential in order of reading in the scaffold
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
                gaps_scaffold[pos]=((int(m),-1),MISSING_GAP_SEQ_LG) # dummy gap, to be filled when reading the gaps
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
# entry ((adj1,adj2),gap_sequence,occurrences) means there is a gap between markers adj1 and adj2, 
# adj1 appearing before adj2 in the scaffold order and gap_sequence being associated
# to the smaller marker; occurrences record the gaps extant occurrences
# entries are ordered in order of appearance in the scaffold order
gaps_file=open(sys.argv[4],"r").readlines()
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

        if lg1!=0 and lg2!=0 and not os.path.isfile(sys.argv[5]+gap+".fasta_ancestral"):
            lg1=-1
            lg2=-1

        if lg1>0 or lg2>0: # the gap is not empty, so we need to find  sequence
            gap_sequence_file=open(sys.argv[5]+gap+".fasta_ancestral","r").readlines()
            for ll in gap_sequence_file:
                if ll[0]!=">":
                    gap_sequence=ll #.lower()
        elif lg1==0 and lg2==0: # the gap is empty
            gap_sequence=""
        elif lg1==-1 and lg2==-1: # error
            gap_sequence=MISSING_GAP_SEQ
            print "MISSING GAP SEQUENCE: "+l

        if len(gap_line)>5: # gap from a repeat spanning interval
            ext1=halve_marker(int(gap_line[7])) # first non-repeat framing the interval
            ext2=halve_marker(int(gap_line[len(gap_line)-3])) # second non-repeat framing the interval
            pos=int(gap_line[4])/2 # position in the interval
            oext1=order_inv[ext1]  # position in the scaffold sequential order
            oext2=order_inv[ext2]  # same as above
            # assumption: a scaffold does not start by a repeat
            if oext1<oext2:        # the interval was read in the scaffold as written in the interval file
                a1=int(order[oext1+pos][0])   # leftmost marker of the adjacency
                a2=int(order[oext1+pos+1][0]) # rightmost marker of the adjacency
                gaps_scaffold[oext1+pos]=((a1,a2),len(gap_sequence),[])
                gap_index=oext1+pos
            else: # the interval was read reversed
                a1=int(order[oext1-pos-1][0]) # same as above
                a2=int(order[oext1-pos][0])   # same as above
                gaps_scaffold[oext1-pos-1]=((a1,a2),len(gap_sequence),[])
                gap_index=oext1-pos-1

        elif copy_number[adj1]==1 and copy_number[adj2]==1: # unique markers in the gap
            oadj1=order_inv[adj1] # position of smaller marker of adjacency in scaffold
            oadj2=order_inv[adj2] # position of larger marker of adjacency in scaffold
            if oadj1==oadj2-1:
                gaps_scaffold[oadj1]=((adj1,adj2),len(gap_sequence),[])
                gap_index=oadj1
            elif oadj2==oadj1-1:
                gaps_scaffold[oadj2]=((adj2,adj1),len(gap_sequence),[])
                gap_index=oadj2
            elif oadj1<oadj2: #circular chromosome finishing at oadj2
                gaps_scaffold[oadj2]=((adj2,adj1),len(gap_sequence),[])
                gap_index=oadj2
            elif  oadj1>oadj2: #circular chromosome finishing at oadj1
                gaps_scaffold[oadj1]=((adj1,adj2),len(gap_sequence),[])
                gap_index=oadj1

        else: # Mix of one repeat and one non-repeat
            if copy_number[adj1]==1:
                non_repeat=adj1                
                repeat=adj2
            else:
                non_repeat=adj2
                repeat=adj1
            onon_repeat=order_inv[non_repeat]
            if int(order[onon_repeat-1][0])==repeat and int(order[onon_repeat+1][0])!=repeat:
                gaps_scaffold[onon_repeat-1]=((repeat,non_repeat),len(gap_sequence),[])
                gap_index=onon_repeat-1
            elif int(order[onon_repeat+1][0])==repeat and int(order[onon_repeat-1][0])!=repeat:
                gaps_scaffold[onon_repeat]=((non_repeat,repeat),len(gap_sequence),[])
                gap_index=onon_repeat
            elif int(order[onon_repeat-1][0])==repeat and int(order[onon_repeat+1][0])==repeat:
                print "NON-REPEAT FRAMED BY REPEATS OF THE SAME MARKER "+str(repeat)+" "+str(non_repeat)

        # else: # unique markers gap
        #     oadj1=order_inv[adj1] # position of smaller marker of adjacency in scaffold
        #     oadj2=order_inv[adj2] # position of larger marker of adjacency in scaffold
        #     if oadj1==oadj2-1:
        #         gaps_scaffold[oadj1]=((adj1,adj2),len(gap_sequence),[])
        #         gap_index=oadj1
        #     else: # includes cases of oadj2==oadj1-1 and oadj2 is the ending position of a circular scaffold started in oadj1
        #         gaps_scaffold[oadj2]=((adj2,adj1),len(gap_sequence),[])
        #         gap_index=oadj2

    elif len(l)>1:
        gaps_scaffold[gap_index][2].append(l.rstrip('\n'))


output_map=open(sys.argv[6],"w")
current_car=""
pos=1 # current physical position in ancestral sequence
current_car="" # current CAR whose sequence is being constructed
for g in gaps_scaffold:  # reading gaps by order of appearance in the scaffold
    gap=gaps_scaffold[g] # current gap
    m1=gap[0][0]         # current marker to write
    m1_info=order[g]     # info of current first marker of the gap
    car=m1_info[1]       # CAR of the current marker
    sign=m1_info[2]      # sign, in the scaffold order, of the current marker
    m2=gap[0][1]         # next marker, -1 if not in a gap seen in the gaps file

    if car!=current_car: # starting a new CAR
        output_map.write(genome_name+"."+car)
    current_car=car
    start_current_car=g # position of the first marker of the current CAR

    # adding the current marker
    pos_end=pos+families_sequences_length[int(m1)]-1
    if sign==0: # current marker not reversed
        output_map.write("MARKER "+str(m1)+" "+str(pos)+"-"+str(pos_end)+" + \n")
    else:
        output_map.write("MARKER "+str(m1)+" "+str(pos)+"-"+str(pos_end)+" - \n")
    for occurrence in families_contigs[int(m1)]:
        output_map.write("   "+occurrence+"\n")
    families_extant[int(m1)].sort()
    for occurrence in families_extant[int(m1)]:
        output_map.write("   "+occurrence+"\n")
    pos=pos_end+1

    # adding the following gap if not at the end of a linear CAR
    if m2!=-1:
        pos_end=pos+gap[1]-1
        if g<len(order) and int(order[g+1][0])==m2: # linear CAR or middle of a CAR
            sign=order[g+1][2]
        else: # m2 is not the following marker: end of circular CAR
            sign=order[start_current_car][2]
        if pos<=pos_end:
            if sign==0:
                output_map.write("GAP "+str(m1)+"-"+str(m2)+" "+str(gap[1])+" "+str(pos)+"-"+str(pos_end)+" + \n")
            else:
                output_map.write("GAP "+str(m1)+"-"+str(m2)+" "+str(gap[1])+" "+str(pos)+"-"+str(pos_end)+" - \n")
        else:
            output_map.write("GAP "+str(m1)+"-"+str(m2)+" "+str(gap[1])+"\n")        
        gap[2].sort()
        for occurrence in gap[2]:
            start=occurrence.split(':')[1].split('-')[0]
            end=occurrence.split(':')[1].split('-')[1]
            output_map.write("   "+occurrence+"\n")
        pos=pos_end+1
    else:
        pos_end=pos+gap[1]-1
        output_map.write("GAP MISSING_GAP "+str(gap[1])+" "+str(pos)+"-"+str(pos_end)+" MISSING_GAP \n")
        pos=pos_end+1
