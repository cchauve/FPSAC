# FPSAC: Fast Phylogenetic Scaffolding of Ancient Contigs.
# May 2013.
# Cedric Chauve: cedric.chauve@sfu.ca

# Integrate a set of chosen adjacencies between CARs extremities to
# the existing gaps and scaffold order

# argument 1: Input  Chosen CARs extremities
# argument 2: Input  Scaffold order doubled
# argument 3: Input  Dollo gaps file
# argument 4: Input  gaps FASTA files prefixes
# argument 5: Output Combined scaffold order file
# argument 6: Output Combined gaps file

import sys

# command_line="python "
# for a in sys.argv:
#     command_line=command_line+' '+a
# print "------> COMMAND: "+command_line

GAP_CAR="NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"

updated_gaps_file=open(sys.argv[6],"w")
Dollo_gaps_file=open(sys.argv[3],"r").readlines()
gap_id=1
for g in Dollo_gaps_file:
    updated_gaps_file.write(g)
    if g[0]==">":
        gap_id+=1
CARS_adj_file=open(sys.argv[1],"r").readlines()
for a in CARS_adj_file:
    a1=int(a.split(':')[1].split(' ')[0])
    a2=int(a.split(':')[1].split(' ')[1])
    gap=a.rstrip('\n').split('#')[1:]
    updated_gaps_file.write(">gap_"+str(gap_id)+" adjacency "+str(a1)+"-"+str(a2)+" CAR_extremities_adjacency 50-50\n")
    gap_sequence_file=open(sys.argv[4]+"/gap_"+str(gap_id)+".fasta_ancestral","w")
    gap_sequence_file.write(GAP_CAR)
    for g in gap:
        orientation=g.split(' ')[1].rstrip( ) # orientation of the adjacenc
        # + means it was seen as a1 a2 in the genome
        # - means it was seen as a2 a1 in the genome
        sgn="+" # indicates on which strand to get the gap sequence
        # rule: gap is taken on on the same strand than a1, the smaller marker
        if (orientation=="+" and a1%2==1) or (orientation=="-" and a1%2==0):
            sgn="-"
        updated_gaps_file.write(g.split(' ')[0]+" "+sgn+"\n")
    updated_gaps_file.write("\n")
    gap_id+=1

updated_order_file=open(sys.argv[5],"w")
scaffold_order_file=open(sys.argv[2],"r").readlines()
CARS=[]
for s in scaffold_order_file:
    if s[0]==">":
        updated_order_file.write(s)
    elif s[0]=="#":
        car_name=s.rstrip('\n')
    else:
        s1=s.split(' ')
        ext1=s1[1]
        ext2=s1[-2]
        car=s1[1:-1]
        CARS.append((car_name,ext1,ext2,car,s1[0],s1[-1]))

for a in CARS_adj_file:
    a1=a.split(':')[1].split(' ')[0]
    a2=a.split(':')[1].split(' ')[1]
    modified_CARS=[]
    for c in CARS:
        found=False
        if c[1]==a1:
            c1=c
            rev1=True
            found=True
        if c[2]==a1:
            c1=c
            rev1=False
            found=True
        if c[1]==a2:
            c2=c
            rev2=False
            found=True
        if c[2]==a2:
            c2=c
            rev2=True
            found=True
        if found==False:
            modified_CARS.append(c)
    if c1!=c2:
        if rev1:
            c1[3].reverse()
        if rev2:
            c2[3].reverse()
        new_car=c1[3]+c2[3]
        modified_CARS.append((c1[0],new_car[0],new_car[-1],new_car,"_Q","Q_\n"))
    else:
        modified_CARS.append((c1[0],c1[1],c1[2],c[3],"_C","C_\n"))
    CARS=[]
    for c in modified_CARS:
        CARS.append(c)

car_id=1
for c in CARS:
    updated_order_file.write("#CAR"+str(car_id)+"\n")
    updated_order_file.write(c[4]+" "+' '.join(c[3])+" "+c[5])
    car_id+=1
