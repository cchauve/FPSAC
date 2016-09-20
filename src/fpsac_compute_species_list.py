# FPSAC: Fast Phylogenetic Scaffolding of Ancient Contigs.
# May 2013.
# Cedric Chauve: cedric.chauve@sfu.ca

# List all species in aphylogenetic tree
 
# argument 1 = Input  phylogenetic tree
# argument 2 = Input  species selected (ALL, INGROUPS, OUTGROUPS)
# argument 3 = Output species file

import sys,string,time,math

# command_line="python "
# for a in sys.argv:
#     command_line=command_line+' '+a
# print "------> COMMAND: "+command_line

tree_file=open(sys.argv[1],"r").readline()
type=sys.argv[2]
output_file=open(sys.argv[3],"w")

##########################################
##### READING SPECIES TREE FILE ##########
##########################################

ancestor_id = -1
duplication_id = -1

id=0
tree={}
pile=[]
t=0
while t<len(tree_file):
    if tree_file[t]=="(":
        id=id+1
        tree[id]=["N"+str(id),-1,-1,-1]
        if len(pile)>0:
            tree[id][1]=pile[-1]
            if(tree[pile[-1]][2]==-1):
                tree[pile[-1]][2] = id
            else:
                tree[pile[-1]][3] = id
	
        pile.append(id)
        t=t+1
    elif tree_file[t]==")":
        t=t+1
        if tree_file[t]=="@":
            ancestor_id=pile[-1]
            t=t+1
        if tree_file[t]==":":
            t=t+1
            begin=t
            while tree_file[t]!="," and tree_file[t]!=")":
                t=t+1
            length=float(tree_file[begin:t])
        if tree_file[t]==";" or tree_file[t]=="\n" or tree_file == " ":
            t=t+2            
        del pile[-1]
    elif tree_file[t]==",":
        t=t+1
    else:
        id=id+1
        tree[id]=["",-1,-1,-1]
        if len(pile)>0:
            tree[id][1]=pile[-1]
            if(tree[pile[-1]][2]==-1):
                tree[pile[-1]][2] = id
            else:
                tree[pile[-1]][3] = id
        pile.append(id)
        begin=t
        while  t < len(tree_file) and tree_file[t]!="," and tree_file[t]!=")" and tree_file[t]!=":":
            t=t+1
        nom=tree_file[begin:t]
        tree[pile[-1]][0]=nom
        if t < len(tree_file) and tree_file[t]==":":
            t=t+1
            begin=t
            while tree_file[t]!="," and tree_file[t]!=")":
                t=t+1
            length=float(tree_file[begin:t])
        del pile[-1]


def partition_genomes(genome1,genome2,genome3,id,part_id):
	if(tree[id][2]==-1 and tree[id][3]==-1):
		if(part_id==1):
			genome1.append(tree[id][0])
		if(part_id==2):
			genome2.append(tree[id][0])
		if(part_id==3):
			genome3.append(tree[id][0])

	elif(id==ancestor_id):
		partition_genomes (genome1,genome2,genome3,tree[id][2],2)
		partition_genomes (genome1,genome2,genome3,tree[id][3],3)
	elif(id!=duplication_id):
		partition_genomes (genome1,genome2,genome3,tree[id][2],part_id)
		partition_genomes (genome1,genome2,genome3,tree[id][3],part_id)

genomes1=[]
genomes2=[]
genomes3=[]

partition_genomes (genomes1,genomes2,genomes3,1,1)

nb_species = len(genomes1) + len(genomes2) + len(genomes3)

genome1_id=[]
genome2_id=[]
genome3_id=[]

genome_name=[]

i = 0
m={}
for s in genomes1:
	m[s]=i
	genome1_id.append(i)
	genome_name.append(s)
	i+=1
for s in genomes2:
	m[s]=i
	genome2_id.append(i)
	genome_name.append(s)
	i+=1
for s in genomes3:
	m[s]=i
	genome3_id.append(i)
	genome_name.append(s)
	i+=1

if type=="ALL":
	output_file.write("#outgroup\n")
#endif

if type=="ALL" or type=="OUTGROUPS":
	for spe in genomes1:
		output_file.write(spe+"\n")
#endif

if type=="ALL":
	output_file.write("#ingroup\n")
#endif

if type=="ALL" or type=="INGROUPS":
	for spe in genomes2:
		output_file.write(spe+"\n")
	for spe in genomes3:
		output_file.write(spe+"\n")
#endif

output_file.close()
