# FPSAC: Fast Phylogenetic Scaffolding of Ancient Contigs.
# May 2013.
# Cedric Chauve: cedric.chauve@sfu.ca

# Takes a phylogenetic tree and a set of profiles and produces the
# ancestral numbers by minimizing the number of duplications and
# losses along the phylogenetic tree
 
# argument 1 = Input  species tree with marked ancestor location
# argument 2 = Input  homologous families profile file
# argument 3 = Output ancestral content file

import sys,string,math,time
from tree import *

# command_line="python "
# for a in sys.argv:
#     command_line=command_line+' '+a
# print "------> COMMAND: "+command_line

time0 = time.time()

BAR = True
BAR = False

species_tree = readTree(open(sys.argv[1],"r").readline())
families_file = open(sys.argv[2],"r").readlines()
output_file = open(sys.argv[3],"w")


leaves = getLeaves(species_tree,getRoot(species_tree))
extant_species = []
for l in leaves:
    spe = getName(species_tree,l).lower()
    extant_species.append(spe)
    
ancestral_species = []
for node in getNodes(species_tree):
    if not isLeaf(species_tree,node):
        leaves = getLeaves(species_tree,node)
        species = []
        for l in leaves:
            species.append(getName(species_tree,l).lower())
        species.sort()
        ancestral_species.append([getName(species_tree,node),species])

postorder = getNodes(species_tree)
index = {}
for p in range(len(postorder)):
    index[postorder[p]] = p
isapostorder = False
while not isapostorder:
    isapostorder = True
    for p in range(len(postorder)):
        if not isRoot(species_tree,postorder[p]):
            i = index[getParent(species_tree,postorder[p])]
            if p > i:
                isapostorder = False
                index[postorder[p]] = i
                index[postorder[i]] = p
                temp = postorder[p]
                postorder[p] = postorder[i]
                postorder[i] = temp

genome = {}
for s in extant_species:
    genome[s] = 0
for s in ancestral_species:
    genome[s[0]] = 0
numero_famille = 0
indice = 0
nbdup = 0
while indice < len(families_file):
    if families_file[indice][0] == ">":
        
        nom_famille = families_file[indice][1:].split()[0]
        profils = {}
        for s in extant_species:
            profils[s] = 0

        indice = indice + 1
        while indice < len(families_file) and len(families_file[indice].split()) >= 2:
            spe = families_file[indice].split()[0].lower()
            nombre = int(families_file[indice].split()[1])
            if spe in extant_species:
                profils[spe] = profils[spe] + nombre
                genome[spe] = genome[spe] + nombre
            indice = indice + 1
        if indice < len(families_file) and len(families_file[indice])>1:
            indice = indice-1
            
        value = {}
        for current in postorder:
            if isLeaf(species_tree,current):
                spe = getName(species_tree,current).lower()
                value[current] = [profils[spe],profils[spe]]
            else:
                mindesmax = min(value[getChildren(species_tree,current)[0]][1],
                                value[getChildren(species_tree,current)[1]][1])
                maxdesmin = max(value[getChildren(species_tree,current)[0]][0],
                                value[getChildren(species_tree,current)[1]][0])
                debut = min(mindesmax,maxdesmin)
                fin = max(mindesmax,maxdesmin)
                if (fin == 0 and (value[getChildren(species_tree,current)[0]][1] > 0 or
                                  value[getChildren(species_tree,current)[1]][1] > 0)):
                    fin = 1
                value[current] = [debut,fin]

        postorder.reverse()
        for current in postorder:
            if isLeaf(species_tree,current):
                spe = getName(species_tree,current).lower()
                value[current] = value[current][0]
                val = value[getParent(species_tree,current)]
                if val < value[current]:
                    nbdup = nbdup + value[current] - max(val,1)
            else:
                spe = getName(species_tree,current)
                if isRoot(species_tree,current):
                    if value[current][0] == 0:
                        valc1 = value[getChildren(species_tree,current)[0]][1]
                        valc2 = value[getChildren(species_tree,current)[1]][1]
                        if valc1 == 0 or valc2 == 0:
                            value[current] = 0
                        else:
                            value[current] = 1
                    else:
                        value[current] = value[current][0]
                else:
                    val = value[getParent(species_tree,current)]
                    if val == 0 and value[current][1] > 0:
                        valc1 = value[getChildren(species_tree,current)[0]][1]
                        valc2 = value[getChildren(species_tree,current)[1]][1]
                        if valc1 == 0 or valc2 == 0:
                            value[current] = 0
                        else:
                            value[current] = max(1,value[current][0])
                            nbdup = nbdup + value[current]
                    else:
                        if val >= value[current][0] and val <= value[current][1]:
                            value[current] = val
                        elif val < value[current][0]:
                            value[current] = value[current][0]
                            nbdup = nbdup + value[current] - max(val,1)
                        else:
                            value[current] = value[current][1]
                profils[spe] = value[current]
                genome[spe] = genome[spe] + value[current]

        if getAncestor(species_tree) != -1:
            output_file.write(">"+nom_famille+" "+str(profils[getName(species_tree,getAncestor(species_tree))])+"\n")
        else:
            output_file.write(">"+nom_famille+"\n")
            for n in postorder:
                if not isLeaf(species_tree,n):
                    output_file.write(getName(species_tree,n)+" "+str(profils[getName(species_tree,n)])+"\n")
            output_file.write("\n")
        postorder.reverse()
        

    indice = indice + 1

#print
#print nbdup , "duplications"


