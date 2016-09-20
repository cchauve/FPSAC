# FPSAC: Fast Phylogenetic Scaffolding of Ancient Contigs.
# May 2013.
# Cedric Chauve: cedric.chauve@sfu.ca

# Mapping extant annotations to ancestral genome

import sys
from operator import itemgetter

# arg 1: Input  ancestral map
# arg 2: Input  extant annotations
# arg 3: Output annotated ancestral map

annotations_file=open(sys.argv[2],"r").readlines()
genomes=[]
annotations={}
for a in annotations_file:
    print a
    genome=a.split(' ')[1].split('.')[0]
    if  genome not in genomes:
        genomes.append(genomes)
        annotations[genome]=[]
    genome=a.split(' ')[1].split('.')[1]
    start=a.split(' ')[2]
    end=a.split(' ')[3]
    annotation=a.rstrip( ).split(' ')[5:-1]
    annotations[genome].append((chr,start,end,annotation))

for g in genomes:
    annotations[g].sort(key=itemgetter(0,1))

