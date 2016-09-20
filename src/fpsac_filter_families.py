# FPSAC: Fast Phylogenetic Scaffolding of Ancient Contigs.
# May 2013.
# Cedric Chauve: cedric.chauve@sfu.ca

# Compute some statistics on multiplicities and filter high
# multiplicity families
 
# argument 1 = Input  homologous families occurrences file
# argument 2 = Input  ancestral content file
# argument 3 = Input  multiplicity threshold for filtering families
# argument 4 = Output filtered families occurrences file

import string,sys

# command_line="python "
# for a in sys.argv:
#     command_line=command_line+' '+a
# print "------> COMMAND: "+command_line

families_file = open(sys.argv[1],"r").readlines()
content_file = open(sys.argv[2],"r").readlines()
THRESHOLD=int(sys.argv[3])
output = open(sys.argv[4],"w")

family = {}
indice = 0
while indice < len(families_file):
    if families_file[indice][0] == ">":
        nom_famille = families_file[indice].split()[0]
        family[nom_famille] = []
        indice = indice + 1
        while indice < len(families_file) and len(families_file[indice].split(":")) >= 2:
            family[nom_famille].append(families_file[indice])
            indice = indice + 1
    else:
        indice = indice + 1
    
content = {}
for ligne in content_file:
     mots = ligne.split()
     content[mots[0]] = int(mots[1])

print "STAT: ",len(content.keys()),len(family.keys()),"families"
families = content.keys()

rep = {}
for f in families:
    if content[f] <= THRESHOLD and content[f]>0:
        output.write(f+"\n")
        for h in family[f]:
            output.write(string.join(h.split()[:2]," ")+"\n")
	    contig = h.split()[2].split(",")
	    for c in contig:
		    rep[c.split(":")[0].split("(")[0]] = 1
        output.write("\n")
#	if content[f] > 1:
#		print "STAT: family with multiplicity ",f,"\t",content[f]
    else:
	   print "STAT: family with multiplicity ",f,"\t",content[f],"filtered"
      

print "STAT: ",len(rep.keys()),"represented contigs"
