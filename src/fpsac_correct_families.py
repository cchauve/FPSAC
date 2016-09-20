# FPSAC: Fast Phylogenetic Scaffolding of Ancient Contigs.
# May 2013.
# Cedric Chauve: cedric.chauve@sfu.ca

# Correcting homologous families to discard occurrences with no orientation

# argument 1 = Input  homologous families occurrences file
# argument 2 = Input  homologous families profile file
# argument 3 = Output corrected families occurrences
# argument 4 = Output corrected families profiles
# argument 5 = Output list of removed occurrences

import sys

# command_line="python "
# for a in sys.argv:
#     command_line=command_line+' '+a
# print "------> COMMAND: "+command_line

families_file=open(sys.argv[1],"r").readlines()
profiles_file=open(sys.argv[2],"r").readlines()
output_families=open(sys.argv[3],"w")
output_profiles=open(sys.argv[4],"w")
output_errors=open(sys.argv[5],"w")

# Reading profiles
profiles={}
for l in profiles_file:
    if len(l)>1:
        if l[0]==">":
            fam=l.split('>')[1].rstrip('\n')
            profiles[fam]={}
        else:
            species=l.split(' ')[0]
            copynumber=int(l.split(' ')[1].rstrip('\n'))
            profiles[fam][species]=copynumber

# Reading families and writing output
fam=">"
for l in families_file:
    if len(l)>1:
        if l[0]==">":
            fam=l.split('>')[1].rstrip('\n')
            output_families.write(l)
            output_profiles.write(l)
        elif l.split(' ')[1]=='+' or  l.split(' ')[1]=='-':
            output_families.write(l)
        else:
            species=l.split('.')[0]
            profiles[fam][species]-=1
            output_errors.write(fam+" "+l)
    else:
        for s in profiles[fam]:
            output_profiles.write(s+" "+str(profiles[fam][s])+"\n")
        output_profiles.write('\n')
        output_families.write('\n')
for s in profiles[fam]:
    output_profiles.write(s+" "+str(profiles[fam][s])+"\n")
