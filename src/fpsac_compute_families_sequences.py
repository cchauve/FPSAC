# FPSAC: Fast Phylogenetic Scaffolding of Ancient Contigs.
# May 2013.
# Cedric Chauve: cedric.chauve@sfu.ca

# Computing the DNA sequence of ancestral markers in families

# argument 1 = Input  families with contig names file
# argument 2 = Input  ancestral contigs FASTA file (sequence on one line)
# argument 3 = Output ancestral markers FASTA file
# argument 4 = Output ancestral markers length file

import sys
import graph_tools

# command_line="python "
# for a in sys.argv:
#     command_line=command_line+' '+a
# print "------> COMMAND: "+command_line

contigsfile=open(sys.argv[2],"r").readlines()
ctgs_seqs={}
for l in contigsfile:
    if l[0]==">":
        ctg_name=l.rstrip('\n').split('>')[1].split(' ')[0]
    else:
         ctgs_seqs[ctg_name]=l.rstrip('\n')

familiesfile=open(sys.argv[1],"r").readlines()
output_seq=open(sys.argv[3],"w")
output_lgt=open(sys.argv[4],"w")
for l in familiesfile:
    if l[0]==">":
        fam=l.rstrip('\n')
        seq_done=False
    else:
        if len(l)>1 and seq_done==False:
            ctg_name=l.split(' ')[2].rstrip('\n').split(',')[0].split(':')[0]
            coords=l.split(' ')[2].rstrip('\n').split(',')[0].split(':')[1]
            start=int(coords.split('-')[0])-1
            end=int(coords.split('-')[1])
            output_lgt.write(fam+' '+str(end-start)+'\n')
            if len(ctg_name.rstrip('(rev)'))==len(ctg_name):
                output_seq.write(fam+'\n'+ctgs_seqs[ctg_name][start:end]+'\n')
            else:
                ctg_name=ctg_name.rstrip('(rev)')
                output_seq.write(fam+'\n'+ctgs_seqs[ctg_name][start:end][::-1]+'\n')
            seq_done=True
