# FPSAC: Fast Phylogenetic Scaffolding of Ancient Contigs.
# May 2013.
# Cedric Chauve: cedric.chauve@sfu.ca

# Reformatting a BLAST hits file ij order to have the following format for each hit:
# target_species_name target_chromosome target_start target_end orientation(1/-1) query_name query_start query_end query_length

# argument 1: Input  BLAST hits file
#             obtained with option -outfmt 6 (output format 6)
# argument 2: Input  file recording the length of the query sequences
# argument 3: Output file of reformatted hits

import sys,string

# command_line="python "
# for a in sys.argv:
#     command_line=command_line+' '+a
# print "------> COMMAND: "+command_line

blast_file=open(sys.argv[1],"r").readlines()
length_file=open(sys.argv[2],"r").readlines()
output_file=open(sys.argv[3],"w")

###########################################################################
###############  READING THE LENGTH FILE###################################
###########################################################################
ctgs_lengths={}
for l in length_file:
    line=l.split(' ')
    ctgs_lengths[line[0].split('>')[1]]=line[2].rstrip('\n')


###########################################################################
###############  READING THE HITS FILE ####################################
###########################################################################

hits = {}
for l in range(len(blast_file)):
    line = blast_file[l].split()
    name_contig = line[0]
    length_contig = ctgs_lengths[name_contig]
    species = line[1].split('.')[0]
    chromosome = line[1].split('.')[1]
    start_contig = int(line[6])
    stop_contig = int(line[7])
    start_genome = int(line[8])
    stop_genome = int(line[9])
    evalue = float(line[-2])
    if stop_genome >  start_genome:
        direction = 1
    else:
        direction = -1
        start_genome = int(line[9])
        stop_genome = int(line[8])
            
    hits[l+1] = [species,chromosome, # species_name, chromosome_name
                 start_genome,stop_genome,direction,  # start, stop, dir on the genome
                 name_contig,start_contig,stop_contig,length_contig,False] # name, start, stop, dir on the contig


###########################################################################
###############  WRITING FORMATTED HITS ###################################
###########################################################################

ident_hits = hits.keys()
for h in ident_hits:
    for g in hits[h][:-1]:
        output_file.write(str(g)+" ")
    output_file.write("\n")

