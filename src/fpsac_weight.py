# FPSAC: Fast Phylogenetic Scaffolding of Ancient Contigs.
# May 2013.
# Cedric Chauve: cedric.chauve@sfu.ca

# Phylogenetic weighting
 
# argument 1 = Input  unweighted information file
# argument 2 = Input  species tree with marked ancestor
# argument 3 = Output weighted information file
# argument 4 = Input  doubled markers flag (d)

import sys,bm,weight

# command_line="python "
# for a in sys.argv:
#     command_line=command_line+' '+a
# print "------> COMMAND: "+command_line

double = len(sys.argv) > 4 and sys.argv[4] == 'd'
# True if the matrix has
# doubled markers to marker
# single marker weights
# correctly
														
if len(sys.argv) > 4 and sys.argv[4] != 'd':
	offset = int(sys.argv[4])
offset = 0
#endif
	
m=bm.BinaryMatrix()		# matrix
m.from_file(sys.argv[1])
weight.weight_matrix(m, sys.argv[2], double, offset)

f=file(sys.argv[3], 'w')		# weighted matrix file
f.write(str(m))
f.close()
