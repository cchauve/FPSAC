# FPSAC: Fast Phylogenetic Scaffolding of Ancient Contigs.
# May 2013.
# Cedric Chauve: cedric.chauve@sfu.ca

# Filtering adjacencies using a maximum weighted matching algorithm

# argument 1 = Input  weighted adjacencies file
# argument 2 = Input  ancestral content file
# argument 3 = Input  doubled marker flag (0/1)
# argument 4 = Output kept adjacencies file
# argument 4 = Output discarded adjacencies file

import sys,bm,math,decimal
import linearization as lz

# command_line="python "
# for a in sys.argv:
#     command_line=command_line+' '+a
# print "------> COMMAND: "+command_line

m = bm.BinaryMatrix()		# matrix
m.from_file(sys.argv[1])

# Reading a multiplicity file giving the multiplicity of each marker	
f=open(sys.argv[2])
lines=f.read().split('\n')
f.close()
multiplicity={}
i=0
while i<len(lines):
	if len(lines[i])>0:
		if sys.argv[3]=='0':
			splitpair=lines[i].split(' ')
			splitpair[0]=int(splitpair[0][1:])
			splitpair[1]=int(splitpair[1])
			multiplicity[splitpair[0]]=splitpair[1]
# If doubled, assign the same multiplicity for the head and the tail of a doubled marker
		elif sys.argv[3]=='1':
                        splitpair=lines[i].split(' ')
                        splitpair[0]=int(splitpair[0][1:])
                        splitpair[1]=int(splitpair[1])
                        multiplicity[2*splitpair[0]-1]=splitpair[1]
			multiplicity[2*splitpair[0]]=splitpair[1]
	i+=1
# If doubled, discard adjacencies between different ends of the same marker. Store them in another matrix so we can use them later.
trivial_adjacencies=bm.BinaryMatrix()
if sys.argv[3]=='1':
	i=0
	while i<m._height:
		list_of_markers=list(m.get_row(i))
		if len(list_of_markers)>2:
			print 'Adjacencies only!'
			sys.exit()
		if math.ceil(decimal.Decimal(list_of_markers[0])/2)==math.ceil(decimal.Decimal(list_of_markers[1])/2):
			trivial_adjacencies.add_row_info(m.get_row_info(i))
			m.remove_row(i)
			i-=1
		i+=1
# Call the linearization routine with suitable arguments
if sys.argv[3]=='0':
	[new_matrix,discarded_rows]=lz.linearization(m,multiplicity,2)
elif sys.argv[3]=='1':
        [new_matrix,discarded_rows]=lz.linearization(m,multiplicity,1)
# Reintroduce adjacencies between ends of the same marker
	for i in xrange(0,trivial_adjacencies._height):
		new_matrix.add_row_info(trivial_adjacencies.get_row_info(i))
# Store retained rows and discarded rows 
o=file(sys.argv[4],'w')
o.write(str(new_matrix))
o.close()
o=open(sys.argv[5],'w')
discarded_rows.write(o.write)
o.close()

