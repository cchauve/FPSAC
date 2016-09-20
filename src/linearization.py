import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/src/NETWORKX/lib/python2.7/site-packages/')

import math
import bm
import networkx as nx

def linearization_directed(m,multiplicity): # multiplicity is a dictionary, with markers and their reversals as keys.
	# Creating auxiliary adjacency graph

	# Step 1: Adding 2 vertices for each copy of a marker
	G=nx.Graph()
	for col in m.get_support():
		for k in xrange(0,multiplicity[col]):
			G.add_node(str(col)+'_'+str(k))
	
	
	# Step 2: Adding 2 vertices for each adjacency.
	# Also adding corresponding edges and giving them weights as given in the matrix
	
	for i in xrange(0,m._height):#row in m:
		row_weight=float(m.get_row_info(i)._weight)
		row_id=m.get_row_info(i)._id
		if len(m.get_row(i))>2:
			print 'Adjacencies only!'
			sys.exit(0)
		for col in m.get_row(i):
			G.add_node('edge_'+str(i)+'_'+str(col))
			for k in xrange(0,multiplicity[col]):
				G.add_edge(str(col)+'_'+str(k),'edge_'+str(i)+'_'+str(col),weight=row_weight)
		column_set=list(m.get_row(i))
		G.add_edge('edge_'+str(i)+'_'+str(column_set[0]),'edge_'+str(i)+'_'+str(column_set[1]),weight=row_weight)
	
	
	
	original_weight=m.get_weight()
	
	# Maximum weight matching algorithm called on the auxiliary graph.
	# This gives us A maximum matching, and, correspondingly, a maximum weight 2m-matching in the matrix
	# This DOES NOT maximize the number of adjacencies chosen.
	matching=nx.max_weight_matching(G)
	present=0

	new_matrix=bm.BinaryMatrix()
	discarded_rows=bm.BinaryMatrix()
	
	# Using the matching to find a component mCi1P respecting spanning subgraph of maximum weight.
	# The corresponding matrices- for the adjacencies retained, and those discarded- are the outputs.
	for vertex in matching.keys():
		if (not vertex[0:5]=='edge_' or not matching[vertex][0:5]=='edge_'):
			# Condition for retaining an adjacency: both ends of the corresponding adjacency must be matched, but not to each other.
			row=0
			first_end=None
			other_end=None
			if vertex[0:5]=='edge_':
				k=5
				while vertex[k]!='_':
					k+=1
				row=int(vertex[5:k])
				first_end=int(vertex[k+1:])
	                        for col in m.get_row_info(row)._set:
	                                if col!=first_end:
						other_end=col
			elif matching[vertex][0:5]=='edge_':
				k=5
				while matching[vertex][k]!='_':
					k+=1
				row=int(matching[vertex][5:k])

				first_end=int(matching[vertex][k+1:])

				for col in m.get_row_info(row)._set:
					if col!=first_end:
						other_end=col
			# Check if both ends are matched
			if  'edge_'+str(row)+'_'+str(other_end) in matching.keys():
				if m.get_row_info(row) not in new_matrix:
					new_matrix.add_row_info(m.get_row_info(row))
			else:
			# Case when one edge-vertex is matched to a copy of a column, but the other is unmatched.
				if m.get_row_info(row) not in discarded_rows:
					discarded_rows.add_row_info(m.get_row_info(row))
		elif (vertex[0:5]=='edge_' and  matching[vertex][0:5]=='edge_'):
			k=5
		# The two edge-vertices of an adjacency are matched to each other.
			while vertex[k]!='_':
				k+=1
			if m.get_row_info(int(vertex[5:k])) not in discarded_rows:
				discarded_rows.add_row_info(m.get_row_info(int(vertex[5:k])))
	 
	return [new_matrix,discarded_rows]


def linearization(m,multiplicity,f): # multiplicity is a dictionary, with markers and their reversals as keys.
	# Creating auxiliary adjacency graph

	# Step 1: Adding 2 vertices for each copy of a marker
	G=nx.Graph()
	for col in m.get_support():
		for k in xrange(0,f*multiplicity[col]):
			G.add_node(str(col)+'_'+str(k))
	
	
	# Step 2: Adding 2 vertices for each adjacency.
	# Also adding corresponding edges and giving them weights as given in the matrix
	
	for i in xrange(0,m._height):#row in m:
		row_weight=float(m.get_row_info(i)._weight)
		row_id=m.get_row_info(i)._id
		if len(m.get_row(i))>2:
			print 'Adjacencies only!'
			sys.exit(0)
		for col in m.get_row(i):
			G.add_node('edge_'+str(i)+'_'+str(col))
			for k in xrange(0,f*multiplicity[col]):
				G.add_edge(str(col)+'_'+str(k),'edge_'+str(i)+'_'+str(col),weight=row_weight)
		column_set=list(m.get_row(i))
		G.add_edge('edge_'+str(i)+'_'+str(column_set[0]),'edge_'+str(i)+'_'+str(column_set[1]),weight=row_weight)
	
	
	
	original_weight=m.get_weight()
	
	# Maximum weight matching algorithm called on the auxiliary graph.
	# This gives us A maximum matching, and, correspondingly, a maximum weight 2m-matching in the matrix
	# This DOES NOT maximize the number of adjacencies chosen.
	matching=nx.max_weight_matching(G)
	present=0

	new_matrix=bm.BinaryMatrix()
	discarded_rows=bm.BinaryMatrix()
	
	# Using the matching to find a component mCi1P respecting spanning subgraph of maximum weight.
	# The corresponding matrices- for the adjacencies retained, and those discarded- are the outputs.
	for vertex in matching.keys():
		if (not vertex[0:5]=='edge_' or not matching[vertex][0:5]=='edge_'):
			# Condition for retaining an adjacency: both ends of the corresponding adjacency must be matched, but not to each other.
			row=0
			first_end=None
			other_end=None
			if vertex[0:5]=='edge_':
				k=5
				while vertex[k]!='_':
					k+=1
				row=int(vertex[5:k])
				first_end=int(vertex[k+1:])
	                        for col in m.get_row_info(row)._set:
	                                if col!=first_end:
						other_end=col
			elif matching[vertex][0:5]=='edge_':
				k=5
				while matching[vertex][k]!='_':
					k+=1
				row=int(matching[vertex][5:k])

				first_end=int(matching[vertex][k+1:])

				for col in m.get_row_info(row)._set:
					if col!=first_end:
						other_end=col
			# Check if both ends are matched
			if  'edge_'+str(row)+'_'+str(other_end) in matching.keys():
				if m.get_row_info(row) not in new_matrix:
					new_matrix.add_row_info(m.get_row_info(row))
			else:
			# Case when one edge-vertex is matched to a copy of a column, but the other is unmatched.
				if m.get_row_info(row) not in discarded_rows:
					discarded_rows.add_row_info(m.get_row_info(row))
		elif (vertex[0:5]=='edge_' and  matching[vertex][0:5]=='edge_'):
			k=5
		# The two edge-vertices of an adjacency are matched to each other.
			while vertex[k]!='_':
				k+=1
			if m.get_row_info(int(vertex[5:k])) not in discarded_rows:
				discarded_rows.add_row_info(m.get_row_info(int(vertex[5:k])))
	 
	return [new_matrix,discarded_rows]
	
def linearization_undirected(m,multiplicity): # multiplicity is a dictionary, with markers and their reversals as keys.
	# Creating auxiliary adjacency graph

	# Step 1: Adding 2 vertices for each copy of a marker
	G=nx.Graph()
	for col in m.get_support():
		for k in xrange(0,2*multiplicity[col]):
			G.add_node(str(col)+'_'+str(k))
	
	
	# Step 2: Adding 2 vertices for each adjacency.
	# Also adding corresponding edges and giving them weights as given in the matrix
	
	for i in xrange(0,m._height):#row in m:
		row_weight=float(m.get_row_info(i)._weight)
		row_id=m.get_row_info(i)._id
		if len(m.get_row(i))>2:
			print 'Adjacencies only!'
			sys.exit(0)
		for col in m.get_row(i):
			G.add_node('edge_'+str(i)+'_'+str(col))
			for k in xrange(0,2*multiplicity[col]):
				G.add_edge(str(col)+'_'+str(k),'edge_'+str(i)+'_'+str(col),weight=row_weight)
		column_set=list(m.get_row(i))
		G.add_edge('edge_'+str(i)+'_'+str(column_set[0]),'edge_'+str(i)+'_'+str(column_set[1]),weight=row_weight)
	
	
	
	original_weight=m.get_weight()
	
	# Maximum weight matching algorithm called on the auxiliary graph.
	# This gives us A maximum matching, and, correspondingly, a maximum weight 2m-matching in the matrix
	# This DOES NOT maximize the number of adjacencies chosen.
	matching=nx.max_weight_matching(G)
	present=0

	new_matrix=bm.BinaryMatrix()
	discarded_rows=bm.BinaryMatrix()
	
	# Using the matching to find a component mCi1P respecting spanning subgraph of maximum weight.
	# The corresponding matrices- for the adjacencies retained, and those discarded- are the outputs.
	for vertex in matching.keys():
		if (not vertex[0:5]=='edge_' or not matching[vertex][0:5]=='edge_'):
			# Condition for retaining an adjacency: both ends of the corresponding adjacency must be matched, but not to each other.
			row=0
			first_end=None
			other_end=None
			if vertex[0:5]=='edge_':
				k=5
				while vertex[k]!='_':
					k+=1
				row=int(vertex[5:k])
				first_end=int(vertex[k+1:])
	                        for col in m.get_row_info(row)._set:
	                                if col!=first_end:
						other_end=col
			elif matching[vertex][0:5]=='edge_':
				k=5
				while matching[vertex][k]!='_':
					k+=1
				row=int(matching[vertex][5:k])
				first_end=int(matching[vertex][k+1:])
				for col in m.get_row_info(row)._set:
					if col!=first_end:
						other_end=col
			# Check if both ends are matched
			if  'edge_'+str(row)+'_'+str(other_end) in matching.keys():
				if m.get_row_info(row) not in new_matrix:
					new_matrix.add_row_info(m.get_row_info(row))
			else:
			# Case when one edge-vertex is matched to a copy of a column, but the other is unmatched.
				if m.get_row_info(row) not in discarded_rows:
					discarded_rows.add_row_info(m.get_row_info(row))
		elif (vertex[0:5]=='edge_' and  matching[vertex][0:5]=='edge_'):
			k=5
		# The two edge-vertices of an adjacency are matched to each other.
			while vertex[k]!='_':
				k+=1
			if m.get_row_info(int(vertex[5:k])) not in discarded_rows:
				discarded_rows.add_row_info(m.get_row_info(int(vertex[5:k])))
	 
	return [new_matrix,discarded_rows]
