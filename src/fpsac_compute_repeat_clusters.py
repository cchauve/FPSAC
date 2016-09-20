# FPSAC: Fast Phylogenetic Scaffolding of Ancient Contigs.
# May 2013.
# Cedric Chauve: cedric.chauve@sfu.ca

# Computing maximal repeat clusters in an assembly graph

# argument 1 = Input  edges file (doubled)
# argument 2 = Input  copy numbers file (undoubled)
# argument 3 = Output repeat clusters file

import sys
import graph_tools

# command_line="python "
# for a in sys.argv:
#     command_line=command_line+' '+a
# print "------> COMMAND: "+command_line

copynumbers={}
graph={}
marker_adjacencies=[]
cn_file=open(sys.argv[2],"r").readlines()
i=0
while i<len(cn_file):
    undoubled_marker=int(cn_file[i].split('>')[1].split(' ')[0])
    copynumber=int(cn_file[i].split(' ')[1])
    copynumbers[2*undoubled_marker]=copynumber
    copynumbers[2*undoubled_marker-1]=copynumber
    marker_adjacencies.append((2*undoubled_marker-1,2*undoubled_marker))
    graph[2*undoubled_marker]=[]
    graph[2*undoubled_marker-1]=[]
    i=i+1

edges_file=open(sys.argv[1],"r").readlines()
i=0
while i<len(edges_file):
    e1=int(edges_file[i].split(':')[1].split(' ')[0])
    e2=int(edges_file[i].split(':')[1].split(' ')[1])
    if copynumbers[e1]>1 and copynumbers[e2]>1:
        graph[e1].append(e2)
        graph[e2].append(e1)
    i=i+1
for m in marker_adjacencies:
    if copynumbers[m[0]]>1:
        graph[m[0]].append(m[1])
        graph[m[1]].append(m[0])

comp = graph_tools.composantes(graph)

clusters=open(sys.argv[3],"w")
id=1
for c in comp:
    c.sort()
    if len(c)>1:
        clusters.write(str(id)+"|1;REPEAT_CLUSTER:")
        for m in c:
            clusters.write(str(m)+" ")
        id+=1
        clusters.write("\n")
