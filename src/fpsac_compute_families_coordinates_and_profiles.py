# FPSAC: Fast Phylogenetic Scaffolding of Ancient Contigs.
# May 2013.
# Cedric Chauve: cedric.chauve@sfu.ca

# Computing families of orthologous markers from ancient contigs/extant genomes alignments hits

# argument 1 = Input  hits file
# argument 2 = Output homologous families occurrences file
# argument 3 = Output homologous families profile file
# argument 4 = Input  length threshold for filtering occurrences and families
# argument 5 = Input  dissimilarity threshold for clustering hits

import sys,string,graph_tools,time

# command_line="python "
# for a in sys.argv:
#     command_line=command_line+' '+a
# print "------> COMMAND: "+command_line

time0=time.time()

hit_file = open(sys.argv[1]).readlines()
family_file = open(sys.argv[2],"w")
profile_file = open(sys.argv[3],"w")
FILTER_HITS = int(sys.argv[4]) #100
FILTER_FAMILIES = int(FILTER_HITS + 0.1*FILTER_HITS)
ERROR = (100.0-float(sys.argv[5]))/100.0 #0.05

TRACE_LEVEL = 0

def trace(seq,level):
	if level <= TRACE_LEVEL:
		print seq

#########################
##### 1/ READING THE HITS
#########################

#########  READ FUNCTION

# READS A HIT FILE AND FILLS A HIT TABLE IN WHICH EACH HIT IS
# [species,chromosome,start,stop,direction     # EXTANT SEGMENT
#  contig_name,start,stop]         # CONTIG 
def read_hits(hit_file):       
    hits = []
    for ligne in range(len(hit_file)):
        mots = hit_file[ligne].split()
        if mots[4] == "1":
            direction = "+"
        else:
            direction = "-"
        hits.append([mots[0],mots[1],int(mots[2]),int(mots[3]),direction,
			     mots[5],int(mots[6]),int(mots[7])])
    return hits


#########  OUTPUT FUNCTION

def write_hits(hits,file_name):
    sortie = open(file_name,"w")
    for h in hits:
        sortie.write(h[0]+" "+h[1]+" "+str(h[2])+" "+str(h[3])+" "+str(h[4])+" "+
			  h[5]+" "+str(h[6])+" "+str(h[7])+" "+str(h[8])+"\n")



###########################################################
##### 2/ CUTTING THE HITS
###########################################################

#########  COMPARISON FUNCTIONS

def compare_hits_chromosome(x,y):
    if x[0] != y[0]:
        return cmp(x[0],y[0])
    elif  x[1] != y[1]:
        return cmp(x[1],y[1])
    elif  x[2] != y[2]:
        return cmp(x[2],y[2])
    else:
        return cmp(y[3],x[3])

def compare_hits_contig(x,y):
    if x[5] != y[5]:
        return cmp(x[5],y[5])
    elif  x[6] != y[6]:
        return cmp(x[6],y[6])
    elif x[7] != y[7]:
        return cmp(y[7],x[7])
    elif x[0] != y[0]:
        return cmp(x[0],y[0])
    elif  x[1] != y[1]:
        return cmp(x[1],y[1])
    elif  x[2] != y[2]:
        return cmp(x[2],y[2])
    else:
        return cmp(y[3],x[3])

def oppose(x):
    if x == "1":
        return "-1"
    elif x == "+":
        return "-"
    elif x == "-1":
        return "1"
    elif x == "-":
        return "+"
    elif x == 1:
        return -1
    elif x == -1:
        return 1




# DO TWO HITS IN THE LIST OVERLAP ON CONTIGS (STRICT OVERLAP, RETURN FALSE IF EQUAL OR DISJOINT)
#~ def intersect_contig(hits):
    #~ result = False
    #~ hits.sort(lambda x,y: compare_hits_contig(x,y))
    #~ h = 0
    #~ while h < len(hits)-1 and not result:
        #~ if (hits[h][5] == hits[h+1][5] and # same name
            #~ (hits[h][6] != hits[h+1][6] or hits[h][7] != hits[h+1][7]) and # different coordinates
            #~ hits[h][7] >= hits[h+1][6]): # the second begins before the end of the first
            #~ result = True
        #~ else:
            #~ h = h + 1
    #~ return result


# RETURNS THE AMOUNT OF INTERSECTION ON CHROMOSOME GIVEN
# THE TWO HITS ARE ON THE SAME CHROMOSOME
def intersection(hits,a,b):
    return (max(0,(min(hits[a][3],hits[b][3])-max(hits[a][2],hits[b][2])+1)/
                         float(max(hits[a][3],hits[b][3])-min(hits[a][2],hits[b][2])+1)))

# DO TWO HITS IN THE LIST OVERLAP ON THE CHROMOSOMES? (FLEXIBLE OVERLAP, WITH ERROR)
def intersect_chromosome(hits,error):
    result = False
    hits.sort(lambda x,y: compare_hits_chromosome(x,y))
    h = 0
    while h < len(hits) and not result:
        h1 = h + 1
        while (h1 < len(hits) and
               hits[h][0] == hits[h1][0] and
               hits[h][1] == hits[h1][1] and
               hits[h][3] >= hits[h1][2]):
            inter = intersection(hits,h,h1)
            result = (inter > error and inter < 1-error)
            h1 = h1 + 1
        h = h + 1
    return result


#########  STAT FUNCTION FOR DEVELOPMENT
#~ def write_intersections(hits,file_name):
    #~ hits.sort(lambda x,y: compare_hits_chromosome(x,y))
    #~ sortie = open(file_name,"w")
    #~ h = 0
    #~ while h < len(hits):
        #~ mindesmax = hits[h][3]
        #~ maxdesmin = hits[h][2]
        #~ h1 = h + 1
        #~ while (h1 < len(hits) and
               #~ hits[h][0] == hits[h1][0] and
               #~ hits[h][1] == hits[h1][1] and
               #~ hits[h1][2] <= hits[h][3]):
            #~ mindesmax = min(mindesmax,hits[h1][3])
            #~ maxdesmin = max(maxdesmin,hits[h1][2])
            #~ sortie.write(str(float(min(hits[h][3],hits[h1][3])-
                                   #~ max(hits[h][2],hits[h1][2])+1)/
                             #~ max(hits[h][3]-hits[h][2]+1,hits[h1][3]-hits[h1][2]+1))+"\n")
            #~ h1 = h1 + 1
        #~ if maxdesmin > mindesmax:
            #~ print hits[h]
    
        #~ h = h + 1

########  FILTER FUNCTIONS

# RETURNS ALL HITS IN THE LIST WITCH ARE BIGGER THAN MIN_SIZE
def filter(hits,MIN_SIZE):
    result = []
    for h in range(len(hits)):
        if hits[h][3] - hits[h][2] >= MIN_SIZE and hits[h][7] - hits[h][6] >= MIN_SIZE:
            result.append(list(hits[h]))
    return result
    

########  CUT FUNCTIONS

# NEW_HITS ARE THE HITS CUT WHEN INTERSECTIONS ON CONTIGS
def cut_along_contigs(hits):
    new_hits = []
    hits.sort(lambda x,y: compare_hits_contig(x,y))
    h = 0
    a_decouper = []
    debut = 0
    while h < len(hits):
        if (not hits[h][6] in a_decouper):
            a_decouper.append(hits[h][6])
        if (not hits[h][7]+1 in a_decouper):
            a_decouper.append(hits[h][7]+1)
        if h == len(hits) - 1 or hits[h][5] != hits[h+1][5]:
            a_decouper.sort()
            i = 0
            while i < len(a_decouper) - 1:
                segment = [a_decouper[i],a_decouper[i+1]-1]
                j = debut
                while j <= h:
                    if min(segment[1],hits[j][7]) >= max(segment[0],hits[j][6]):
                        direction = hits[j][4]
                        proportion = float(hits[j][3]-hits[j][2]+1)/(hits[j][7]-hits[j][6]+1)
                        if direction == "1" or direction == "+":
                            start = hits[j][2]+int(round((segment[0]-hits[j][6])*proportion))
                            stop = hits[j][3]+int(round((segment[1]-hits[j][7])*proportion))
                        else:
                            start = hits[j][2]-int(round((segment[1]-hits[j][7])*proportion))
                            stop = hits[j][3]-int(round((segment[0]-hits[j][6])*proportion))
                        new_hits.append([hits[j][0],hits[j][1],min(start,stop),max(stop,start),direction,
						    hits[j][5],segment[0],segment[1]])
                    j = j + 1
                i = i + 1
            a_decouper = []
            debut = h + 1
        h = h + 1
    return new_hits

def cut_along_chromosomes(hits):
    new_hits = []
    hits.sort(lambda x,y: compare_hits_chromosome(x,y))
    h = 0
    while h < len(hits):
	species = hits[h][0]
	chromosome = hits[h][1]
	direction = hits[h][4]

	a_decouper = [hits[h][2],hits[h][3]+1]

	h1 = h + 1
	while (h1 < len(hits) and
	       hits[h1][0] == species and
	       hits[h1][1] == chromosome and
	       hits[h1][2] <= hits[h][3]):
	    if (not hits[h1][2] in a_decouper):
		a_decouper.append(hits[h1][2])
	    if (hits[h1][3] < hits[h][3]) and (not hits[h1][3]+1 in a_decouper):
		a_decouper.append(hits[h1][3]+1)
	    h1 = h1 + 1
	    
	h1 = h - 1
	while (h1 >= 0 and
	       hits[h1][0] == species and
	       hits[h1][1] == chromosome and
	       hits[h1][3] >= hits[h][2]):
	    if (hits[h1][3] <= hits[h][3]) and (not hits[h1][3]+1 in a_decouper):
		a_decouper.append(hits[h1][3]+1)
	    if (hits[h1][2] > hits[h][2]) and (not hits[h1][2] in a_decouper):
		a_decouper.append(hits[h1][2])
	    h1 = h1 - 1

	a_decouper.sort()
	i = 0
	while i < len(a_decouper) - 1:
	    segment = [a_decouper[i],a_decouper[i+1]-1]
	    proportion = float(hits[h][3]-hits[h][2]+1)/(hits[h][7]-hits[h][6]+1)
	    if direction == "1" or direction == "+":
		new_hits.append([species,chromosome,segment[0],segment[1],direction,
				 hits[h][5],
				 hits[h][6]+int(round((segment[0]-hits[h][2])/proportion)),
				 hits[h][7]+int(round((segment[1]-hits[h][3])/proportion))])
	    else:
		new_hits.append([species,chromosome,segment[0],segment[1],direction,
				 hits[h][5],
				 hits[h][6]-int(round((segment[1]-hits[h][3])/proportion)),
				 hits[h][7]-int(round((segment[0]-hits[h][2])/proportion))])                
	    i = i + 1
			      
	h = h + 1
    return new_hits


###########################################################
##### MAIN PROCEDURE FOR SEGMENTATION
###########################################################


hits = read_hits(hit_file)
trace(str(len(hits))+" hits in the input",0)

hits = cut_along_contigs(hits)
#~ print "intersection on contigs",intersect_contig(hits)
#~ print len(hits),"hits"
#~ print "intersection on chromosomes",intersect_chromosome(hits,ERROR)
hits = filter(hits,1)
numero = 1
again = intersect_chromosome(hits,ERROR)
while again:
    
    print "iteration",len(hits),"hits"
    hits = cut_along_chromosomes(hits)
    hits = filter(hits,1)

    #~ write_intersections(hits,"intersections")

    hits = cut_along_contigs(hits)
    hits = filter(hits,FILTER_HITS)

    again = intersect_chromosome(hits,ERROR)

    numero = numero + 1
    

trace(str(len(hits))+" hits after cutting",0)



# ###########################################################
# ##### 3/ CLUSTERING HITS INTO FAMILIES
# ###########################################################


####3.1/ CONSTRUCTING THE GRAPH

# 3.1.1 VERTICES
graph = {}
hits.sort(lambda x,y: compare_hits_chromosome(x,y))
for h in range(len(hits)):
    hits[h].append(h)   #  CAREFUL : THE INDEX OF THE HIT WITH SORTING ALONG CHROMOSOME IS ADDED
    graph[h] = []

###3.1.2 EDGES CORRESPONDING TO INTERSECTIONS ON CHROMOSOMES
    
hits.sort(lambda x,y: compare_hits_chromosome(x,y))
h = 0
while h < len(hits):
#    sys.stdout.write("\r"+str(h)+"  in "+str(int(time.time()-time0))+" sec")
    h1 = h + 1
    while (h1 < len(hits) and
           hits[h][0] == hits[h1][0] and
           hits[h][1] == hits[h1][1] and
           hits[h][3] >= hits[h1][2]):
        inter = intersection(hits,h,h1)
        if inter > 1-ERROR:
            graph[h].append(h1)
            graph[h1].append(h)
            hits[h1][2] = hits[h][2]
            hits[h1][3] = hits[h][3]
        else:
            temp = hits[h][3]
            hits[h1][2] = hits[h][3]+1
        h1 = h1 + 1

    h = h + 1

#~ print intersect_contig(hits),intersect_chromosome(hits,0)

#####3.1.3 EDGES CORRESPONDING TO INTERSECTIONS ON CONTIGS

hits.sort(lambda x,y: compare_hits_contig(x,y))
h = 0
while h < len(hits):
#    sys.stdout.write("\r"+str(h)+"  in "+str(int(time.time()-time0))+" sec")
    h1 = h + 1
    while (h1 < len(hits) and
           hits[h][5] == hits[h1][5] and
           hits[h1][6] <= hits[h][7]):
        graph[hits[h][-1]].append(hits[h1][-1])
        graph[hits[h1][-1]].append(hits[h][-1])
        h1 = h1 + 1

    h = h + 1



##### 3.2 COMPONENTS OF THE GRAPH

comp = graph_tools.composantes(graph)
#~ print len(comp),"familles"

#### FILTERING THE FAMILIES

hits.sort(lambda x,y: compare_hits_chromosome(x,y)) # IMPORTANT: COMP REFERS TO INDICES OF SORTED HITS
c = 0
while c < len(comp):
	sizes = 5000000
	for i in comp[c]:
		sizes = min(sizes,hits[i][3]-hits[i][2])
	if sizes < FILTER_FAMILIES:
		del comp[c]
	else:
		c = c + 1

print len(comp),"familles"


##### 3.3 COVERAGE statistics
hits.sort(lambda x,y: compare_hits_chromosome(x,y))
coverage = 0
rep = {}
for c in comp:
	for h in c:
		rep[hits[h][5]] = 0
		#~ if ((hits[h][0] == "Yersinia_pestis_CO92") and
			#~ (h == len(hits) - 1 or
			#~ hits[h][0] != hits[h+1][0] or
			#~ hits[h][1] != hits[h+1][1] or 
			#~ hits[h][2] != hits[h+1][2] or
			#~ hits[h][3] != hits[h+1][3])):
			#~ coverage = coverage + hits[h][3] - hits[h][2] + 1

print "number of represented contigs",len(rep.keys())

####  WRITING FAMILIES (COMPLICATED FOR THE DIRECTION)

hits.sort(lambda x,y: compare_hits_chromosome(x,y)) # IMPORTANT: COMP REFER TO INDICES OF SORTED HITS
for c in range(len(comp)):
    trace("family "+str(c)+" over"+str(len(comp)),1)
    family_file.write(">"+str(c+1)+"\n")
    component = comp[c]
    component.sort(lambda x,y: compare_hits_chromosome(hits[x],hits[y]))
    graph_contig = {}
    for i in component:
        graph_contig[hits[i][5]] = []
        trace(hits[i],1)
    ref_contig = hits[component[0]][5]
    trace(ref_contig,1)
    
    # group hits having the same coordinates on the chromosome
    # current_hit : [species,chromosome,start,stop,
    #                       direction,contig_name,start,stop,
    #			      ....
    #                       direction,contig_name,start,stop]]
    
    i = 0
    while i < len(component):
        h = component[i]
	trace("a l'interieur de la composante "+str(hits[h]),1)
        if i > 0:
            hmoins = component[i-1]
        if (i == 0 or
            hits[h][0] != hits[hmoins][0] or
            hits[h][1] != hits[hmoins][1] or
            hits[h][2] != hits[hmoins][2] or 
            hits[h][3] != hits[hmoins][3]):
            current_hit = list(hits[h][:8])
        else:
            current_hit = current_hit + list(hits[h][4:8])
            if current_hit[4] == current_hit[-4]:
                if [current_hit[-3],1] not in graph_contig[current_hit[5]]:
                    graph_contig[current_hit[5]].append([current_hit[-3],1])
                if [current_hit[5],1] not in graph_contig[current_hit[-3]]:
                    graph_contig[current_hit[-3]].append([current_hit[5],1])
            else:
                if [current_hit[-3],-1] not in graph_contig[current_hit[5]]:
                    graph_contig[current_hit[5]].append([current_hit[-3],-1])
                if [current_hit[5],-1] not in graph_contig[current_hit[-3]]:
                    graph_contig[current_hit[-3]].append([current_hit[5],-1])
	i = i + 1
	
    i = 0
    while i < len(component):
        h = component[i]
	if i > 0:
            hmoins = component[i-1]
        if i < len(component) - 1:
            hplus = component[i+1]
        if (i == 0 or
            hits[h][0] != hits[hmoins][0] or
            hits[h][1] != hits[hmoins][1] or
            hits[h][2] != hits[hmoins][2] or 
            hits[h][3] != hits[hmoins][3]):
            current_hit = list(hits[h][:8])
        else:
            current_hit = current_hit + list(hits[h][4:8])
   	
                
        if (i == len(component) - 1 or
            hits[h][0] != hits[hplus][0] or
            hits[h][1] != hits[hplus][1] or
            hits[h][2] != hits[hplus][2] or 
            hits[h][3] != hits[hplus][3]):
            family_file.write(current_hit[0]+"."+current_hit[1]+":"+
                         str(current_hit[2])+"-"+str(current_hit[3])+" ")
            if current_hit[5] == ref_contig:
                family_file.write(current_hit[4])
            elif [current_hit[5],1] in graph_contig[ref_contig]:
                family_file.write(current_hit[4])
            elif [current_hit[5],-1] in graph_contig[ref_contig]:
                family_file.write(oppose(current_hit[4]))
            else:
                print "BUGBUGBUG",graph_contig[ref_contig],current_hit[5]
            family_file.write(" "+current_hit[5])
	    if not current_hit[5] == ref_contig and [current_hit[5],-1] in graph_contig[ref_contig]:
		    family_file.write("(rev)")
	    family_file.write(":"+str(current_hit[6])+"-"+str(current_hit[7]))
            if len(current_hit) > 8:
                family_file.write(","+current_hit[9])
		if not current_hit[9] == ref_contig and [current_hit[9],-1] in graph_contig[ref_contig]:
			family_file.write("(rev)")
		family_file.write(":"+str(current_hit[10])+"-"+str(current_hit[11]))
            family_file.write("\n")
            
        i = i + 1
    family_file.write("\n")


##### 3.4 PROFILES

for c in range(len(comp)):
    profile_file.write(">"+str(c+1)+"\n")
    local_hits = []
    especes = {}
    for h in comp[c]:
        if [hits[h][0],hits[h][1],hits[h][2],hits[h][3]] not in local_hits:
            local_hits.append([hits[h][0],hits[h][1],hits[h][2],hits[h][3]])
            if especes.has_key(hits[h][0]):
                especes[hits[h][0]] = especes[hits[h][0]] + 1
            else:
                especes[hits[h][0]] = 1
    for e in especes.keys():
        profile_file.write(e+" "+str(especes[e])+"\n")
    profile_file.write("\n")

