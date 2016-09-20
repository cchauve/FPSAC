# FPSAC: Fast Phylogenetic Scaffolding of Ancient Contigs.
# May 2013.
# Cedric Chauve: cedric.chauve@sfu.ca

# Computing adjacencies and repeat-spanning intervals according to a Dollo parsimony criterion

# argument 1 = Input  markers file, for  doubled markers
# argument 2 = Input  profile file for undoubled markers families
# argument 3 = Input  species list
# argument 4 = Input  species pairs
# argument 5 = Input  0/1 flag for circular genomes
# argument 6 = Input  type of syntenies (ADJ, RSI)
# argument 7 = Output conserved syntenies file

import sys,copy

# command_line="python "
# for a in sys.argv:
#     command_line=command_line+' '+a
# print "------> COMMAND: "+command_line

# ----- Reading species ----------------------------------------------------------
specieslist = []
speciesfile = open(sys.argv[3],"r").readlines()
i=0
while i<len(speciesfile):
    if speciesfile[i][0] != "#":
        specieslist.append(speciesfile[i].rstrip('\n'))
    i=i+1
        
# ----- Reading species pairs----------------------------------------------------------
speciespairslist = []
speciespairsfile = open(sys.argv[4],"r").readlines()
i=0
while i<len(speciespairsfile):
    if speciespairsfile[i][0] != "#":
        species1=speciespairsfile[i].split(' ')[0]
        species2=speciespairsfile[i].split(' ')[1].rstrip('\n')
        speciespairslist.append((species1,species2))
    i=i+1

# ----- Reading profiles ----------------------------------------------------------
copynumbers={}
profilesfile = open(sys.argv[2],"r").readlines()
i=0
while i<len(profilesfile):
    undoubled_marker=int(profilesfile[i].split('>')[1].split(' ')[0])
    copynumber=int(profilesfile[i].split('>')[1].split(' ')[1])
    copynumbers[2*undoubled_marker]=copynumber
    copynumbers[2*undoubled_marker-1]=copynumber
    i=i+1

# ----- Reading circular Flag----------------------------------------------------------

if sys.argv[5]=="0":
    flag_circ=False
else:
    flag_circ=True
    
# ----- Reading markers ----------------------------------------------------------

# defines a marker
class Marker:
	def __init__(self, ident, species, ch, st, end):
		self._id = ident
		self._species = species
		self._ch = ch
		self._st = st
		self._end = end
	#enddef	
	def __lt__(self, a):
		return self._st < a._st
	#enddef	
	def __le__(self, a):
		return self._st <= a._st
	#enddef
#endclass

# extracts the markers from a markers file
# return: markers - dict of list of list of Marker
def extract_markers(markers_file):
    markers = {}	# markers
    markersfile = open(markers_file,"r").readlines()  # the file of the markers
    sp_chm={}	# species and chromosome dictionary mapping chromosomes to indices #the chromosomes in the markers
    marker_list=[]  # the read markers
    i=0		# file index
    id = ""	# current marker id
	
    while i<len(markersfile):
        if markersfile[i][0]=='>':
            id=markersfile[i][1:-1].split(' ')[0]
            i=i+1
            marker_list.append(int(id))
            while i<len(markersfile) and len(markersfile[i])>1:
                spe=markersfile[i][:markersfile[i].find(".")]		# current species
                chr=markersfile[i][markersfile[i].find(".")+1:markersfile[i].find(":")]		# current chromosome
                if not spe in markers:
                    markers[spe] = []
                    sp_chm[spe]={}
                start=int(markersfile[i][markersfile[i].find(":")+1:markersfile[i].find("-")])		# start pos
                stop=int(markersfile[i][markersfile[i].find("-")+1:markersfile[i].find(" ")])		# end pos
                if(sp_chm[spe].has_key(chr) == False):
                    if sp_chm[spe]=={}:
                        sp_chm[spe][chr] = 0
                    else:
                        sp_chm[spe][chr] = max(sp_chm[spe].values())+1
                    if (sp_chm[spe][chr] >= len(markers[spe])):
                        for k in range(sp_chm[spe][chr]-len(markers[spe])+1):
                            markers[spe].append([])
                markers[spe][sp_chm[spe][chr]].append(Marker(int(id), spe, chr, int(start), int(stop)))
                i=i+1
        i+=1
    for i in markers:
        for j in range(len(markers[i])):
            markers[i][j].sort(key=lambda m:m._st)	# sort markers by starting positions on the chromosome. markers are treated as linear for now
    return markers,marker_list

markers, marker_list = extract_markers(sys.argv[1])

# ----- Auxiliary functions ----------------------------------------------------------
# def join_chromosomes(blocks, spe, markers,circ):
#     extended_markers = copy.copy(markers)
#     seq = []
#     fake = max(markers) + 1
#     copynumbers[fake]=0
#     for ch in blocks[spe]:
#         if(len(ch) > 0):
#             seq += [m._id for m in ch]
#             if flag_circ:
#                 seq += [m._id for m in ch]                
#             seq.append(fake)
#     return seq

# Tests if two doubled markers originate from the same undoubled one
# Assumption: marker 1 < marker2
def test_marker_adj(marker1, marker2):
    if marker1%2!=0 and marker2==marker1+1:
        return True
    else:
        return False

# ----- Detecting adjacencies ----------------------------------------------------------
# Computing the list of pairs (adjacencies,species)
# Adjacencies might not be conserved: all non-marker adjacencies are recorded

# Modify here: need to remember marker positions in chromosomes.
# Suggested solution: store Marker data structures in sequence rather than marker ids.
# Problem: Too much space used.
# Alternatives? Suggest we directly compute adjacencies, without creating sequence. 

if sys.argv[6]=="ADJ":
    offset=1
    adjacencieslist=[]
    for species in specieslist:
	fake=Marker(max(marker_list)+offset,species,'dummy_ch',0,0)
	is_first=0
	for chr in markers[species]:
            if len(chr)>0:
                fake._ch=chr[0]._ch
                if flag_circ:
                    chr.extend(chr) # extend chromosome by a copy of itself if circular.
                for pos,marker1 in enumerate(chr[:-1]): # need to optimize if circular. if so, only need to go to len(chr)<<1+1
                    if pos==0 and is_first>0 and not (flag_circ):
                        adjacencieslist.append((str(marker1._id)+' '+str(fake._id),species,marker1._ch,str(marker1._st-1),str(marker1._st-1),'-'))
                    marker2=chr[pos+1]
                    gap_sequence=species+'.'+chr[pos]._ch+':'+str(marker1._end+1)+'-'+str(marker2._st-1)
                    if marker1._id>marker2._id and not test_marker_adj(marker2._id,marker1._id):
                        sgn="-" # The adjacency has been seen reversed
                        adjacencieslist.append((str(marker2._id)+" "+str(marker1._id),species,gap_sequence,sgn))
                    elif marker1._id<marker2._id and not test_marker_adj(marker1._id,marker2._id):
                        sgn="+" # The adjacency has been seen as written in the output file
                        adjacencieslist.append((str(marker1._id)+" "+str(marker2._id),species,gap_sequence,sgn))
            is_first+=1
	offset+=1
    adjacencieslist.sort()
    syntenieslist=adjacencieslist

# ----- Detecting repeat spanning intervals----------------------------------------------
# Computing the list of pairs (intervals,species)    
elif sys.argv[6]=="RSI":
    fake = Marker(max(marker_list)+1,'dummy','dummy',0,0)
    intervalslist=[]
    for species in specieslist:
	for chr in markers[species]:
            if flag_circ:
                chr.extend(chr)		# if the chromosome is circular, extend it with another copy to detenct rsi
            if len(chr)>0:
                start=-1
                for pos,marker1 in enumerate(chr[:-1]):
                    marker2=chr[pos+1]
                    if copynumbers[marker1._id]==1 and copynumbers[marker2._id]>1:
                        start = pos
                    if copynumbers[marker1._id]>1 and copynumbers[marker2._id]==1:
                        end=pos+1
                        if pos>1 and start>-1 and chr[start]._id > chr[end]._id:
                            p=end
                            interval=''
                            gap_sequence=species+'.'+chr[p]._ch+':'
                            while p>start:
                                interval = interval+' '+str(chr[p]._id)
                                p-=1
                            while p<end:
                                gap_sequence = gap_sequence+str(chr[p]._end+1)+'-'+str(chr[p+1]._st-1)+','	# add string containing gap information
                                p+=1
                            gap_sequence=gap_sequence.rstrip(',')
                            interval=interval+' '+str(chr[start]._id)
                            interval=interval.lstrip()
                            intervalslist.append((interval,species,gap_sequence,'-')) # The interval has been seen reversed
                        elif pos>1 and start>-1 and chr[start]._id < chr[end]._id:
                            p=start
                            interval=''
                            gap_sequence=species+'.'+chr[p]._ch+':'
                            while p<end:
                                interval = interval+' '+str(chr[p]._id)
                                gap_sequence = gap_sequence+str(chr[p]._end+1)+'-'+str(chr[p+1]._st-1)+','	# add string containing gap information
                                p+=1
                            gap_sequence=gap_sequence.rstrip(',')
                            interval=interval+' '+str(chr[end]._id)
                            interval=interval.lstrip()
                            intervalslist.append((interval,species,gap_sequence,'+'))  # The interval has been seen as written in the output file
    intervalslist.sort()
    syntenieslist=intervalslist

# ----- Filtering and writing----------------------------------------------
# Filtering syntenies not conserved in informative species pairs, fusing and writing similar syntenies
output=open(sys.argv[7],"w")
previoussynteny=""
previousgap=""
foundspecies=set([])
id=-1
for i in syntenieslist:
    synteny=i[0]
    species=i[1]
    if synteny != previoussynteny:
        if id>-1:
            toconserve=0
            for sp in speciespairslist:
                if sp[0] in foundspecies and sp[1] in foundspecies:
                    toconserve=1
            if toconserve==1:
                stringsynteny=stringsynteny+','.join(foundspecies)+":"+previoussynteny+' #'.join(gap_sequence)+"\n"
                output.write(stringsynteny)
        id=id+1
        foundspecies=set([])
        stringsynteny=str(id)+"|0;"#+species
        gap_sequence=[' #'+i[2]+' '+i[3]]
    elif i[2]!=previousgap:
        gap_sequence.append(i[2]+' '+i[3])
    foundspecies.add(species)
    previoussynteny=synteny
    previousgap=i[2]

# Last synteny
toconserve=0
for sp in speciespairslist:
    if sp[0] in foundspecies and sp[1] in foundspecies:
        toconserve=1
if toconserve==1:
    stringsynteny=stringsynteny+','.join(foundspecies)+":"+previoussynteny+' #'.join(gap_sequence)+"\n"
    output.write(stringsynteny)
output.close()
