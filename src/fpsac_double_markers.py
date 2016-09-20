# FPSAC: Fast Phylogenetic Scaffolding of Ancient Contigs.
# May 2013.
# Cedric Chauve: cedric.chauve@sfu.ca

# Doubling markers
 
# argument 1 = Input  undoubled families occurrences file
# argument 2 = Output doubled families occurrences file

import sys,string

# command_line="python "
# for a in sys.argv:
#     command_line=command_line+' '+a
# print "------> COMMAND: "+command_line

families_file=open(sys.argv[1],"r").readlines()
output_file=open(sys.argv[2],"w")

i=0
number=1
while i<len(families_file):
    if families_file[i][0]=='>':
        id=families_file[i][1:families_file[i].find(" ")]
	comment=families_file[i][families_file[i].find(" ")+1:]	
        marker=[]
        j=i
        i=i+1
        while i<len(families_file) and len(families_file[i])>1:
            spe=families_file[i][:families_file[i].find(".")]
            chr=families_file[i][families_file[i].find(".")+1:families_file[i].find(":")]
            start=int(families_file[i][families_file[i].find(":")+1:families_file[i].find("-")])
            stop=int(families_file[i][families_file[i].find("-")+1:families_file[i].find(" ")])
            dir=families_file[i][families_file[i].find(" ")+1:families_file[i].find(" ")+2]
            marker.append([spe,chr,start,stop,dir])
            i=i+1
        output_file.write(">"+str(2*int(id)-1)+" "+comment)
        for m in marker:
            if m[4]=="+":
                output_file.write(m[0]+"."+m[1]+":"+str(m[2])+"-"+str(m[2]+1)+" +\n")
            else:
                output_file.write(m[0]+"."+m[1]+":"+str(m[3]-1)+"-"+str(m[3])+" +\n")
        output_file.write("\n")        
        output_file.write(">"+str(2*int(id))+" "+comment)
        for m in marker:
            if m[4]=="+":
                output_file.write(m[0]+"."+m[1]+":"+str(m[3]-1)+"-"+str(m[3])+" +\n")
            else:
                output_file.write(m[0]+"."+m[1]+":"+str(m[2])+"-"+str(m[2]+1)+" +\n")
        output_file.write("\n")        
        i=i+1
        number=number+2

