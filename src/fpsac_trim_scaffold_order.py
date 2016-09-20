# FPSAC: Fast Phylogenetic Scaffolding of Ancient Contigs.
# May 2013.
# Cedric Chauve: cedric.chauve@sfu.ca

# Trim CARs to remove extremities composed of repeats. This is intended to deal with the case where
# a repeat r forms a repeat cluster and is adjacent to to non-repeats u and v but there is no repeat
# spanning interval u.r.v, in which case the repeat r appears two times instead of once. This point needs
# to be improved in the next version.

import sys

# argument 1 = CARs file
# argument 2 = copy numbers
# argument 3 = trimmed CARs
# argument 4 = trimmed extremities

# command_line="python "
# for a in sys.argv:
#     command_line=command_line+' '+a
# print "------> COMMAND: "+command_line

# Read profiles
copynumbers={}
mult_file=open(sys.argv[2]).readlines()
for line in mult_file:
    marker=int(line[1:].split(' ')[0])
    multiplicity=int(line[1:].split(' ')[1])
    copynumbers[2*marker]=multiplicity
    copynumbers[2*marker-1]=multiplicity

# Read scaffold
CARS=open(sys.argv[1]).readlines()

# Process CARs
strips=[]
new_CARS=open(sys.argv[3],'w')
trimmed_CAR=open(sys.argv[4],'w')
for car in CARS:
    # Read reconstructed species name or CAR name
    if car[0]!='_':
        new_CARS.write(car)
        trimmed_CAR.write(car)
    else:
        car=car.rstrip( ).split(' ')
        if (car[0]=="_Q"): # Nothing need to be done if we have a circular CAR
            stop=0
            i=-2
            trimmed=''
            while stop==0:
                if copynumbers[int(car[i])]==1: # If we encounter a marker of multiplicity 1, add it as a reference to the trimmed string and terminate
                    #trimmed=car[i]+' '+trimmed  # Comment out if we do not need a reference
                    stop+=1  # Terminate loop
                else:
                    trimmed=' '.join(car[i-1:i+1])+' '+trimmed # For copy number>1, add it and its implicity pair to trimmed string
                    car[i-1]=car[-1]  # Set last marker to Q_
                    car=car[:i]      # Truncate CAR
            trimmed_CAR.write(trimmed+'\n')
        new_CARS.write(' '.join(car)+'\n')
        
new_CARS.close()
trimmed_CAR.close()
