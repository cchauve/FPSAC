# FPSAC: Fast Phylogenetic Scaffolding of Ancient Contigs.
# May 2013.
# Cedric Chauve: cedric.chauve@sfu.ca

# From a scaffold file with doubled markers to the same scaffold with undoubled oriented markers

# argument 1 = Input  doubled scaffold file
# argument 2 = Output undoubled scaffold file

import sys,string

# command_line="python "
# for a in sys.argv:
#     command_line=command_line+' '+a
# print "------> COMMAND: "+command_line

name_doubled=sys.argv[1]
name_halved=sys.argv[2]

input_file=open(name_doubled,"r").readlines()
output_file=open(name_halved,"w")

for i in range(len(input_file)):
	if input_file[i][0]==">" or input_file[i][0]=="#":
		output_file.write(input_file[i])
	else:
		mots=input_file[i].split()
		m=0
		while m<len(mots):
			if mots[m].find("_")>=0 or mots[m] == "T":
				output_file.write(mots[m]+" ")
				m=m+1
			else:
				m1=int(mots[m])
				try:
					m2=int(mots[m+1])
				except:
					if mots[m-1].find("C")>=0:
						m2=-100
					else:
						m=m+1

						continue
					#endif
				#endtry

				if abs(m1-m2)!=1:
					if mots[m-1].find("C")>=0:
						if m1%2==0:
							output_file.write("-"+str(m1/2)+" ")
						else:
							output_file.write(str((m1+1)/2)+" ")
						#endif

						m=m+1

						continue
					else:
						print("ERROR:  C1P_halve_PQRtree " + str(m1) + " " + str(m2))
						sys.exit(0)
					#endif
				#endif

				if m1<m2 and m2%2==0:
					output_file.write(str(m2/2)+" ")
				elif m2<m1 and m1%2==0:
					output_file.write("-"+str(m1/2)+" ")
				else:
					if mots[m-1].find("C")>=0:
						if m1%2==0:
							output_file.write("-"+str(m1/2)+" ")
						else:
							output_file.write(str((m1+1)/2)+" ")
						#endif

						m=m+1

						continue
					else:
						print("ERROR:  C1P_halve_PQRtree " + str(m1) + " " + str(m2))
						sys.exit(0)
					#endif
				#endif

				m=m+2
			#endif
		#endwhile

		output_file.write("\n")
	#endif
#endwhile

output_file.close()
