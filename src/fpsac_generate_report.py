# FPSAC: Fast Phylogenetic Scaffolding of Ancient Contigs.
# May 2013.
# Cedric Chauve: cedric.chauve@sfu.ca

# Generate a report on a run of the FPSAC pipeline

import sys

report_file=open(sys.argv[1],"w")

# CONTIGS
ctgs_length_file=open(sys.argv[2],"r").readlines()
nb_ctgs=0
lg_ctgs=0
lg_ctgs_100=0
lg_ctgs_500=0
lg_ctgs_1000=0
lg_ctgs_5000=0
nb_ctgs_100=0
nb_ctgs_500=0
nb_ctgs_1000=0
nb_ctgs_5000=0
for ctg in ctgs_length_file:
    nb_ctgs+=1
    lg=int(ctg.split(' ')[2].rstrip( ))
    lg_ctgs+=lg
    if lg>=100:
        nb_ctgs_100+=1
        lg_ctgs_100+=lg
    if lg>=500:
        nb_ctgs_500+=1
        lg_ctgs_500+=lg
    if lg>=1000:
        nb_ctgs_1000+=1
        lg_ctgs_1000+=lg
    if lg>=5000:
        nb_ctgs_5000+=1
        lg_ctgs_5000+=lg
report_file.write("CONTIGS\n")
report_file.write("\tall\t\t"+str(nb_ctgs)+"\t"+str(lg_ctgs)+"nt\n")
report_file.write("\t>=100nt\t\t"+str(nb_ctgs_100)+"\t"+str(lg_ctgs_100)+"nt\n")
report_file.write("\t>=500nt\t\t"+str(nb_ctgs_500)+"\t"+str(lg_ctgs_500)+"nt\n")
report_file.write("\t>=1000nt\t"+str(nb_ctgs_1000)+"\t"+str(lg_ctgs_1000)+"nt\n")
report_file.write("\t>=5000nt\t"+str(nb_ctgs_5000)+"\t"+str(lg_ctgs_5000)+"nt\n")

# EXTANT GENOMES
genomes_length_file=open(sys.argv[3],"r").readlines()
genomes=[]
lg_genomes={}
kar_genomes={} # karyotype
coverage_genomes={}
nb_genomes=0
for extant in genomes_length_file:
    molecule=extant.split(' ')[0].split('>')[1]
    genome=molecule.split('.')[0]
    lg=int(extant.split(' ')[2].rstrip( ))
    if genome not in genomes:
        genomes.append(genome)
        lg_genomes[genome]=lg
        nb_genomes+=1
        kar_genomes[genome]=1
        coverage_genomes[genome]=0
    else:
        lg_genomes[genome]+=lg
        kar_genomes[genome]+=1
lg_genomes_all=0
kar_genomes_all=0
for g in genomes:
    lg_genomes_all+=lg_genomes[g]
    kar_genomes_all+=kar_genomes[g]
report_file.write("EXTANT GENOMES\n")
report_file.write("\t"+str(nb_genomes)+" genomes\t"+str(kar_genomes_all)+" chrs/plasmids\t"+str(lg_genomes_all)+"nt\n")
for g in genomes:
    report_file.write("\t"+g+"\t\t\t"+str(kar_genomes[g])+" chrs/plasmids\t"+str(lg_genomes[g])+"nt\n")

# FAMILIES AND HITS
hits_file=open(sys.argv[4],"r").readlines()
families_all_file=open(sys.argv[5],"r").readlines()
families_mult_file=open(sys.argv[6],"r").readlines()
nb_mult_fam=0
copies_mult_fam=0
fam_mult={}
for f in families_mult_file:
    fam=f.split(' ')[0].split('>')[1]
    mult=int(f.split(' ')[1].rstrip( ))
    copies_mult_fam+=mult
    fam_mult[fam]=mult
    if mult>1:
        nb_mult_fam+=1
families_filtered_file=open(sys.argv[7],"r").readlines()
filtered_fams=[]
for f in families_filtered_file:
    if f[0]==">":
        fam=f.split(' ')[0].split('>')[1].rstrip( )
        filtered_fams.append(fam)
families_length_file=open(sys.argv[8],"r").readlines()
nb_filtered_fam=len(families_length_file)-len(filtered_fams)
lg_fam=0
lg_fam_filtered=0
for f in families_length_file:
    fam=f.split(' ')[0].split('>')[1]
    lg=int(f.split(' ')[1].rstrip( ))
    lg_fam+=fam_mult[fam]*lg
    if fam in filtered_fams:
        lg_fam_filtered+=fam_mult[fam]*lg
report_file.write("HOMOLOGOUS FAMILIES\n")
report_file.write("\tcontigs/genomes hits\t"+str(len(hits_file))+"\n")
report_file.write("\tall\t\t"+str(len(families_length_file))+"\tmult.>1 "+str(nb_mult_fam)+"\t"+str(lg_fam)+"nt\n")
report_file.write("\tfiltered\t"+str(len(families_length_file)-nb_filtered_fam)+"\tmult.>1 "+str(nb_mult_fam-nb_filtered_fam)+"\t"+str(lg_fam_filtered)+"nt\n")

for f in families_filtered_file:
    if f[0]==">":
        fam=f.split(' ')[0].split('>')[1].rstrip( )
    elif len(f)>1:
        genome=f.split('.')[0]
        start=int(f.split(':')[1].split('-')[0])
        end=int(f.split(':')[1].split('-')[1].split(' ')[0])
        coverage_genomes[genome]+=end-start+1
        
report_file.write("COVERAGE OF EXTANT GENOMES BY HOMOLOGOUS FAMILIES\n")
for g in genomes:
    report_file.write("\t"+g+"\t\t\t"+str(coverage_genomes[g])+"nt/"+str(lg_genomes[g])+"nt\t"+str(float(coverage_genomes[g])/(float(lg_genomes[g])))+"\n")

# ANCESTRAL ADJACENCIES
adj_file=open(sys.argv[9],"r").readlines()
kept_adj_file=open(sys.argv[10],"r").readlines()
disc_adj_file=open(sys.argv[11],"r").readlines()
adj_05=0
adj_1=0
kadj_05=0
kadj_1=0
for a in adj_file:
    w=float(a.split('|')[1].split(';')[0])
    if w>=0.5:
        adj_05+=1
    if w>=1:
        adj_1+=1
for a in kept_adj_file:
    w=float(a.split('|')[1].split(';')[0])
    if w>=0.5:
        kadj_05+=1
    if w>=1:
        kadj_1+=1
dadj_05=adj_05-kadj_05
dadj_1=adj_1-kadj_1
report_file.write("ANCESTRAL ADJACENCIES\n")
report_file.write("\tall\t\t"+str(len(adj_file))+"\tweight>=0.5 "+str(adj_05)+"\tweight=1 "+str(adj_1)+"\n")
report_file.write("\tkept\t\t"+str(len(kept_adj_file))+"\tweight>=0.5 "+str(kadj_05)+"\tweight=1 "+str(kadj_1)+"\n")
report_file.write("\tdiscarded\t"+str(len(disc_adj_file))+"\tweight>=0.5 "+str(dadj_05)+"\tweight=1 "+str(dadj_1)+"\n")

# REPEAT CLUSTERS
rep_clust_file=open(sys.argv[12],"r").readlines()
report_file.write("REPEAT CLUSTERS\n\t"+str(len(rep_clust_file))+"\n")

# REPEAT SPANNING INTERVALS
rsi_file=open(sys.argv[13],"r").readlines()
kept_rsi_file=open(sys.argv[14],"r").readlines()
disc_rsi_file=open(sys.argv[15],"r").readlines()
rsi_05=0
rsi_1=0
krsi_05=0
krsi_1=0
for a in rsi_file:
    w=float(a.split('|')[1].split(';')[0])
    if w>=0.5:
        rsi_05+=1
    if w>=1:
        rsi_1+=1
for a in kept_rsi_file:
    w=float(a.split('|')[1].split(';')[0])
    if w>=0.5:
        krsi_05+=1
    if w>=1:
        krsi_1+=1
drsi_05=rsi_05-krsi_05
drsi_1=rsi_1-krsi_1
report_file.write("REPEAT SPANNING INTERVALS\n")
report_file.write("\tall\t\t"+str(len(rsi_file))+"\tweight>=0.5 "+str(rsi_05)+"\tweight=1 "+str(rsi_1)+"\n")
report_file.write("\tkept\t\t"+str(len(kept_rsi_file))+"\tweight>=0.5 "+str(krsi_05)+"\tweight=1 "+str(krsi_1)+"\n")
report_file.write("\tdiscarded\t"+str(len(disc_rsi_file))+"\tweight>=0.5 "+str(drsi_05)+"\tweight=1 "+str(drsi_1)+"\n")

# SCAFFOLD
cars_file=open(sys.argv[19],"r").readlines()
nb_cars=0
nb_circ_cars=0
for c in cars_file:
    if c[0]=="_":
        nb_cars+=1
    if c[0]=="_C":
        nb_circ_cars+=1
nb_lin_cars=nb_cars-nb_circ_cars
unassigned_adjacencies=open(sys.argv[20],"r").readlines()
report_file.write("SCAFFOLDS\n\tnumber "+str(nb_cars)+"\tcircular "+str(nb_circ_cars)+"\tlinear "+str(nb_lin_cars)+"\tunassigned adjacencies "+str(len(unassigned_adjacencies))+"\n")


# CARS EXTREMITIES ADJACENCIES
caradj_file=open(sys.argv[16],"r").readlines()
#kept_caradj_file=open(sys.argv[17],"r").readlines()
#disc_caradj_file=open(sys.argv[18],"r").readlines()
caradj_05=0
#kcaradj_05=0
nb_caradj=0
#nb_kcaradj=0
#nb_dcaradj=0
for a in caradj_file:
    w=float(a.split('|')[1].split(';')[0])
    if w<1:
        nb_caradj+=1
    if w>=0.5 and w<1:
        caradj_05+=1
#for a in kept_caradj_file:
#    w=float(a.split('|')[1].split(';')[0])
#    if w<1:
#        nb_kcaradj+=1
#    if w>=0.5 and w<1:
#        kcaradj_05+=1
#dcaradj_05=caradj_05-kcaradj_05
#nb_dcaradj= nb_caradj- nb_kcaradj
report_file.write("SCAFFOLD EXTREMITIES ADJACENCIES\n")
report_file.write("\tall\t\t"+str(nb_caradj)+"\tweight>=0.5 "+str(caradj_05)+"\n")
#report_file.write("\tkept\t\t"+str(nb_kcaradj)+"\tweight>=0.5 "+str(kcaradj_05)+"\n")
#report_file.write("\tdiscarded\t"+str(nb_dcaradj)+"\tweight>=0.5 "+str(dcaradj_05)+"\n")

# GAPS
gaps_file=open(sys.argv[18],"r").readlines()
nb_supp_gaps=0
gap_diff=0
nb_gap_diff=0
nb_gap_diff_10=0
gap_diff_10=0
nb_unsupp_gaps=0
for g in gaps_file:
    if g[0]==">":
        gap=g.rstrip( ).split(' ')
        if gap[-2]=="supported_length_interval":
            nb_supp_gaps+=1
            interval=gap[-1].split('-')
            int1=int(interval[0])
            int2=int(interval[1])
            gap_diff+=int2-int1
            if int2-int1>0:
                nb_gap_diff+=1
            if int2-int1>10:
                nb_gap_diff_10+=1
                gap_diff_10+=int2-int1
        else:
            nb_unsupp_gaps+=1
report_file.write("INTER MARKERS GAPS\n")       
report_file.write("\twith supported length interval\t\t"+str(nb_supp_gaps)+"\tambig. "+str(nb_gap_diff)+" "+str(gap_diff)+"nt\tambig.>10nt "+str(nb_gap_diff_10)+" "+str(gap_diff_10)+"nt\n")
report_file.write("\twith unsupported length interval\t"+str(nb_unsupp_gaps)+"\n")

#MAP
lg_map=0
lg_ctgs_map=0
lg_gaps_map=0
map_file=open(sys.argv[17],"r").readlines()
for l in map_file:
    entry=l.split(' ')
    if len(entry)>0 and entry[0]=="MARKER":
        start=int(entry[2].split('-')[0])
        end=int(entry[2].split('-')[1])
        lg_map+=end-start+1
        lg_ctgs_map+=end-start+1
    elif len(entry)>0 and entry[0]=="GAP":
        lg=int(entry[2])
        lg_map+=lg
        lg_gaps_map+=lg
report_file.write("ANCESTRAL SEQUENCE\n")
report_file.write("\tlength "+str(lg_map)+"nt\tfrom contigs "+str(lg_ctgs_map)+"nt\tfrom gaps "+str(lg_gaps_map)+"nt\n")
