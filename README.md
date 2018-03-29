# FPSAC, Fast Phylogenetic Scaffolding of Ancient Contigs.
## June 2013. 
## Contact: Cedric Chauve (Dept. Mathematics, Simon Fraser University), cedric.chauve@sfu.ca

##1. WHAT IS FPSAC?

FPSAC is a suite of shell and python scripts used to scaffold ancient
contigs, for an ancestral genome of interest, using a phylogenetic
comparative approach.

The methodological principles of FPSAC are described in the paper
"FPSAC: Fast Phylogenetic Scaffolding of Ancient Contigs", by Ashok
Rajaraman, Eric Tannier and Cedric Chauve, currently in revision for
Bioinformatics. More details are provided in section 4 below.

## 2. INSTALLING FPSAC?

FPSAC is composed of a set of Python scripts, located in the src
directory, together with a master shell script, fpsac.sh. So there is
no compilation required.

FPSAC requires the following softwares to be available:  
makembindex (from the BLAST distribution 2.2.27+ and above)  
makeblastdb (from the BLAST distribution 2.2.27+ and above)  
blastn      (from the BLAST distribution 2.2.27+ and above)  
muscle      (version v3.8.31 and above)

FPSAC has been developed and tested on Unix (including Linux and
MacOS) systems, using Python (http://www.python.org/) version 2.7.3,
with the numpy library (http://numpy.scipy.org/).

## 3. USING FPSAC.

FPSAC takes as input   
(1) ancient_contigs.fa: a FASTA file of ancient contigs (ancestral
genome),  
(2) extant_genomes.fa: a FASTA file of extant genomes,  
(3) species_tree: a species tree describing the relationship between
the extant genomes and with a marked ancestral node indicating trhe
position of the ancestor of interest,  
(4) three parameters s, l, c that will be described in section 4
below.  
 
l is an integer  
s is an integer between 0 and 100  
c is an integer greater than 1  

FPSAC can be run with the following line command:

./fpsac.sh output_directory_name species_tree ancient_contigs.fa extant_genomes.fa l s c ancestor_name

Contraints on the input data are the following:  
- species names have to be the same in the extant genomes file and in the species tree  
- the species tree needs to be in simple nexus format with no comments and branch length  
- the species tree needs to be fully resolved but for the outgroups  
- the species names should not contain any of the symbols '.', ':', '-'  

The results of running FPSAC as above are available in the directory  
DIR=output_directory_name/results_ancestor_name
  
The important results files are the following:  
- output_directory_name/data			   
  formatted data  
- $DIR/contigs/families_with_contig_names   
  homologous families  
- $DIR/scaffold/scaffold_order  
  order of ancestral markers along the scaffolds  
- $DIR/scaffold/scaffold.fasta  
  scaffold sequence, with gaps filled by Ns  
- $DIR/finishing/ancestral_sequence.fasta  
  scaffold sequence with gaps filled by estimated sequences  
- $DIR/finishing/ancestral_sequence_map  
  map of the scaffold detailling the link between ancestral segments  
  and extant segments  
- output_directory_name/report_ancestor_name  
  report on the FPSAC run  

Many other intermediate files are generated that are available in
$DIR/contigs/, $DIR/scaffold/, $DIR/finishing/, as described in the
following section.

An example is provided using the Black Deat agent data set, see
directory black_death and file run_black_death.sh.

## 4. FPSAC METHODOLOGY.
 
FPSAC relates to a generic scheme for reconstructing ancestral genomes
organization and is composed of four phases:

- Computing homologous families. A homologous family is composed of at
  least one contig segment (ancestral marker) and several extant
  genomes segments (extant markers), that pairwise align, with high
  similarity, along their whole length. Each homologous family is
  assigned a multiplicity bounding the number of occurrences (copy
  number) of ancestral marker(s) from this family in the ancestral
  genome.

- Computing putative ancestral adjacencies. An ancestral adjacency is
  composed of two ancestral markers that are believed to have been
  contiguous in the ancestral genome. We predict them using a Dollo
  parsimony principle. All adjacencies are weighted according to their
  phylogenetic conservation, defining a weighted adjacency graph.

- Scaffolding from ancestral adjacencies. If the set of all
  ancestral adjacencies is not compatible with a circular chromosomal
  structure that respects the multiplicity constraints of homologous
  families, we compute a maximum weight subset of adjacencies that is
  compatible with such a circular chromosomal structure.  Next, as
  adjacencies alone can define several contig orders, due to repeated
  ancestral markers forming tangles in the adjacency graph, conserved
  intervals spanning repeats are used to clear the ambiguities, in a
  ancestral markers forming tangles in the adjacency graph, conserved
  intervals spanning repeats are used to clear the ambiguities, in a
  way similar to the use of mate-pairs to scaffold extant genomes.

- Estimating inter-markers gap lengths and sequences. For each
  ancestral adjacency, the length of the ancestral gap between the two
  involved markers is estimated, again using a Dollo parsimony
  criterion, from the length of the gap between the corresponding
  extant adjacencies (extant gaps). The sequences of the extant gaps
  whose length agrees with the estimated ancestral gap length are
  aligned into a multiple sequence alignment, that is used to
  reconstruct a putative ancestral gap sequence.

### 4.1. Computing homologous markers families.  
 
We map the ancient contigs onto the extant genomes using megablast
with default parameters.  Every significant hit, defined here by a
length of at least l nucleotides and a sequence identity above s%,
where l and s are two of the three parameters entered in the comand
line, indicates two homologous sequences, one located on a contig, and
one located on an extant genome. Due to rearrangements and repeats,
some contigs do not align over their whole length to every extant
genome, indicating potential evolutionary breakpoints. In order to
detect families of homologous segments, we apply an iterative
segmentation procedure, which produces contig and extant genome
segments such that (1) contig segments align over their whole length
to extant genomes segments, and (2) pairs of extant genome segments do
not overlap ( i.e. either they have the same coordinates, or they are
completely disjoint). Contigs/genomes hits are iteratively segmented,
with segments below length l being discarded as they appear, until
both (1) and (2) are statisfied. Then all aligned sequences can be
clustered into families of highly similar ancient and extant
sequences.

Next, we assign to each homologous family a multiplicity, that is the
expected number of occurrences of the ancestral marker of the family
in the ancestral genome.  The multiplicity of a family is computed
from the number of occurrences of the extant markers in the extant
genomes the family profile, in order to minimize the number of
evolutionary gain/loss along the branches of the considered
phylogenetic tree. It is computed by a linear time dynamic programming
algorithm.

Finally, families are filtered to discard all families of multiplicity
above c, the third parameter passed to the master script fpsac.sh.


The files produced by this stage are in the directory $DIR/contigs and
are the following:

- megablast_hits: raw hits obtained with megablast  
- ancient_extant_hits: hits between ancient contigs and extant genomes
  the format of this file is the following: each line describes a hit
  extant_genome_name chromosome start end orientation(1/-1) contig_name start end orientation (1/-1)

- families_ancestral_content: ancestral multiplicity of each
  homologous family  
- families_profiles: extant profile (multiplicity) of each family

- families.fasta: ancestral sequence associated to each family  
- families_length: length of the above sequences

- families_with_contig_names: coordinates of the segments that define
  families, both in the extant genomes and in the ancient contigs.

- families_filtered: coordinates of families whose ancestral
  mutliplicity is at most c.

- families_filtered_doubled: filtered families coordinates after
  doubling each marker by introducing one doubled marker for each
  undoubled marker extremity (for family I, the doubled markers are
  2*I-1 and 2*I)

### 4.2. Scaffolding.
 
All files obtained during the scaffolding stage are in $DIR/scaffold.

Predicting ancestral adjacencies.  To account for the orientation of
markers in predicted ancestral syntenic features (adjacencies and
intervals), we decompose each marker (ancestral or extant) into two
marker extremities, its head and its tail, and strore these doubled
families in $DIR/contigs/families_filtered_doubled.

Adjacencies are then defined in terms of marker extremities instead of
markers, and are computed following a Dollo parsimony principle: two
ancestral marker extremities form an ancestral adjacency if they are
contiguous (no other marker is between them in the chromosome) in at
least two extant genomes whose evolutionary path contains the ancestor
of interest.  Adjacencies are weighted according to their patterns of
phylogenetic conservation. The weighted adjacency graph is defined as
follows: its vertices are the markers extremities and its edges are
the weighted adjacencies.

The adjacencies are stored in the files

- adjacencies: unweighted adjacencies  
- adjacencies_weighted: weighted adjacencies  

The format of the adjacencies file is as follows:  
id|weight;species_containing_the_adjacency:marker1 marker2 list_of_coordinates_of_extant_gaps_between_marker1_and_marker2

Circularization into ancestral scaffolds.  An ancestral scaffold is a
linear or circular order of ancestral markers.  The set of ancestral
adjacencies might not translate into an unambiguous set of ancestral
scaffolds for (1) there might not exist a set of circular or linear
markers orders that contain all adjacencies and respects the
multiplicity of each marker, and (2) even if ancestral adjacencies can
be organized in ancestral scaffolds it can be ambiguous, i.e.  several
dsets of scaffolds can exist, because of marker multiplicities.

To address point (1), we compute a maximum weight subset of ancestral
adjacencies such that every marker extremity belongs to a number of
adjacencies that is at most the multiplicity of the marker family.
Such a selected subset of ancestral adjacencies is compatible with an
order of the markers into a set of linear and/or circular scaffolds
that respects the copy number constraint given by the ancestral marker
multiplicities. 

The results of this stage are in the files  
- adjacencies_kept: adjacencies kept for thre scaffolding  
- adjacencies_discarded: adjacencies removed from the adjacency graph  

The maximum weight subset is computed using a maximum weight matching
algorithm implemented in the NetworkX python library:
http://networkx.github.io/.

However, there can be many markers orders compatible with the
adjacency set, due to markers with multiplicity >1 creating tangles in
the adjacency graph, leading to ambiguities regarding the structure of
the considered genome. To address point (2), we rely on small
conserved intervals that span markers with multiplicity >1, called
repeat spanning intervals. Repeat spanning intervals are weighted
using the same method than ancestral adjacencies, and for each repeat
cluster, we greedily select repeat spanning intervals that both are
compatible with the adjacencies selected during the previous step and
containing markers of the cluster and satisfy the multiplicity
constraints of the markers.

Repeat spanning intervals (RSI)  are in the files  
- repeat_spanning_intervals: all RSI, unweighted  
- repeat_spanning_intervals_weighted: all RSI, weighted  
- repeat_spanning_intervals_kept: RSI kept for clearing ambiguities in the scaffolds  
- repeat_spanning_intervals_discarded: RSI discarded by the greedy heuristic  
The repeat clusters are stored in the file  
- repeat_clusters  

All adjacencies, repeat spanning intervals, repeat clusters are
described with the doubled markers.

Provided all repeats are spanned by enough conserved intervals, this
results into an unambiguous scaffolding that includes all ancestral
markers, including repeated ones. Otherwise adjacencies composed of
two repeats and that do not belong to a repeat spanning interval are
discarded, resulting in a more fragmented, but unambiguous, set of
scaffolds.

The resulting markers orders that define the backbone of the scaffolds
are available in two files:  
- scaffold_order_doubled: orders defined on doubled markers  
- scaffold_order: orders on undoubled markers  

In the scaffold files, scaffolds are called CAR (for Contiguous
Ancestral Region), and a CAR starting by _Q and ending by Q_ is a
linear CAR, while a CAR starting by _C and ending by C_ is a circular
CAR.

### 4.3. Estimating inter-contig gaps lengths and sequences.
 
An ancestral gap in the ancestral marker ordering is the sequence
located between two consecutive ancestral markers. For each ancestral
gap, we consider the extant genomes in which the occurrences of the
corresponding markers are consecutive and in the same respective
orientations as in the CAR, thus defining an extant gap.  We define a
conserved extant gap as an extant gap whose length is equal in two
extant genomes whose evolutionary path contains the ancestor,
following again a Dollo criterion. The lengths of conserved extant
gaps define an ancestral gap length interval. If there is no conserved
extant gap, the ancestral gap length interval is defined by the
minimum and maximum extant gaps length between the markers in the
extant genomes.

The ancestral gaps, together with their corresponding extant gaps
coordinates are in the file  
- $DIR/scaffold/gaps_coordinates  
The same, but with the gap length interval are in the file  
- $DIR/scaffold/gaps_coordinates_and_length

Next, for each ancestral gap, we align all sequences of extant gaps
whose length is in this interval into a multiple sequence alignment,
obtained with the software muscle. 

The sequences of the extant gaps are the files   
- *fasta_extant in the directory $DIR/finishing/alignments  
The muscle alignments are in the files  
- *fasta_muscle in the directory $DIR/finishing/alignments

A parsimonious estimation of each ancestral gap sequence is obtained
from the corresponding alignment of extant gaps sequences using the
classical Fitch algorithm. 

The inferred ancestral sequences are in the files   
- *fasta_ancestral in the directory $DIR/finishing/alignments

Finally, we combine all sequences from contigs and gaps into a set of
scaffolds sequences. There are two scaffold sequence files  
- $DIR/scaffold/scaffold.fasta: gaps are filled by Ns (number chosen
  as the midpoint of the gap length interval)  
- $DIR/finishing/ancestral_sequence.fasta: gaps replaced by the
  inferred sequences

There is also a map of the ancetsral sequence, that describes the
extant genome segments that support every ancetsral sequence segment,
in the file  
- $DIR/finishing/ancestral_sequence_map

### 4.4. Refining scaffolding with low support adjacencies.
 
It is possible to refine the scaffolding by considering all pairs of
markers located at scaffolds extremities and looking for support
(however not under a Dollo criterion) in extant genomes for such
potential adjacencies. This is done by default by FPSAC, and
availmable in the file  
- $DIR/scaffold/adjacencies_CARS  
These adjacencies ar not weighted.

Next, after inspection, if a subset of such CARs extremities
adjacencies is deemed reliable, as it was on the Black Death data set,
they can be integrated to the scaffolding using the script
./fpsac_extra_1.sh.

For the Black Death data set, for example, three of the four CARs
extremities were deemed reliable ands were recorded in the file
black_death/final/results_Black_Death_Agent_8291/scaffold/adjacencies_CARS_chosen
and then integrated to the previous scaffolding by the command

./fpsac_extra_1.sh \  
   black_death/final \  
   black_death/final/results_Black_Death_Agent_8291/scaffold/adjacencies_CARS_chosen \  
   Black_Death_Agent_8291

The resulting files, in the directories $DIR/scaffold and
$DIR/finishing are all suffixed by "_CARS_combined". 
The joined scaffolds, in the file
$DIR/finishing/ancestral_sequence_CARS_combined.fasta are separated by
a sequence of 50 Ns.

