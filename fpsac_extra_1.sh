# FPSAC, Fast Phylogenetic Scaffolding of Ancient Contigs.
# June 2013. 
# cedric.chauve@sfu.ca

# Integrating a set of chosen extremities between CARs/scaffolds

# PARAMETERS
#
# $1 = directory where data and results are stored, directory name (created if required)
# $2 = chosen CARS extremities file (doubled markers)
# $3 = ancestor name, single word text

DIR=$1
DATA=$DIR/input
RESULTS=$DIR/results_$3
CONTIGS=$RESULTS/contigs
SCAFFOLD=$RESULTS/scaffold
FINISHING=$RESULTS/finishing
FINISHING_ALG=$FINISHING/alignments
ANCESTOR_NAME=$3

ANCESTRAL_CONTIGS_FA=$DATA/ancestral_contigs.fasta 

FAMILIES_COORDINATES_ALL=$CONTIGS/families_with_contig_names
FAMILIES_FA=$CONTIGS/families.fasta
FAMILIES_LG=$CONTIGS/families_length

ADJACENCIES_CARS=$2
GAPS_COORDINATES_LG=$SCAFFOLD/gaps_coordinates_and_length
GAPS_COORDINATES_LG_COMBINED=$SCAFFOLD/gaps_coordinates_and_length_CARS_combined
SCAFFOLD_ORDER_DOUBLED=$SCAFFOLD/scaffold_order_doubled
SCAFFOLD_ORDER_DOUBLED_COMBINED=$SCAFFOLD/scaffold_order_doubled_CARS_combined
SCAFFOLD_ORDER_COMBINED=$SCAFFOLD/scaffold_order_CARS_combined
SCAFFOLD_FA_COMBINED=$SCAFFOLD/scaffold_CARS_combined.fasta
SCAFFOLD_FINISHED_FA_COMBINED=$FINISHING/ancestral_sequence_CARS_combined.fasta
SCAFFOLD_MAP_COMBINED=$FINISHING/ancestral_sequence_map_CARS_combined
RSI_KEPT=$SCAFFOLD/repeat_spanning_intervals_kept
FAMILIES_ANCESTRAL_CONTENT=$CONTIGS/families_ancestral_content

echo "---> Integrating CARS adjacencies with ketp adjacencies and gaps"
python src/fpsac_integrate_CARS_extremities_adjacencies.py \
    ${ADJACENCIES_CARS} \
    ${SCAFFOLD_ORDER_DOUBLED} \
    ${GAPS_COORDINATES_LG} \
    ${FINISHING_ALG}/ \
    ${SCAFFOLD_ORDER_DOUBLED_COMBINED} \
    ${GAPS_COORDINATES_LG_COMBINED}
# Output: $SCAFFOLD_ORDER_DOUBLED_COMBINED - scaffold order with doubled markers
# Output: $GAPS_COORDINATES_LG_COMBINED - gaps with extra gaps added

echo "---> Undoubling markers"
python src/fpsac_halve_scaffold.py \
    ${SCAFFOLD_ORDER_DOUBLED_COMBINED} \
    ${SCAFFOLD_ORDER_COMBINED}
# Output: $SCAFFOLD_ORDER_COMBINED - scaffold order

echo "---> Computing scaffold sequence with Ns in gaps"
python src/fpsac_compute_scaffold_sequence_N.py \
    ${FAMILIES_FA} \
    ${SCAFFOLD_ORDER_COMBINED} \
    ${GAPS_COORDINATES_LG_COMBINED} \
    ${SCAFFOLD_FA_COMBINED}
# Output: $SCAFFOLD_FA_COMBINED - scaffold sequences with Ns in gaps

echo "---> Computing ancestral sequence with gaps"
python src/fpsac_compute_scaffold_sequence.py \
    ${FAMILIES_FA} \
    ${SCAFFOLD_ORDER_COMBINED} \
    ${GAPS_COORDINATES_LG_COMBINED} \
    ${FINISHING_ALG}/ \
    ${SCAFFOLD_FINISHED_FA_COMBINED}
# Output: $SCAFFOLD_FINISHED_FA_COMBINED - FASTA sequence of the ancestral scaffold, including gaps

echo "---> Computing ancestral map"
python src/fpsac_compute_scaffold_map.py \
    ${FAMILIES_COORDINATES_ALL} \
    ${FAMILIES_FA} \
    ${SCAFFOLD_ORDER_COMBINED} \
    ${GAPS_COORDINATES_LG_COMBINED} \
    ${FINISHING_ALG}/ \
    ${SCAFFOLD_MAP_COMBINED}
# Output: $SCAFFOLD_MAP_COMBINED - map of the ancestral scaffold
