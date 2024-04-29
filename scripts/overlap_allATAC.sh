
MAIN_DIR=$1 # Main project directory, with Rmds
OUTPUT_DIR=$2 # Output for files
DMRDIRECTION=$3 # Up- or down-regulated DMRs
RERUN_ADDFEATUREIDS=$4 # Whether to regenerate non-DMR features (ATAC etc)

# Run for stage/developmental
echo "##### WORKING ON STAGE DMRS"
bash $MAIN_DIR/scripts/overlapDMRs.sh \
     $MAIN_DIR/results/dmr/dmrff_results.stage_as_continuous.treat_log2FC0.1.$DMRDIRECTION.tsv \
     $MAIN_DIR/results/dmr/foetal_knee_hg19.bed \
     $OUTPUT_DIR \
     ${DMRDIRECTION}_dDMR ATAC \
     True $RERUN_ADDFEATUREIDS \
     0.001

# Run for sex
echo "##### WORKING ON SEX DMRS"
bash $MAIN_DIR/scripts/overlapDMRs.sh \
     $MAIN_DIR/results/dmr/dmrff_results.sex.$DMRDIRECTION.tsv \
     $MAIN_DIR/results/dmr/foetal_knee_hg19.bed \
     $OUTPUT_DIR \
     ${DMRDIRECTION}_sDMR ATAC \
     True False \
     0.01
