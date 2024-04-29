
MAIN_DIR=$1 # Main project directory, with Rmds
OUTPUT_DIR=$2 # Output directory
DMRDIRECTION=$3 # Whether up- or down-regulated DMRs
RERUN_ADDFEATUREIDS=$4 # Whether we need to regenerate features IDs


# Loop over features
echo "##### WORKING ON STAGE DMRS \\\\\\\\\\\\\\\\"
for feature in Tss TssFlank Transcribed Enhancer Repressed
do

  # Report feature
  echo "##### WORKING ON :- $feature"

  # Run
  bash $MAIN_DIR/scripts/overlapDMRs.sh \
       $MAIN_DIR/results/dmr/dmrff_results.stage_as_continuous.treat_log2FC0.1.$DMRDIRECTION.tsv \
       $MAIN_DIR/results/dmr/ROADMAPstates/E049_15_coreMarks_dense.${feature}Only.bed \
       $OUTPUT_DIR \
       ${DMRDIRECTION}_dDMR $feature \
       True $RERUN_ADDFEATUREIDS \
       1

done

# Loop over features
echo "##### WORKING ON SEX DMRS \\\\\\\\\\\\\\\\"
for feature in Tss TssFlank Transcribed Enhancer Repressed
do

  # Report feature
  echo "##### WORKING ON :- $feature"

  # Run
  bash $MAIN_DIR/scripts/overlapDMRs.sh \
       $MAIN_DIR/results/dmr/dmrff_results.sex.$DMRDIRECTION.tsv \
       $MAIN_DIR/results/dmr/ROADMAPstates/E049_15_coreMarks_dense.${feature}Only.bed \
       $OUTPUT_DIR \
       ${DMRDIRECTION}_sDMR $feature \
       True False \
       1

done
