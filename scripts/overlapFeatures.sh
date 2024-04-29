
DMRDIRECTION=$1
PADJTHR=$2

# Loop over features
## Stage
for feature in Tss TssFlank transcribed enhancer repressed
do

  # Report feature
  echo "##### WORKING ON :- $feature"

  # Run
  bash ../../scripts/overlapDMR_ROADMAP.sh \
       dmrff_results.stage_as_continuous.treat_log2FC0.1.$DMRDIRECTION.tsv \
       ROADMAPstates/E049_15_coreMarks_dense.${feature}Only.bed \
       ${DMRDIRECTION}_dDMR chondroMSC_$feature \
       True False \
       $PADJTHR

done

# Loop over features
## Sex
for feature in Tss TssFlank transcribed enhancer repressed
do

  # Report feature
  echo "##### WORKING ON :- $feature"

  # Run
  bash ../../scripts/overlapDMR_ROADMAP.sh \
       dmrff_results.sex.$DMRDIRECTION.tsv \
       ROADMAPstates/E049_15_coreMarks_dense.${feature}Only.bed \
       ${DMRDIRECTION}_sDMR chondroMSC_$feature \
       True False \
       $PADJTHR

done
