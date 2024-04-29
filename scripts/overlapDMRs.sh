




#################### CONFIGURATION ####################

# DEFINE MAIN VARIABLES
AFILE=$1 # Main intervals of interest
BFILE=$2 # Reference intervals
OUTPUT_DIR=$3
DMRTYPE=$4 # What is A? This will be the created as an output directory and will be appended as the suffix of the DMR IDs
FEATUREID=$5 # What is B?
IDSFORA=$6 # Add IDs to A? True if so, else anything (ie False) for no
IDSFORB=$7 # Add IDs to B? True if so, else anything (ie False) for no
PVALTHR=$8

# Get file names without suffix
AFILENAME=$( basename $AFILE .tsv )
BFILENAME=$( basename $BFILE .bed )





######################## CODE #########################

# Convert regions to unix formatting, else AWK gets confused
# dos2unix $AFILE
# dos2unix $BFILE

# Make file
mkdir -p ${OUTPUT_DIR}/${DMRTYPE}/overlap${FEATUREID}/

# Make DMR file as .bed
## Non significance filtered
cat $AFILE | awk 'NR > 1 {print $3"\t"$4"\t"$5"\t"$1}' > ${OUTPUT_DIR}/${DMRTYPE}/$AFILENAME.withids.bed
## Significance filtered
cat $AFILE | awk -v PVALTHR=$PVALTHR '$14 < PVALTHR {print $3"\t"$4"\t"$5"\t"$1}' > /tmp/k && \
awk '{print $1"\t"$2"\t"$3"\t"$4}' /tmp/k > ${OUTPUT_DIR}/${DMRTYPE}/$AFILENAME.withids.pFilt${PVALTHR}.bed

# If we want to add IDs to the ATAC regions
if [ $IDSFORB = True ]
then

	# Add unique DMR IDs
	echo "#### ADDING IDs to ${FEATUREID} regions"
	paste <( awk '{print $1"\t"$2"\t"$3}' $BFILE ) \
     	      <( awk -v FEATURE="$FEATUREID" '{ for (i = 1; i <= NR; i++); print FEATURE"_"i-1 }' $BFILE ) \
      	      > ${OUTPUT_DIR}/$BFILENAME.withids.bed

fi

# Paramters for bedtools
# -wo 	  : Report overlaps in A and B
# -a, -b  : Input regions

# Overlap intervals
echo "#### Overlapping intervals"
bedtools intersect \
		  -wo \
		  -a ${OUTPUT_DIR}/${DMRTYPE}/$AFILENAME.withids.pFilt${PVALTHR}.bed \
		  -b ${OUTPUT_DIR}/$BFILENAME.withids.bed \
		  > ${OUTPUT_DIR}/${DMRTYPE}/overlap${FEATUREID}/${DMRTYPE}_${FEATUREID}_all_overlaps.with_coords.bed
# Drop coords for B, retain only ID
awk '{print $1"\t"$2"\t"$3"\t"$4";"$8"\t"$9}' \
		  ${OUTPUT_DIR}/${DMRTYPE}/overlap${FEATUREID}/${DMRTYPE}_${FEATUREID}_all_overlaps.with_coords.bed \
		  > ${OUTPUT_DIR}/${DMRTYPE}/overlap${FEATUREID}/${DMRTYPE}_${FEATUREID}_all_overlaps.bed

# Keep only unique DMRs that overlap
echo "#### Keeping only unique DMRs"
grep \
    -w -f \
    <( cat ${OUTPUT_DIR}/${DMRTYPE}/overlap${FEATUREID}/${DMRTYPE}_${FEATUREID}_all_overlaps.bed | awk '{print $4 }' | sed 's/[;].*//g' ) \
    ${OUTPUT_DIR}/${DMRTYPE}/$AFILENAME.withids.pFilt${PVALTHR}.bed |\
    uniq \
    > ${OUTPUT_DIR}/${DMRTYPE}/overlap${FEATUREID}/unique_${DMRTYPE}_${FEATUREID}_overlaps.bed

# Get details on ATAC-overlapping DMRs
grep \
    -w -f \
    <( cat ${OUTPUT_DIR}/${DMRTYPE}/overlap${FEATUREID}/unique_${DMRTYPE}_${FEATUREID}_overlaps.bed | awk '{print $4}' ) \
    $AFILE \
    > ${OUTPUT_DIR}/${DMRTYPE}/overlap${FEATUREID}/$AFILENAME.${FEATUREID}overlap.tsv
cat \
   <( head -n 1 $AFILE ) \
   ${OUTPUT_DIR}/${DMRTYPE}/overlap${FEATUREID}/$AFILENAME.${FEATUREID}overlap.tsv \
   > /tmp/k && mv /tmp/k \
   ${OUTPUT_DIR}/${DMRTYPE}/overlap${FEATUREID}/$AFILENAME.${FEATUREID}overlap.tsv

grep \
    -w -f \
    <( cat ${OUTPUT_DIR}/${DMRTYPE}/overlap${FEATUREID}/${DMRTYPE}_${FEATUREID}_all_overlaps.with_coords.bed | awk '{print $8}' ) \
    ${OUTPUT_DIR}/$BFILENAME.withids.bed \
    > ${OUTPUT_DIR}/${DMRTYPE}/overlap${FEATUREID}/$BFILENAME.${FEATUREID}overlap.bed

# Copy input to output folders
cp $AFILE ${OUTPUT_DIR}/${DMRTYPE}/$AFILENAME.tsv
cp ${OUTPUT_DIR}/$BFILENAME.withids.bed ${OUTPUT_DIR}/${DMRTYPE}/$BFILENAME.withids.bed
