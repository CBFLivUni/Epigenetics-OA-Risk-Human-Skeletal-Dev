




#################### CONFIGURATION ####################

# DEFINE MAIN VARIABLES
AFILE=$1 # Main intervals of interest
BFILE=$2 # Reference intervals
DMRTYPE=$3 # What is A? This will be the created as an output directory and will be appended as the suffix of the DMR IDs
FEATUREID=$4 # What is B?
IDSFORA=$5 # Add IDs to A? True if so, else anything (ie False) for no
IDSFORB=$6 # Add IDs to B? True if so, else anything (ie False) for no
PVALTHR=$7

# Get file names without suffix
AFILENAME=$( basename $AFILE .tsv )
BFILENAME=$( basename $BFILE .bed )





######################## CODE #########################

# Convert regions to unix formatting, else AWK gets confused
# dos2unix $AFILE
# dos2unix $BFILE

# Make output
mkdir -p ${DMRTYPE}/overlap_${FEATUREID}

# Make DMR file as .bed
## Non significance filtered
cat $AFILE | awk 'NR > 1 {print $3"\t"$4"\t"$5"\t"$1}' > ${DMRTYPE}/$AFILENAME.withids.bed
## Significance filtered
cat $AFILE | awk -v PVALTHR=$PVALTHR '$14 < PVALTHR {print $3"\t"$4"\t"$5"\t"$1}' > /tmp/k && \
awk '{print $1"\t"$2"\t"$3"\t"$4}' /tmp/k > ${DMRTYPE}/$AFILENAME.withids.pFilt${PVALTHR}.bed

# If we want to add IDs to the ATAC regions
if [ $IDSFORB = True ]
then

	# Add unique DMR IDs
	echo "#### ADDING IDs to ${FEATUREID} regions"
	paste <( awk 'NR > 1 {print $1"\t"$2"\t"$3}' $BFILE ) \
     	      <( awk -v FEATURE="$FEATUREID" 'NR > 1 { for (i = 1; i < NR; i++); print FEATURE"_"i-1 }' $BFILE ) \
      	      > $BFILENAME.withids.bed

fi

# Paramters for bedtools
# -wo 	  : Report overlaps in A and B
# -a, -b  : Input regions

# Overlap intervals
echo "#### Overlapping intervals"
bedtools intersect \
		  -wo \
		  -a ${DMRTYPE}/$AFILENAME.withids.pFilt${PVALTHR}.bed \
		  -b $BFILENAME.withids.bed |\
    		  sed -z 's/\t_/\t/g' \
		  > ${DMRTYPE}/overlap_${FEATUREID}/${DMRTYPE}_${FEATUREID}_all_overlaps.with_coords.bed
# Drop coords for B, retain only ID
awk '{print $1"\t"$2"\t"$3"\t"$4";"$8"\t"$9}' \
		  ${DMRTYPE}/overlap_${FEATUREID}/${DMRTYPE}_${FEATUREID}_all_overlaps.with_coords.bed \
		  > ${DMRTYPE}/overlap_${FEATUREID}/${DMRTYPE}_${FEATUREID}_all_overlaps.bed

# Keep only unique DMRs that overlap
echo "#### Keeping only unique DMRs"
grep \
    -w -f \
    <( cat ${DMRTYPE}/overlap_${FEATUREID}/${DMRTYPE}_${FEATUREID}_all_overlaps.bed | awk '{print $4 }' | sed 's/[;].*//g' ) \
    ${DMRTYPE}/$AFILENAME.withids.pFilt${PVALTHR}.bed |\
    uniq \
    > ${DMRTYPE}/overlap_${FEATUREID}/unique_${DMRTYPE}_${FEATUREID}_overlaps.bed

# Get details on ATAC-overlapping DMRs
grep \
    -w -f \
    <( cat ${DMRTYPE}/overlap_${FEATUREID}/unique_${DMRTYPE}_${FEATUREID}_overlaps.bed | awk '{print $4}' ) \
    $AFILE \
    > ${DMRTYPE}/overlap_${FEATUREID}/${AFILENAME}_${FEATUREID}_overlap.tsv
cat \
   <( head -n 1 $AFILE ) \
   ${DMRTYPE}/overlap_${FEATUREID}/${AFILENAME}_${FEATUREID}_overlap.tsv \
   > /tmp/k && mv /tmp/k \
   ${DMRTYPE}/overlap_${FEATUREID}/${AFILENAME}_${FEATUREID}_overlap.tsv
grep \
    -w -f \
    <( cat ${DMRTYPE}/overlap_${FEATUREID}/${DMRTYPE}_${FEATUREID}_all_overlaps.with_coords.bed | awk '{print $8}' ) \
    $BFILENAME.withids.bed \
    > ${DMRTYPE}/overlap_${FEATUREID}/${BFILENAME}_${FEATUREID}overlap.bed

# Copy input to output folders
cp $AFILE ${DMRTYPE}/$AFILENAME.tsv
cp $BFILENAME.withids.bed ${DMRTYPE}/$BFILENAME.withids.bed
