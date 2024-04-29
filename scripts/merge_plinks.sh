# Define variables
RUNDIR=$( pwd )

PATTERN=$1
OUTMERGE=$2
# SPLITX_NEWVAL=$2
INDIR=${RUNDIR}/${3}
OUTDIR=${RUNDIR}/${4}
OUTMERGE=${OUTDIR}${OUTMERGE}

# Get universal pattern to match all files of interest irregardless of file suffix
UNIVERSAL_PATTERN=$(echo $PATTERN | sed 's/[.]bed//g'); echo $UNIVERSAL_PATTERN





# Make output directory
mkdir -p $OUTDIR
mkdir -p ${INDIR}../strand_flip/

# Check number of files
ls -lh ${INDIR}${PATTERN} | wc -l





# Try and merge
echo " ##### Performing initial merge (to get nuissance SNPS with >3 alleles) ....."

# Remove old files
cd $INDIR
rm *${OUTMERGE}* *tmp* *temporary*




# Most of the operations in this loop are derived from: https://www.biostars.org/p/101191/
# Attempt to merge then perform strand flip
for plink in ${INDIR}${PATTERN}
do
  # Get sample ID
  plinkname=$( basename $plink .bed )
  echo "################### WORKING ON : $plinkname ###################"

  # Make new merge file without the sample of interest
  ls -lh ${INDIR}${UNIVERSAL_PATTERN}*.bim | awk '{print $9}' > bims.txt
  ls -lh ${INDIR}${UNIVERSAL_PATTERN}*.fam | awk '{print $9}' > fams.txt
  ls -lh ${INDIR}${UNIVERSAL_PATTERN}*.bed | awk '{print $9}' > beds.txt
  paste beds.txt bims.txt fams.txt | grep -v $plinkname > ${OUTDIR}files_to_merge.txt
  wc -l ${OUTDIR}files_to_merge.txt

  # Attempt merge (to get .missnp file)
  echo " ##### ATTEMPTING FIRST MERGE"
  plink \
       --bfile ${INDIR}${plinkname} \
       --keep-allele-order \
       --merge-list ${OUTDIR}files_to_merge.txt \
       --out /tmp/${plinkname}.1st_merge
  cd $RUNDIR

  # If any missnps
  if [[ -f ${OUTDIR}${plinkname}.1st_merge.missnp ]];
  then

    # Get number of missnps
    n_missnp=$( wc -l ${OUTDIR}${plinkname}.1st_merge.missnp )
    echo $n_misssnp
    echo "   ## There are " ${n_missnp:1:0} " awkward SNPs"

    # Attempt to flip SNPs
    echo " ##### ATTEMPTING TO FLIP NUISSANCE SNPS"
    plink \
       --bfile ${INDIR}${plinkname} \
       --keep-allele-order \
       --flip /tmp/${plinkname}.1st_merge.missnp \
       --make-bed \
       --out ${INDIR}../strand_flip/${plinkname}.strand_flip


    echo ${INDIR}../strand_flip/${plinkname}.strand_flip

    # Attempt merge to identify current nuissance SNPs
    echo " ##### MERGE TO IDENTIFY REMAINING NUISSANCE SNPS (PROBABLY > bi allelic)"
    plink \
       --bfile ${INDIR}../strand_flip/${plinkname}.strand_flip \
       --keep-allele-order \
       --merge-list ${OUTDIR}files_to_merge.txt \
       --make-bed \
       --out ${OUTDIR}${plinkname}.final_merge
    cd $RUNDIR

    # Exclude nuissance SNPs
    # merge-merge is there because PLINK's doing some strange wizardry
    echo " ##### REMOVING NUISSANCE SNPS NOT FIXED BY STRAND FLIPPING"
    plink \
       --bfile ${INDIR}../strand_flip/${plinkname}.strand_flip \
       --keep-allele-order \
       --exclude ${OUTDIR}${plinkname}.final_merge-merge.missnp \
       --make-bed \
       --out ${OUTDIR}${plinkname}.no_multi_snps
  else

    #
    echo "   ## There are no awkward SNPs"

    # Else if no missing SNPs, just rename temp file to the output
    mv /tmp/${plinkname}.1st_merge.bed ${OUTDIR}${plinkname}.no_multi_snps.bed
    mv /tmp/${plinkname}.1st_merge.bim ${OUTDIR}${plinkname}.no_multi_snps.bim
    mv /tmp/${plinkname}.1st_merge.fam ${OUTDIR}${plinkname}.no_multi_snps.fam

  fi
done





# Make new merge file
echo " ##### PERFORMING FINAL MERGE.........."
ls -lh ${OUTDIR}*${UNIVERSAL_PATTERN}*no_multi_snps*.bim | awk '{print $9}' > bims.txt
ls -lh ${OUTDIR}*${UNIVERSAL_PATTERN}*no_multi_snps*.fam | awk '{print $9}' > fams.txt
ls -lh ${OUTDIR}*${UNIVERSAL_PATTERN}*no_multi_snps*.bed | awk '{print $9}' > beds.txt
paste beds.txt bims.txt fams.txt > ${OUTDIR}files_to_merge.txt

# Merge
plink --keep-allele-order --merge-list ${OUTDIR}files_to_merge.txt --out $OUTMERGE
cd $RUNDIR










