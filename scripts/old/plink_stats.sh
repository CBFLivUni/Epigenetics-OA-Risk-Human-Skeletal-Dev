
# Define directories
INDIR=$(pwd)/$1
OUTDIR=$(pwd)/$2

# Make output directory
mkdir -p $OUTDIR
echo $OUTDIR
cd $OUTDIR

# Loop over
for plink in $INDIR*.bed
do
  # Get basename
  plinkname=$( basename $plink .bed )
  echo "##################################" $plinkname "##################################"

  # Calculate frequency
  plink \
       --nonfounders \
       --bfile ${INDIR}${plinkname} \
       --freq --out ${OUTDIR}${plinkname}.freq_out

  # Calculate frequency
  plink \
       --nonfounders \
       --bfile ${INDIR}${plinkname} \
       --missing --out ${OUTDIR}${plinkname}.missing

  # Calculate frequency
  plink \
       --nonfounders \
       --bfile ${INDIR}${plinkname} \
       --hardy --out ${OUTDIR}${plinkname}.hwe

done
