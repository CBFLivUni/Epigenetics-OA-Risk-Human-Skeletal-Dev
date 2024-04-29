
# Define variables
INDIR=$1
OUTDIR=$2

# Make output directory
mkdir -p $OUTDIR

# Loop over
for plink in $INDIR*.bed
do

  # Get basename
  plinkname=$( basename $plink .bed )
  echo $plinkname

  # Calculate frequency
  plink \
       --nonfounders \
       --list-duplicate-vars suppress-first \
       --bfile ${INDIR}${plinkname} \
       --out ${OUTDIR}${plinkname}

  # Deduplicate
  plink \
       --bfile ${INDIR}${plinkname} \
       --exclude ${OUTDIR}${plinkname}.dupvar \
       --make-bed --out ${OUTDIR}${plinkname}.nodup

done
