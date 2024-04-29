
# Define variables
FILE_PATTERN=$1
INDIR=$2
OUTDIR=$3

# Make output directory
mkdir -p $OUTDIR

# Loop over
for plink in ${INDIR}*${FILE_PATTERN}*.bed
do

  # Get basename
  plinkname=$( basename $plink .bed )
  echo $plinkname

  # Deduplicate
  plink \
       --check-sex \
       --bfile ${INDIR}${plinkname} \
       --make-bed --out ${OUTDIR}${plinkname}.checksex

done
