
# Define variables
MAF_VAL=$1
MIND_VAL=$2
SNP_VAL=$3
HWE_VAL=$4
INDIR=$5
OUTDIR=$6

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
       --maf $MAF_VAL \
       --mind $MIND_VAL \
       --geno $SNP_VAL \
       --hwe $HWE_VAL \
       --bfile ${INDIR}${plinkname} \
       --make-bed --out ${OUTDIR}${plinkname}.q_filt
done
