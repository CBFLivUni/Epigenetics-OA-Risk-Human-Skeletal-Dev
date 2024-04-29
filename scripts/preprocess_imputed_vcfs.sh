
# Variables
PATTERN=$1
FASTA=$2
INDIR=$3
REF=$4
OUTDIR=$INDIR

# Get list of chromosome files matching the pattern
CHRM_FILES=$( ls ${INDIR}${PATTERN} )
echo $CHRM_FILES

# Concatenate VCFs: Based off of https://github.com/huw-morris-lab/imputation
bcftools concat ${CHRM_FILES} -Ou |\
bcftools view -Ou -i 'R2>0.3' |\
bcftools norm -Ou -m -any |\
bcftools norm -Ou -f $FASTA |\
bcftools annotate -Oz -x ID -I +'%CHROM:%POS:%REF:%ALT' -o ${OUTDIR}imputed.$REF.vcf.gz
