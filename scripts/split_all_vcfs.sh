
# Define variables
PATTERN=$1
INDIR=$2
OUTDIR=$3
CURRENT_DIR=$(pwd)

# Make output directory
mkdir -p $OUTDIR

# Loop over files
for vcf in ${INDIR}${PATTERN}
do
  # Log current
  vcfname=$( basename $vcf .vcf.gz )
  echo "### ON :- $vcfname"

  # Index VCF
  cd $INDIR
  bcftools index -f $vcfname.vcf.gz
  cd $CURRENT_DIR

  # Split vcf
  bash splitvcf_by_chrm.sh $vcf $OUTDIR

done










