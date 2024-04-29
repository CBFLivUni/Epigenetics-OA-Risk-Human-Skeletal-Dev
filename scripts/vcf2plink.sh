
# Define variables
FILE_PATTERN=$1
INDIR=$2
OUTDIR=$3

# Make output directory
mkdir -p $OUTDIR

# Loop over
for vcf in $INDIR*$FILE_PATTERN
do
  echo $vcf
  # Get name
  vcfname=$( basename $vcf .vcf )
  echo $vcfname

  # Remove variants with multiple alleles
  bcftools view \
               --max-alleles 2 \
               --exclude-types indels \
               $vcf > /tmp/${vcfname}.max2allele.vcf

  # VCF to Plink
  plink \
       --recode \
       --keep-allele-order \
       --vcf /tmp/${vcfname}.max2allele.vcf \
       --make-bed \
       --out ${OUTDIR}${vcfname}.plink

done
