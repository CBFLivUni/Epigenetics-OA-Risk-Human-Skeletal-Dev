
# Variables
IN_VCF=$1
OUTDIR=$2

# Loop over chromosomes
for chrm in {1..22} X Y MT
do
  # Get VCF name
  vcfname=$( basename $IN_VCF .vcf.gz)

  # Split by chromosomes
  bcftools view $IN_VCF -Oz --regions ${chrm} -o ${OUTDIR}/${vcfname}.chr${chrm}.vcf.gz

  # Unzip
  bcftools view ${OUTDIR}/${vcfname}.chr${chrm}.vcf.gz -Ob -o ${OUTDIR}/${vcfname}.chr${chrm}.vcf
done
