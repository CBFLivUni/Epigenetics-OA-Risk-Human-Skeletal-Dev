# Variables
PATTERN=$1
INDIR=$2

# Loop over
for vcf in ${INDIR}${PATTERN}
do
  echo $vcf

  # Get vcf name
  vcfname=$( basename $vcf .vcf )
  echo $vcfname

  # Replace non-numeric chrm codes with respective numbers
  sed 's/X/23/g; s/Y/24/g; s/MT/26/g' \
                                     $vcf > $INDIR/$vcfname.xymt_numeric.vcf

done
