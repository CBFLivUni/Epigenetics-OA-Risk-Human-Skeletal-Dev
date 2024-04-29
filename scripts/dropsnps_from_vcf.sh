
# Variables
PATTERN=$1
INDIR=$2
OUTDIR=$INDIR

# Loop over vcf files
for vcf in ${INDIR}/*${PATTERN}*.vcf
do

  # Get vcf filenames
  vcfname=$( basename $vcf .vcf ); echo "##### $vcfname"

  # Get unwanted variants
  awk '$6 != "." && $1 !~ "#" {print $3}' $vcf > $OUTDIR/$vcfname.unwanted_variants.txt

  # Get and report number of dropped alleles
  n_snps=$( wc -l $OUTDIR/$vcfname.unwanted_variants.txt | awk '{print $1}' )
  echo "## There are " $n_snps " nuissance variants"

  # Remove
  grep -v -f \
      $OUTDIR/$vcfname.unwanted_variants.txt \
      $vcf > $OUTDIR/$vcfname.no_unwantedsnps.vcf

done
