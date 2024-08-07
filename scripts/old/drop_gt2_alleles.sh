
# Define directories
MAXALLELE=$1
MINALLELE=$2
INDIR=$3
OUTDIR=$4



# Loop over vcf files
for vcf in $INDIR*.vcf
do
  # Get sample id
  vcfname=$( basename $vcf .vcf )
  echo $vcfname  

  # Run 
  bcftools view \
	   -Oz \
	   --max-alleles $MAXALLELE \
	   --min-alleles $MINALLELE $vcf \
	   -o $OUTDIR/$vcfname.2allele.vcf

done
