
# Define directory variables
REF=$1
IN_DIR=$2
OUT_DIR=$3

# Make
mkdir $IN_DIR $OUT_DIR

# Loop over VCF files 
for vcf in $IN_DIR/*.vcf
do

  # Get name
  vcfname=$( basename $vcf .vcf )
  echo $vcfname

  # Run 
  perl vcf2maf.pl \
            --ref-fasta $REF \
            --input-vcf $vcf \
	    --output-maf $OUT_DIR/$vcfname.maf

done

