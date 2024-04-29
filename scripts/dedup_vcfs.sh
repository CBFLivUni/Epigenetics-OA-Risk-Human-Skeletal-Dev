# Based on: https://www.biostars.org/p/335605/

# Set variables
REFNAME=$1
INDIR=$2
OUTDIR=$3

# Loop
for vcf in ${INDIR}*.vcf
do
    # Get name
    vcfname=$( basename $vcf .vcf )
    echo $vcfname 

    # Run BCFtools
    bcftools norm -m-any --check-ref w -f $REFNAME $vcf | \
      bcftools annotate -I +'%CHROM:%POS:%REF:%ALT' | \
        bcftools norm -Ob --rm-dup both \
          > ${OUTDIR}${vcfname}.dedup.bcf

    # Index BCF
    bcftools index ${OUTDIR}${vcfname}.dedup.bcf
done
