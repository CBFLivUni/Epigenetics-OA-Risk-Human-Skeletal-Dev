# Derived from: https://www.biostars.org/p/479135/

# Variables
PATTERN=$1
INDIR=$2
OUTDIR=$3
CURRENTDIR=$(pwd)

# Remove previous merge file list
rm -rf $OUTDIR/file.list

# Loop over VCFs
for F in ${INDIR}${PATTERN}
do
      echo $F
     
      # Sort and convert vcf to bcf
      bcftools sort -O b -o ${F}.bcf $F

      # Index bcf
      cd $OUTDIR bcftools index ${F}.bcf

      # Add to list
      echo "${F}.bcf" >> file.list
done

# Merge listed files
bcftools merge  --file-list file.list  -O z -o merged.vcf.gz
