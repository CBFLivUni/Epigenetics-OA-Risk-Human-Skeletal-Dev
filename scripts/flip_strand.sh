# Variables
PATTERN=$1
FASTA=$2
INDIR=$3
OUTDIR=$4

# Make output
mkdir -p $OUTDIR

# Loop over .bims
for bim in ${INDIR}${PATTERN}
do
  # Get bim name
  bimname=$( basename $bim .bim )
  echo "##### Working on :- $bimname"

  # Edit it .bim file to ensure that it matches up
  awk -F'\t' -vOFS='\t' '{ gsub("X", "23", $1) ; gsub("Y", "24", $1) ; gsub ("MT", "26", $1) ; print }' $bim > ${INDIR}${bimname}.edited.bim

  # Run SNP strand check
  eval "$(conda shell.bash hook)"
  conda activate python2 # snpflip is written in python 2 so need to switch
  echo "## Identifying -ve strand SNPs"
  snpflip \
         -b $bim \
         -f $FASTA \
         -o ${OUTDIR}${bimname}.snpflip

  # Flip strand
  eval "$(conda shell.bash hook)"
  conda activate plink # switch to plink environment
  echo "## Flipping -ve strand SNPs to +ve"
  plink \
       --bfile ${INDIR}${bimname} \
       --flip ${OUTDIR}${bimname}.snpflip.reverse \
       --make-bed \
       --out ${OUTDIR}${bimname}.strand_flip


done
