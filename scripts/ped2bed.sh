
# Define directories
PATTERN=$1
INDIR=$2
OUTDIR=$3

# Loop over
for file in ${INDIR}${PATTERN}
do

  # Get name
  filename=$( basename $file .ped )
  echo $filename

  plink \
       --file ${INDIR}${filename} \
       --keep-allele-order \
       --make-bed \
       --out ${OUTDIR}${filename}
done
