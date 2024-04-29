
# Define variables
PATTERN=$1
INDIR=$2
METAFILE=$3

# Loop over files
for plink in ${INDIR}${PATTERN}
do
  # Get plink filename id
  plinkname=$( basename $plink .bed )
  echo "### On $plinkname"

  # Run R script to convert file ids
  Rscript fix_fam_ids.R \
	  ${INDIR}${plinkname}.fam \
          $METAFILE
done
