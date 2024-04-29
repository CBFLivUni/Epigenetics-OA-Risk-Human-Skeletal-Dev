# Define variables
TOKEN=$( cat $1 )
PATTERN=$2
INDIR=$3

# Check token
echo "##### TOKEN IS :- $TOKEN"

# Loop over files
for vcf in ${INDIR}${PATTERN}
do
  # Get just path
  vcf=$(echo "@${vcf}" | sed '' )

  # Log current
  echo "### ON :- $vcf"

  # Run with API
  curl https://imputationserver.sph.umich.edu/api/v2/jobs/submit/minimac4 \
      -H "X-Auth-Token: $TOKEN" \
      -F "files=$vcf" \
      -F "refpanel=1000g-phase-3-v5" \
      -F "population=eur" \
      -F "build=hg19" \
      -F "r2Filter=0.01"

  # Sleep to let jobs finish
  sleep 600

done
