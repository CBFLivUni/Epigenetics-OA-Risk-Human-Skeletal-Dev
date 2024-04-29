# Define variables
TOKEN=$( cat $1 )
PATHBASENAME=$2

# Check token
echo "##### TOKEN IS :- $TOKEN"

# Run with API
curl https://imputationserver.sph.umich.edu/api/v2/jobs/submit/minimac4 \
      -H "X-Auth-Token: $TOKEN" \
      -F "files=@${PATHBASENAME}.chr1.vcf.gz" \
      -F "files=@${PATHBASENAME}.chr10.vcf.gz" \
      -F "files=@${PATHBASENAME}.chr11.vcf.gz" \
      -F "files=@${PATHBASENAME}.chr12.vcf.gz" \
      -F "files=@${PATHBASENAME}.chr13.vcf.gz" \
      -F "files=@${PATHBASENAME}.chr14.vcf.gz" \
      -F "files=@${PATHBASENAME}.chr15.vcf.gz" \
      -F "files=@${PATHBASENAME}.chr16.vcf.gz" \
      -F "files=@${PATHBASENAME}.chr17.vcf.gz" \
      -F "files=@${PATHBASENAME}.chr18.vcf.gz" \
      -F "files=@${PATHBASENAME}.chr19.vcf.gz" \
      -F "files=@${PATHBASENAME}.chr2.vcf.gz" \
      -F "files=@${PATHBASENAME}.chr20.vcf.gz" \
      -F "files=@${PATHBASENAME}.chr21.vcf.gz" \
      -F "files=@${PATHBASENAME}.chr22.vcf.gz" \
      -F "files=@${PATHBASENAME}.chr3.vcf.gz" \
      -F "files=@${PATHBASENAME}.chr4.vcf.gz" \
      -F "files=@${PATHBASENAME}.chr5.vcf.gz" \
      -F "files=@${PATHBASENAME}.chr6.vcf.gz" \
      -F "files=@${PATHBASENAME}.chr7.vcf.gz" \
      -F "files=@${PATHBASENAME}.chr8.vcf.gz" \
      -F "files=@${PATHBASENAME}.chr9.vcf.gz" \
      -F "refpanel=hrc-r1.1" \
      -F "population=eur" \
      -F "build=hg19"

