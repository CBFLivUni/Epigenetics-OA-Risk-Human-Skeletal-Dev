# Based on: https://github.com/freeseek/gtc2vcf#identifying-chip-type-for-idat-and-cel-files

# Currently doesn't work when running this as a bash script so just copy and paste for now.

# Input files
IDATDIR=$1 # /home/euancrna/projects/CBF/analysis/SarahRice_mQTL/rawdata/GSA_Idats/
GTCDIR=$( echo $IDATDIR | sed 's/GSA_Idats\//gtc_files\//g' )
VCFDIR=$( echo $IDATDIR | sed 's/GSA_Idats\//vcf_files\//g' )
FASTAFILE=$2
BPMFILE=$3 # /home/euancrna/projects/CBF/analysis/SarahRice_mQTL_2/rawdata/E181906_SarahRice_GSA_290323/manifest_data/GSA-24v3-0-A1-manifest-file-bpm/GSA-24v3-0_A1.bpm 
EGTFILE=$4 # /home/euancrna/projects/CBF/analysis/SarahRice_mQTL_2/rawdata/E181906_SarahRice_GSA_290323/manifest_data/GSA-24v3-0-A1-cluster-file/GSA-24v3-0_A1_ClusterFile.egt
GTC2VCFSCRIPT=$5 # /home/euancrna/projects/CBF/analysis/SarahRice_mQTL/utilities/GTCtoVCF/gtc_to_vcf.py

# Make output files
mkdir -p $GTCDIR; echo "##### OUTPUTTING GTC FILES TO: $GTCDIR"
mkdir -p $VCFDIR; echo "##### OUTPUTTING VCF FILES TO: $VCFDIR"

# CLR_ICU_VERSION_OVERRIDE="$(uconv -V | sed 's/.* //g')" LANG="en_US.UTF-8" $HOME/bin/iaap-cli/iaap-cli \
gencall \
	$BPMFILE \
 	$EGTFILE \
	$GTCDIR \
	--idat-folder $IDATDIR \
	--output-gtc \
       --gender-estimate-call-rate-threshold -0.1

# Get VCF files
eval "$(conda shell.bash hook)"
conda activate python2
python2 $GTC2VCFSCRIPT \
       --skip-indels \
       --gtc-paths $GTCDIR \
       --manifest-file $BPMFILE \
       --genome-fasta-file $FASTAFILE \
       --output-vcf-path $VCFDIR

eval "$(conda shell.bash hook)"
source activate plink
