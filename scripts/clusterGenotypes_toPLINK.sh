

# Input variables
MAIN_DIR=/mnt/f/SarahRice_mQTL/
SCRIPT_DIR=$MAIN_DIR/scripts/
DATA_DIR=$MAIN_DIR/rawdata/E181906_SarahRice_GSA_290323/vcf_files/
FASTA_FILE=$MAIN_DIR/genome/GRCh37/human_g1k_v37.fasta
IDAT_DIR=$MAIN_DIR/rawdata/E181906_SarahRice_GSA_290323/GSA_Idats/
BPM_FILE=$MAIN_DIR/rawdata/E181906_SarahRice_GSA_290323/microarray_manifest_data/GSA-24v3-0-A1-manifest-file-bpm/GSA-24v3-0_A1.bpm
EGT_FILE=$MAIN_DIR/rawdata/E181906_SarahRice_GSA_290323/microarray_manifest_data/GSA-24v3-0-A1-cluster-file/GSA-24v3-0_A1_ClusterFile.egt
GTC2VCF=$MAIN_DIR/utilities/genotyping/pre-plink/GTCtoVCF/gtc_to_vcf.py

# Output variables
VCF_DIR=$DATA_DIR/vcf_files/
PLINK_DIR=../qced_data/plink/

# Cluster genotypes
bash $SCRIPT_DIR/sidestep_genomestudio.grch37.sh \
     $IDAT_DIR \
     $FASTA_FILE \
     $BPM_FILE \
     $EGT_FILE \ 
     $GTC2VCF

# Replace non-numeric chromosome codes
bash $SCRIPT_DIR/replace_xymt.sh \
     *.vcf \
     $VCF_DIR

# Generate PLINK format
bash $SCRIPT_DIR/vcf2plink.sh \
     *numeric*.vcf \   
     $VCF_DIR 
     $PLINK_DIR