# Based on: https://github.com/freeseek/gtc2vcf#identifying-chip-type-for-idat-and-cel-files

# Currently doesn't work when running this as a bash script so just copy and paste for now.

CLR_ICU_VERSION_OVERRIDE="$(uconv -V | sed 's/.* //g')" LANG="en_US.UTF-8" $HOME/bin/iaap-cli/iaap-cli \
gencall \
	 /home/euancrna/projects/CBF/analysis/SarahRice_mQTL/rawdata/E181906_SarahRice_GSA_290323/manifest_data/GSA-24v3-0-A2-manifest-file-bpm/GSA-24v3-0_A2.bpm \
	 /home/euancrna/projects/CBF/analysis/SarahRice_mQTL/rawdata/E181906_SarahRice_GSA_290323/manifest_data/GSA-24v3-0-A1-cluster-file/GSA-24v3-0_A1_ClusterFile.egt \
	 /home/euancrna/projects/CBF/analysis/SarahRice_mQTL/rawdata/gtc_files/ \
	 --idat-folder /home/euancrna/projects/CBF/analysis/SarahRice_mQTL/rawdata/GSA_Idats/ \
	 --output-gtc \
       --gender-estimate-call-rate-threshold -0.1


# This doesn't seem to work 
# bcftools +gtc2vcf \
#  	    --no-version -Ou \
#	    --egt /home/euancrna/projects/CBF/analysis/SarahRice_mQTL/rawdata/E181906_SarahRice_GSA_290323/manifest_data/GSA-24v3-0-A1-cluster-file/GSA-24v3-0_A1_ClusterFile.egt \
#	    --csv /home/euancrna/projects/CBF/analysis/SarahRice_mQTL/rawdata/E181906_SarahRice_GSA_290323/manifest_data/GSA-24v3-0-A2-manifest-file-csv/GSA-24v3-0_A2.csv \
#         --gtcs /home/euancrna/projects/CBF/analysis/SarahRice_mQTL/rawdata/gtc_files \
#         --fasta-ref ~/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
#         --extra /home/euancrna/projects/CBF/analysis/SarahRice_mQTL/rawdata/vcf_files/srice_genotyping.gtc2vcf_details.tsv

# Instead use this (which is what the above is a wrapper for)
python /home/euancrna/projects/CBF/analysis/SarahRice_mQTL/utilities/GTCtoVCF/gtc_to_vcf.py \
       --skip-indels \
       --gtc-paths /home/euancrna/projects/CBF/analysis/SarahRice_mQTL/rawdata/gtc_files/ \
       --manifest-file /home/euancrna/projects/CBF/analysis/SarahRice_mQTL/rawdata/E181906_SarahRice_GSA_290323/manifest_data/GSA-24v3-0-A2-manifest-file-bpm/GSA-24v3-0_A2.bpm \
       --genome-fasta-file ../utilities/genome/GRCh38/38 \
       --output-vcf-path /home/euancrna/projects/CBF/analysis/SarahRice_mQTL/rawdata/vcf_files/
