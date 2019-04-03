#!/bin/bash
#########################################################################################
# Author : Vijay Kumar                                                                  #
# Date   : 5/3/2018                                                                     #
# This is the master program for the XHMM pipeline that executes all other steps via    #
# other shell scripts or R/Python scripts.                                              #
#########################################################################################
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=XHMM_pipe
#SBATCH -o XHMM_pipe_out
#SBATCH -e XHMM_pipe_err
#SBATCH --ntasks=1
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --workdir /data/RF_Model_SimonsVIP/logs/extract/xhmm/

echo "Job started on `hostname` at `date`"

source /data/RF_Model_SimonsVIP/config.params

##############################################################################################
# STEP 0: Declare variables, directory locations and other required parameters for the script.
##############################################################################################
for bam_file in `cat ${BAM_LIST}`;
do

echo ${bam_file}
echo ${bam_file} | cut -f5 -d/
sample_name=`echo ${bam_file} | cut -f5 -d/`

#############################################################################
# STEP 1: Generate separate scripts to calculate read counts for each sample
#############################################################################
echo "#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=readcount_${sample_name}.sh
#SBATCH -o readcount_${sample_name}_out
#SBATCH -e readcount_${sample_name}_err
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --workdir ${LOGS_EXT_XHMM_DIR}

java -Xmx3072m -Djava.io.tmpdir=${LOGS_EXT_XHMM_DIR}xhmm_temp -jar ${GATK_SW_DIR}GenomeAnalysisTK.jar -T DepthOfCoverage -I ${bam_file} -L ${TARGET_PROBES} -R ${REF_GENOME} -dt BY_SAMPLE -dcov 5000 -l INFO --omitDepthOutputAtEachBase --omitLocusTable --minBaseQuality 0 --minMappingQuality 20 --start 1 --stop 5000 --nBins 200 --includeRefNSites --countType COUNT_FRAGMENTS -o ${XHMM_RD_DIR}${sample_name}_GATK_OUT" > ${SCRIPTS_EXT_DIR}dynamic_xhmm_rd.sh

if [ ! -f ${XHMM_RD_DIR}${sample_name}_GATK_OUT.sample_summary ];
then
sbatch ${SCRIPTS_EXT_DIR}dynamic_xhmm_rd.sh
fi

done

#rm ${DATA_EXT_XHMM_DIR}xcnv
touch ${DATA_EXT_XHMM_DIR}xcnv

for pattern in 14 S;
do
##################################################################
# STEP 2: Combine all GATK depth-of-coverage outputs into one file
##################################################################
# --GATKdepthsList is a list of all *.list.DATA.sample_interval_summary files to be combined
cd ${XHMM_RD_DIR}
ls | grep GATK_OUT.sample_interval_summary | egrep ^${pattern} > ${DATA_EXT_XHMM_DIR}read_depth_file_list
${XHMM_SW_DIR}xhmm --mergeGATKdepths -o ${DATA_EXT_XHMM_DIR}combined_RD.txt --GATKdepthsList ${DATA_EXT_XHMM_DIR}read_depth_file_list

##############################
# STEP 3: Calculate GC Content
##############################
java -Xmx3072m -jar ${GATK_SW_DIR}GenomeAnalysisTK.jar -T GCContentByInterval -L ${TARGET_PROBES} -R ${REF_GENOME} -o ${DATA_EXT_XHMM_DIR}GC.txt
cat ${DATA_EXT_XHMM_DIR}GC.txt | awk '{if ($2 < 0.1 || $2 > 0.9) print $1}' > ${DATA_EXT_XHMM_DIR}extreme_gc_targets.txt

#########################################################################
# STEP 4: Run PLINK/Seq to filter GATK targets with repeat-maked bases
#########################################################################
#Use custom Python script (would AWK work??) to convert from BED to REG; may need to convert back to Unix-style new lines
#bed_to_reg.py and seqdb.hg19 are stored in PLINK-Seq program directory
python ${PLINK_DIR}bed_to_reg.py ${TARGET_PROBES}
mv ${TARGET_PROBES}.reg ${DATA_EXT_XHMM_DIR}${TARGET_PROBE_FILE}.reg
sed 's///g' ${DATA_EXT_XHMM_DIR}${TARGET_PROBE_FILE}.reg > ${DATA_EXT_XHMM_DIR}targets_no_chr.reg
${PLINK_DIR}pseq . loc-load --locdb ${DATA_EXT_XHMM_DIR}targets.LOCDB --file ${DATA_EXT_XHMM_DIR}targets_no_chr.reg --group targets --out ${DATA_EXT_XHMM_DIR}targets.LOCDB.loc-load

#####################################
# STEP 5: Identify complexity of loci
#####################################
${PLINK_DIR}pseq . loc-stats --locdb ${DATA_EXT_XHMM_DIR}targets.LOCDB --group targets --seqdb ${PLINK_DIR}seqdb.hg19 | awk '{if (NR > 1) print $_}' | sort -k1 -g | awk '{print $10}' | paste ${DATA_EXT_XHMM_DIR}targets_no_chr.reg - > ${DATA_EXT_XHMM_DIR}locus_complexity.txt
cat ${DATA_EXT_XHMM_DIR}locus_complexity.txt | awk '{if ($4 > 0.25) print $1 ":" $2 "-" $3}' > ${DATA_EXT_XHMM_DIR}low_complexity_targets.txt

#############################################################################################################################
# STEP 6: Filter samples by GC content/low complexity, and then find mean-center of targets. Use default parameters for now
#############################################################################################################################
${XHMM_SW_DIR}xhmm --matrix -r ${DATA_EXT_XHMM_DIR}combined_RD.txt --centerData --centerType target -o ${DATA_EXT_XHMM_DIR}centered_RD.txt --outputExcludedTargets ${DATA_EXT_XHMM_DIR}filt_targets.txt --outputExcludedSamples ${DATA_EXT_XHMM_DIR}filt_samples.txt --excludeTargets ${DATA_EXT_XHMM_DIR}extreme_gc_targets.txt --excludeTargets ${DATA_EXT_XHMM_DIR}low_complexity_targets.txt --minTargetSize 10 --maxTargetSize 10000 --minMeanTargetRD 10 --maxMeanTargetRD 500 --minMeanSampleRD 25 --maxMeanSampleRD 200 --maxSdSampleRD 150

###########################################################################
# STEP 7: Run PCA on the mean-centered data & normalize based on PCA values
###########################################################################
${XHMM_SW_DIR}xhmm --PCA -r ${DATA_EXT_XHMM_DIR}centered_RD.txt --PCAfiles ${DATA_EXT_XHMM_DIR}PCA_output.txt

#Normalize data based on PCA values
${XHMM_SW_DIR}xhmm --normalize -r ${DATA_EXT_XHMM_DIR}centered_RD.txt --PCAfiles ${DATA_EXT_XHMM_DIR}PCA_output.txt --normalizeOutput ${DATA_EXT_XHMM_DIR}PCA_normalized.txt --PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7

############################################################################
# STEP 8: Filter and calculate z-scores for targets from PCA-normalized data
############################################################################
${XHMM_SW_DIR}xhmm --matrix -r ${DATA_EXT_XHMM_DIR}PCA_normalized.txt --centerData --centerType sample --zScoreData -o ${DATA_EXT_XHMM_DIR}PCA_zscores.txt --outputExcludedTargets ${DATA_EXT_XHMM_DIR}PCA_z_filt_targets.txt --outputExcludedSamples ${DATA_EXT_XHMM_DIR}PCA_z_filt_samp.txt --maxSdTargetRD 30

#Filter original read-depth data against normalized data
${XHMM_SW_DIR}xhmm --matrix -r ${DATA_EXT_XHMM_DIR}combined_RD.txt --excludeTargets ${DATA_EXT_XHMM_DIR}filt_targets.txt --excludeTargets ${DATA_EXT_XHMM_DIR}PCA_z_filt_targets.txt --excludeSamples ${DATA_EXT_XHMM_DIR}PCA_z_filt_samp.txt -o ${DATA_EXT_XHMM_DIR}same_filtered_RD.txt

######################################
# STEP 9: Call CNVs on normalized data
######################################
#Use default parameters file (params.txt) in XHMM software file; may want to experiment with these settings
${XHMM_SW_DIR}xhmm --discover -p ${XHMM_SW_DIR}params.txt -r ${DATA_EXT_XHMM_DIR}PCA_zscores.txt -R ${DATA_EXT_XHMM_DIR}same_filtered_RD.txt -c ${DATA_EXT_XHMM_DIR}xcnv_${pattern} -a ${DATA_EXT_XHMM_DIR}aux_xcnv -s ${DATA_EXT_XHMM_DIR}

done

rm ${DATA_EXT_XHMM_DIR}xcnv
touch ${DATA_EXT_XHMM_DIR}xcnv

cat ${DATA_EXT_XHMM_DIR}xcnv_14 >> ${DATA_EXT_XHMM_DIR}xcnv
cat ${DATA_EXT_XHMM_DIR}xcnv_S | tail -n+2 >> ${DATA_EXT_XHMM_DIR}xcnv

cp ${DATA_EXT_XHMM_DIR}'xcnv' ${DATA_EXT_CALLS_DIR}'xhmm_calls.txt'

echo "Job ended on `hostname` at `date`"
