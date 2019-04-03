#!/bin/bash
########################################################################
# Author : Vijay Kumar                                                 #
# Date   : 5/3/2018                                                    #
# This is the master program for the CLAMMS pipeline that executes all # 
# other steps via other shell scripts or R/Python scripts.             #
########################################################################
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=clamms_pipeline
#SBATCH -o clamms_pipe.out
#SBATCH -e clamms_pipe.err
#SBATCH --ntasks=1
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=9G
#SBATCH --workdir /data/RF_Model_SimonsVIP/logs/extract/clamms/

echo "Job started on `hostname` at `date`"

source /data/RF_Model_SimonsVIP/config.params

####################################################
# STEP 0: Declare variables, directories and files #
####################################################
export INSERT_SIZE=200
export CLAMMS_DIR='/data/software/clamms/'

#################################################################################
# STEP 1: Convert bigwig to wig file to obtain mappability scores               # 
#################################################################################
cd ${MISC_DOWNLOADS}
if [ ! -f wgEncodeCrgMapabilityAlign75mer.bigWig ];
then
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign75mer.bigWig
fi

#/data/vijay/mappability/bigWigToWig wgEncodeCrgMapabilityAlign75mer.bigWig wgEncodeCrgMapabilityAlign75mer.wig
grep -v '^#' wgEncodeCrgMapabilityAlign75mer.wig | sed 's/^chr//g' | sort -k1,1 -k2,2n> ${DATA_EXT_CLAMMS_DIR}mappability.bed

#################################################################################
# STEP 2: Sort the target probe file and the mappability bed file               #
#################################################################################
sort -k1,1 -k2,2n ${TARGET_PROBES} > ${DATA_EXT_CLAMMS_DIR}sorted_targets.bed
sort -k1,1 -k2,2n ${DATA_EXT_CLAMMS_DIR}mappability.bed > ${DATA_EXT_CLAMMS_DIR}mappability_sorted.bed
sort -k1,1 -k2,2n ${CLAMMS_SW_DIR}data/clamms_special_regions.bed > ${DATA_EXT_CLAMMS_DIR}clamms_special_regions_sorted.bed

#accesses the special_regions.bed in CLAMMS
#chmod +x ${CLAMMS_DIR}/annotate_windows.sh
${CLAMMS_DIR}annotate_windows.sh ${DATA_EXT_CLAMMS_DIR}sorted_targets.bed ${REF_GENOME} ${DATA_EXT_CLAMMS_DIR}mappability_sorted.bed $INSERT_SIZE ${DATA_EXT_CLAMMS_DIR}clamms_special_regions_sorted.bed > ${DATA_EXT_CLAMMS_DIR}windows.bed

##################################################################
# STEP 3: Calculate actual and normalized coverage for each sample 
##################################################################
for bam_file in `cat ${BAM_LIST}`;
do

sample_name=`echo ${bam_file} | cut -f5 -d/`

echo "#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=clamms_cov_${sample_name}
#SBATCH -o clamms_cov_${sample_name}.out
#SBATCH -e clamms_cov_${sample_name}.err
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=32G
#SBATCH --workdir ${LOGS_EXT_CLAMMS_DIR}

source /data/RF_Model_SimonsVIP/config.params

echo "'Job started on `hostname` at `date`'"
${SAMTOOLS_DIR}samtools bedcov -Q 30 ${DATA_EXT_CLAMMS_DIR}windows.bed ${bam_file} | awk '"'{ printf ''"%s\t%d\t%d\t%.6g\n"'', $1, $2, $3, $NF/($3-$2); }'"' > ${CLAMMS_COV_DIR}${sample_name}.coverage.bed
${CLAMMS_SW_DIR}normalize_coverage ${CLAMMS_COV_DIR}${sample_name}.coverage.bed ${DATA_EXT_CLAMMS_DIR}windows.bed > ${CLAMMS_COV_DIR}${sample_name}.norm.cov.bed
echo "'Job ended on `hostname` at `date`'" " > ${SCRIPTS_EXT_DIR}dynamic_clamms.sh

#if [ ! -f ${CLAMMS_COV_DIR}${sample_name}.coverage.bed ];
#then
sbatch ${SCRIPTS_EXT_DIR}dynamic_clamms.sh
#echo ${sample_name}
#fi

done

#################################################################################
# Rename the normalized coverage files.This is done to make sure the participants
# within each family could be differentated
#################################################################################
cd ${CLAMMS_COV_DIR}
for norm_cov_file in `ls *.*.norm.cov.bed`;
do
mv ${norm_cov_file} `echo ${norm_cov_file} | sed 's/\./_/'`
done

#####################################################################################
# STEP 4: Combine normalized coverage into a single panel prior to fitting the model
#####################################################################################
rm ${DATA_EXT_CLAMMS_DIR}clamms_cnvs.bed
touch ${DATA_EXT_CLAMMS_DIR}clamms_cnvs.bed
for sample_group in 1 2;
do
sample_name=${SOURCE_DIR}sample_list_${sample_group}
ls ${CLAMMS_COV_DIR}*norm.cov.bed | grep -f ${sample_name} | while read FILE;
do
    echo "${FILE}"
    #grep "^Y" $FILE | awk '{ x += $4; n++; } END { if (x/(n+1) >= 0.1) print "M"; else print "F"; }'
done > ${DATA_EXT_CLAMMS_DIR}/ref.panel.files.${sample_group}.txt

${CLAMMS_SW_DIR}fit_models ${DATA_EXT_CLAMMS_DIR}ref.panel.files.${sample_group}.txt ${DATA_EXT_CLAMMS_DIR}windows.bed > ${DATA_EXT_CLAMMS_DIR}models.${sample_group}.bed

#####################################################################################
# STEP 5: Use the model file and normalized coverage file to make CNV predictions
#####################################################################################
cd ${CLAMMS_COV_DIR}
for norm_cov_file in `ls *norm.cov.bed | grep -f ${sample_name}`;
do 
${CLAMMS_SW_DIR}call_cnv ${CLAMMS_COV_DIR}${norm_cov_file} ${DATA_EXT_CLAMMS_DIR}models.${sample_group}.bed >> ${DATA_EXT_CLAMMS_DIR}clamms_cnvs.bed
done

done

sed 's/_/\./' ${DATA_EXT_CLAMMS_DIR}clamms_cnvs.bed > ${DATA_EXT_CLAMMS_DIR}clamms_cnvs_final.bed
cp ${DATA_EXT_CLAMMS_DIR}clamms_cnvs_final.bed ${DATA_EXT_CALLS_DIR}clamms_calls.txt

echo "Job ended on `hostname` at `date`"
