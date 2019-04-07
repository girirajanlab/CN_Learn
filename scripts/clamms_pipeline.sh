#!/bin/bash
########################################################################
# Author : Vijay Kumar                                                 #
# Date   : 5/3/2018                                                    #
# This is the master program for the CLAMMS pipeline that executes all # 
# other steps via other shell scripts or R/Python scripts.             #
########################################################################
echo "Job started on `hostname` at `date`"

source /data/test_installation/CN_Learn/config.params

####################################################
# STEP 0: Declare variables, directories and files #
####################################################
export INSERT_SIZE=200
#export CLAMMS_DIR='/data/software/clamms/'

#################################################################################
# STEP 1: Convert bigwig to wig file to obtain mappability scores               # 
#################################################################################
cd ${SOURCE_DIR}
if [ ! -f wgEncodeCrgMapabilityAlign75mer.bigWig ];
then
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign75mer.bigWig
fi

docker run --rm -ti -v ${PROJ_DIR}:${PROJ_DIR} girirajanlab/cnlearn \
            ${SW_DIR}bigWigToWig ${SOURCE_DIR}wgEncodeCrgMapabilityAlign75mer.bigWig \
/data/vijay/mappability/bigWigToWig ${SOURCE_DIR}wgEncodeCrgMapabilityAlign75mer.bigWig \
                                 ${SOURCE_DIR}wgEncodeCrgMapabilityAlign75mer.wig

grep -v '^#' wgEncodeCrgMapabilityAlign75mer.wig | sed 's/^chr//g' | sort -k1,1 -k2,2n \
                                   > ${DATA_CLAMMS_DIR}mappability.bed

#################################################################################
# STEP 2: Sort the target probe file and the mappability bed file               #
#################################################################################
sort -k1,1 -k2,2n ${TARGET_PROBES} > ${DATA_CLAMMS_DIR}sorted_targets.bed
sort -k1,1 -k2,2n ${DATA_CLAMMS_DIR}mappability.bed > ${DATA_CLAMMS_DIR}mappability_sorted.bed

docker run --rm -ti -v ${PROJ_DIR}:${PROJ_DIR} girirajanlab/cnlearn \
sort -k1,1 -k2,2n ${CLAMMS_DIR}data/clamms_special_regions.bed \
             > ${DATA_CLAMMS_DIR}clamms_special_regions_sorted.bed


#accesses the special_regions.bed in CLAMMS
docker run --rm -ti -v ${PROJ_DIR}:${PROJ_DIR} girirajanlab/cnlearn \
${CLAMMS_DIR}annotate_windows.sh ${DATA_CLAMMS_DIR}sorted_targets.bed ${REF_GENOME} \
                                     ${DATA_CLAMMS_DIR}mappability_sorted.bed $INSERT_SIZE \
                                     ${DATA_CLAMMS_DIR}clamms_special_regions_sorted.bed \
                                                     > ${DATA_CLAMMS_DIR}windows.bed

##################################################################
# STEP 3: Calculate actual and normalized coverage for each sample 
##################################################################
for bam_file in `cat ${BAM_FILE_LIST_W_PATH} | head -1`;
do

echo ${bam_file}
sample_name=`echo ${bam_file} | rev | cut -f1 -d/ | rev`

docker run --rm -v ${PROJ_DIR}:${PROJ_DIR} girirajanlab/cnlearn \
${SAMTOOLS_DIR}samtools bedcov -Q 30 ${DATA_CLAMMS_DIR}windows.bed ${bam_file} | \
                   awk '{ printf "%s\t%d\t%d\t%.6g\n", $1, $2, $3, $NF/($3-$2); }' \
                                     > ${DATA_CLAMMS_DIR}${sample_name}.coverage.bed

docker run --rm -t -v ${PROJ_DIR}:${PROJ_DIR} girirajanlab/cnlearn \
${CLAMMS_DIR}normalize_coverage ${DATA_CLAMMS_DIR}${sample_name}.coverage.bed \
                                ${DATA_CLAMMS_DIR}windows.bed \
                                > ${DATA_CLAMMS_DIR}${sample_name}.norm.cov.bed
done

#####################################################################################
# STEP 4: Combine normalized coverage into a single panel prior to fitting the model
#####################################################################################
if [ -f ${DATA_CLAMMS_DIR}clamms_cnvs.bed ];
then
rm ${DATA_CLAMMS_DIR}clamms_cnvs.bed
fi
touch ${DATA_CLAMMS_DIR}clamms_cnvs.bed

sample_name=${SAMPLE_LIST}
ls ${DATA_CLAMMS_DIR}*norm.cov.bed | grep -f ${sample_name} | while read FILE;
do
    echo "${FILE}"
    #grep "^Y" $FILE | awk '{ x += $4; n++; } END { if (x/(n+1) >= 0.1) print "M"; else print "F"; }'
done > ${DATA_CLAMMS_DIR}/ref.panel.files.txt

docker run --rm -v ${PROJ_DIR}:${PROJ_DIR} girirajanlab/cnlearn \
${CLAMMS_DIR}fit_models ${DATA_CLAMMS_DIR}ref.panel.files.txt ${DATA_CLAMMS_DIR}windows.bed \
                                                                > ${DATA_CLAMMS_DIR}models.bed

#####################################################################################
# STEP 5: Use the model file and normalized coverage file to make CNV predictions
#####################################################################################
cd ${DATA_CLAMMS_DIR}
sample_name=${SAMPLE_LIST}
for norm_cov_file in `ls *norm.cov.bed | grep -f ${sample_name}`;
do 
docker run --rm -v ${PROJ_DIR}:${PROJ_DIR} girirajanlab/cnlearn \
${CLAMMS_DIR}call_cnv ${DATA_CLAMMS_DIR}${norm_cov_file} ${DATA_CLAMMS_DIR}models.bed \
                                                       >> ${DATA_CLAMMS_DIR}clamms_cnvs.bed
done

cp ${DATA_CLAMMS_DIR}clamms_cnvs.bed ${DATA_DIR}clamms_calls.txt

echo "Job ended on `hostname` at `date`"
