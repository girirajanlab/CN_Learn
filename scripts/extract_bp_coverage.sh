#!/bin/bash
################################################################################
# Author : Vijay Kumar                                                         #
# Date   : 7/25/2018                                                           #
# This script calculates basepair level coverage information for each sample.  #
# This information is used later to resolve breakpoints of concordant CNVs.    #
#                                                                              #
# IMPORTANT NOTE : It takes several hours for the genomeCoverageBed command to #
# finish for each sample. This script must be parallelized in order to get the #
# outputs for large cohorts in a reasonable amount of time.                    #
################################################################################
echo "Job started on `hostname` at `date`"

source /data/CN_Learn/config.params

for sample in `cat ${SAMPLE_LIST}`;
do

if [ ! -s ${DATA_BPCOV_DIR}${sample}.bpcov.bed ];
then
${BEDTOOLS_DIR}genomeCoverageBed -ibam ${BAM_FILE_DIR}${sample}.bam -bga \
                                     > ${DATA_BPCOV_DIR}${sample}.bpcov.bed
fi

done

echo "Job ended on `hostname` at `date`"

