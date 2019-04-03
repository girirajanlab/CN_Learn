#!/bin/bash
#################################################################
# Author : Vijay Kumar
# Date   : 5/3/2018
# This is the master script for the CODEX pipeline that executes 
# all other steps via other shell scripts or R script.   
#
# Prereqs: None
#################################################################
echo "Job started on `hostname` at `date`"

source TBD/config.params

###############################################
# STEP 0: Declare variables and directory names 
###############################################
TARGET_FILE="targets_auto_no_chr.bed"
SAMPLE_FILE="list_of_bam_files_1"

#########################################################
# STEP 1: Iterate over all samples by chromosome number #
#########################################################
for chr_num in {1..22};do
Rscript ${SCRIPTS_EXT_DIR}codex.r ${chr_num} ${BAM_LIST_1} ${SOURCE_DIR}${SAMPLE_FILE} ${SOURCE_DIR}${TARGET_FILE} ${PROJ_NAME} ${DATA_EXT_CODEX_DIR}
done

#########################################
# STEP 2: Merge the output files together
#########################################
cp ${DATA_EXT_CODEX_DIR}codex_calls.txt ${DATA_EXT_CALLS_DIR}codex_calls.txt

echo "Job ended on `hostname` at `date`"
