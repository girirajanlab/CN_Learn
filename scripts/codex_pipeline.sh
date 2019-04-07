#!/bin/bash
##################################################################
# Author : Vijay Kumar                                           #
# Date   : 4/5/2019                                              #
# This is the master script for the CODEX pipeline that executes #
# all other steps via other shell scripts or R script.           #
#                                                                #
# (c) 2018 - Vijay Kumar	                                 #
# Licenced under the GNU General Public License 3.0.             #
##################################################################
echo "Job started on `hostname` at `date`"

source /data/test_installation/CN_Learn/config.params

PROJ_NAME='cnlearn'
CODEX_OUTPUT_FILE='codex_calls.txt'

for chr_num in {1..20};do
docker run --rm -it -v ${PROJ_DIR}:${PROJ_DIR}:z girirajanlab/cnlearn \
Rscript ${RSCRIPTS_DIR}codex.r ${chr_num} ${BAM_FILE_LIST_W_PATH} \
            ${BAM_FILE_LIST} ${TARGET_PROBES} ${PROJ_NAME} ${DATA_CODEX_DIR}
done

#########################################################
# STEP 1: Iterate over all samples by chromosome number #
#########################################################
for chr_num in {21..21};do

echo "#!/bin/sh
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=codex_pipeline_${chr_num}
#SBATCH -o codex_pipe_${BATCH}_${chr_num}_out
#SBATCH -e codex_pipe_${BATCH}_${chr_num}_err
#SBATCH --ntasks=1
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=32G
#SBATCH --workdir /data/test_installation/CN_Learn/data/logs/
docker run --rm -it -v ${PROJ_DIR}:${PROJ_DIR}:z girirajanlab/cnlearn \
             Rscript ${RSCRIPTS_DIR}codex.r ${chr_num} ${BAM_FILE_LIST_W_PATH} ${BAM_FILE_LIST} ${TARGET_PROBES} ${PROJ_NAME} ${DATA_CODEX_DIR}" \
             > ${SCRIPTS_DIR}dynamic_codex.sh

sbatch ${SCRIPTS_DIR}/dynamic_codex.sh

done

########################################################
# Automated handling of file headers during data merge #
########################################################
if [ -f ${DATA_CODEX_DIR}${CODEX_OUTPUT_FILE} ];
then
rm ${DATA_CODEX_DIR}${CODEX_OUTPUT_FILE}
fi

first_file="yes"
for codex_file in ${DATA_CODEX_DIR}${PROJ_NAME}*;
do
        if [ ${first_file} = "yes" ]; 
        then
                echo "$codex_file copying"
                cp ${codex_file} ${DATA_CODEX_DIR}${CODEX_OUTPUT_FILE}
                first_file="no"
        else
                echo "${codex_file} omit first line"
                tail -n +2 -q ${codex_file} >> ${DATA_CODEX_DIR}${CODEX_OUTPUT_FILE}
        fi
done

#########################################
# STEP 2: Merge the output files together
#########################################
cat ${DATA_CODEX_DIR}codex_calls.txt | sed -e 's/\.bam\b//' > ${DATA_DIR}codex_calls.txt

echo "Job ended on `hostname` at `date`"
