#!/bin/sh
####################################################################
# Script : prechecks.sh                                            #
# Author : Vijay Kumar                                             #
# Date   : 4/5/2019                                                #
#                                                                  #
# This script performs multiple quality checks for the input files #
# required to run CN-Learn. Following are the list of tasks it     #
# performs,                                                        #
# 1) Update the path for the 'source' command in each bash script. #
# 2) Makes sure the paths for reference genome and bam files are   #
#    updated in the config.params file.                            #
# 3) Makes sure the directory with the BAM files is non-empty      #
# 4) Makes sure that all bam files have corresponding index files. #
# 5) Generates lists of samples, bam files and bams with full path #
# 6) Make sure the file with the list of exome capture probes is   #
#    present in the source directory.                              #
# 7) Remove 'CHR' from target probe list and generate a new file   #
#    witn only the autosomes.                                      #
# 8) Make sure the DOCKER_INDICATOR parameter is assigned valid    #
#    values and that updated to SW_DIR align with DOCKER_INDICATOR.#
#                        	                                   #
# (c) 2019 - Vijay Kumar	                                   #
# Licenced under the GNU General Public License 3.0.               #
####################################################################
echo "Task started on `hostname` at `date`"

#######################################################################
# STEP 1: Identify the absolute path of the current working directory #
#         and update the absolute path in all the downloaded scripts. #
#######################################################################
CURRENT_DIR=`pwd`'/'

sed -i "s|PROJ_DIR=TBD|PROJ_DIR=${CURRENT_DIR}|" ${CURRENT_DIR}config.params 
echo "STATUS: PROJ_DIR path in the config.params file has been updated successfully."

for script in `ls ${CURRENT_DIR}scripts/ | grep ".sh$"`; 
do 
sed -i "s|TBD/config.params|${CURRENT_DIR}config.params|" ${CURRENT_DIR}scripts/${script}
done
echo "STATUS: source path in all the bash scripts has been updated successfully."

#########################################################################
# STEP 2: Make sure the directory path with BAM files is updated in the #
#         config.params file and that the directory is not empty.       #
#########################################################################
source ${CURRENT_DIR}config.params
REF_GENOME_PATH=`cat ${CURRENT_DIR}config.params | grep -w "REF_GENOME=TBD" | wc -w`
BAM_FILE_PATH=`cat ${CURRENT_DIR}config.params | grep -w "BAM_FILE_DIR=TBD" | wc -w`

if [ ${REF_GENOME_PATH} -ne 0 ];
then
echo "ERROR: The REF_GENOME variable is not updated in the config.params file in "
echo "${CURRENT_DIR}. Please update this variable with "
echo "the absolute path of the directory hosting the reference genomes and rerun this script."
exit 1
else
if [ ! -f ${REF_GENOME} ];
then
echo ${REF_GENOME} 
echo "ERROR: The REF_GENOME_PATH variable does not seem to correspond to a valid file."
echo "Please make sure that the file path to the reference genome in the config.params file "
echo "is accurate and that the file exists and rerun this script."
exit 1
else
echo "STATUS: REF_GENOME file path has been updated."
fi
fi

if [ ${BAM_FILE_PATH} -ne 0 ]; 
then
echo "ERROR: The BAM_FILE_DIR variable is not updated in the config.params file in "
echo "${CURRENT_DIR}. Please update this variable with "
echo "the absolute path of the directory hosting the BAM files and rerun this script."
exit 1

else
if [ ! -d ${BAM_FILE_DIR} ];
then
echo "ERROR: The BAM_FILE_DIR variable does not seem to correspond to a working directory."
echo "Please make sure that the path to the directory with the BAM files provided in the "
echo "config.params file is accurate and that the directory exists and rerun this script."
exit 1
fi

echo "STATUS: BAM_FILE_DIR file path has been updated"
cd ${BAM_FILE_DIR}
bam_list=`ls | grep ".bam$"`
bam_count=`ls | grep ".bam$" | wc -l`
bai_count=`ls | grep "bam.bai$" | wc -l`

if [ ${bam_count} -eq 0 ];
then
echo "ERROR: The input directory with the list of BAM files does not seem to have "
echo "any BAM files. Please place the BAM files in the directory location provided"
echo "in the config.params file and rerun this script."
exit 1

elif [ ${bam_count} -gt 0 ];
then
echo "STATUS: Input BAM files are available for processing."

if [ ${bam_count} -ne ${bai_count} ];
then
echo "ERROR: The number of index files do not match the number of index files."
echo "Please make sure each bam file has an index file and rerun the script."
exit 1

else
echo "STATUS: Each bam file has a corresponding index file associated with it."
for bam_file in ${bam_list};
do
file_name_wc=`echo ${bam_file} | tr '.' ' ' | wc -w`
if [ ${file_name_wc} -ne 2 ];
then
echo "ERROR: The name of the file ${bam_file} does not follow the expected *.bam format."
echo "Please make sure that all the input bam files follow the format *.bam"
exit 1
fi
done

################################################################################
# STEP 3: Generate the list of samples, bam files and bam files with full path #
################################################################################
echo "STATUS: Creating the required input files with the list of sample names."
cd ${BAM_FILE_DIR}
ls | grep ".bam$" | sed -e 's/.bam//g'> ${SOURCE_DIR}'sample_list.txt'
ls | grep ".bam$" | awk -v path=${BAM_FILE_DIR} '{print path$0}' > ${BAM_FILE_LIST_W_PATH}
ls | grep ".bam$" > ${SOURCE_DIR}'bam_file_list.txt'

fi
fi
fi

#####################################################################################
# STEP 4: Generate the file with the list of samples required for the caller CANOES #
#####################################################################################
echo "STATUS: Creating the required input files with sample names for CANOES processing."
echo 'SAMPLE' > ${DATA_CANOES_DIR}sample_seq
num_samples=`cat ${BAM_FILE_LIST_W_PATH} | wc -l`
for samp in $(seq 1 ${num_samples});
do
echo "S${samp}" >> ${DATA_CANOES_DIR}sample_seq
done

echo 'SAMPLE_NAME' > ${DATA_CANOES_DIR}sample_names
cat ${SAMPLE_LIST} >> ${DATA_CANOES_DIR}sample_names
paste -d , ${DATA_CANOES_DIR}sample_seq ${DATA_CANOES_DIR}sample_names \
                             > ${DATA_CANOES_DIR}sample_map_canoes.csv

cp ${DATA_CANOES_DIR}sample_map_canoes.csv ${DATA_DIR}sample_map_canoes.csv


################################################################
# STEP 5: Create the file with the list of exome target probes
################################################################
if [ -f ${ORIG_PROBES} ];
then
chr_ind=`head -1 ${ORIG_PROBES} | cut -f1 | cut -c1-3 | tr a-z A-Z`
chr_num_ind=`head -1 ${ORIG_PROBES} | cut -c1`

re='^[0-9]+$'
if [ ${chr_ind} = "CHR" ];
then
cat ${ORIG_PROBES} | cut -f1-3 | tr a-z A-Z | \
         awk '{gsub("CHR", ""); if($1 != "X" && $1 != "Y") print $1"\t"$2"\t"$3}' \
                                          > ${SOURCE_DIR}'targets_auto_no_chr.bed'
elif [[ ${chr_ind} =~ $re ]];
then
cat ${ORIG_PROBES} | cut -f1-3 | tr a-z A-Z | \
         awk '{if($1 != "X" && $1 != "Y") print $1"\t"$2"\t"$3}' \
                                               > ${SOURCE_DIR}'targets_auto_no_chr.bed'
else
echo "ERROR: The exome capture prob file is not formatted correctly. Please make "
echo "sure that the input file is tab separated without headers."
exit 1 
fi
else
echo "ERROR: The input file with the list of exome capture probes is not present in the 'source'"
echo "directory. Please place the exome_capture_targets.bed file and rerun this script."
exit 1
fi

########################################################################
# STEP 6: Make sure the software paths align with the Docker indicator #
########################################################################
if [ ${DOCKER_INDICATOR} = 'Y' ] || [ ${DOCKER_INDICATOR} = 'y' ]; 
then
docker_install_indicator=`command -v docker | wc -l`
if [ ${docker_install_indicator} -eq 0 ]; 
then
echo "ERROR: DOCKER_INDICATOR is set to 'Y', but Docker is NOT currently installed."
echo "Please install Docker prior to executing any script that is part of CN-Learn." 
else
echo "STATUS: Docker is installed."
fi
elif [ ${DOCKER_INDICATOR} = 'N' ] || [ ${DOCKER_INDICATOR} = 'n' ]; 
then
docker_path=`cat ${CURRENT_DIR}config.params | grep '/opt/tools/CN_Learn/software/' | wc -l`
if [ ${docker_path} -ne 0 ]; 
then
echo "ERROR: The DOCKER_INDICATOR parameter is set to 'N', but the SW_DIR parameter in the "
echo "config.params file has not been updated yet. It still points to the path inside the "
echo "docker image. Please update this path and rerun this script."
echo ""
echo "NOTE: If not all softwares are installed in the same directory in the host, please review "
echo "and update each of the six paths listed under '4) Software Tools' in the config.params file."
exit 1
fi
else
echo "ERROR: The DOCKER_INDICATOR parameter can only take 'Y' or 'N' as its values. It appears "
echo "that neither values were supplied. Please make sure that the value is set correctly and "
echo "rerun this script."
exit 1
fi

echo "SUCCESS: All prechecks complete. Subsequent scripts can be now executed."

echo "Task ended on `hostname` at `date`"
