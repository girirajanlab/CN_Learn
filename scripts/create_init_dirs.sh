#!/bin/bash

echo "Job started on `hostname` at `date`"

source /data/CN-Learn/config.params
for directory in ${SOURCE_DIR} ${DATA_DIR} ${PRED_DIR} ${VALD_DIR} ${DATA_BPCOV_DIR};
do
if [ ! -d ${directory} ];
then
mkdir ${directory}
fi
done

echo "Job ended on `hostname` at `date`"
