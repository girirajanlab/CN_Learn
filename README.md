# CN-Learn 
CN-Learn is a framework to integrate Copy Number Variant (CNV) predictions made by multiple algorithms using exome sequencing datasets. Unlike traditional integration methods that depend on a single measure of concordance, CN-Learn leverages the extent of concordance among the calls made by multiple methods, in addition to genomic contexts such as GC content, mappability and local read depth fluctuations. Using a wide range of predictors extracted for CNVs from a “gold-standard” CNV call set, CN-Learn builds a ‘Random Forest’ classifier with hundreds of decision trees to estimate the probability of each CNV in a fresh set of samples being true. This strategy provides CN-Learn the ability to segregate true CNV calls from false positives in an informed fashion with extremely high precision.

While CN-Learn has been shown to perform best when built as a ‘Random Forest’ classifier, it can also be built as a ‘Logistic Regression’ or ‘Support Vector Machine’ classifier. CN-Learn can also be seamlessly extended to include newer set of CNV calling algorithms in the future by simply changing a single parameter supplied to CN-Learn. 

## Software Requirements
Given the number of softwares and tools required to run the four individual CNV calling algorithms prior to running CN-Learn, every software/tool required to run CN-Learn end-to-end has been packaged into a Docker image. Following are some of the preinstalled software tools,
1) Python 3.7.3
2) R 3.4.4
3) GATK 3.5
4) bedtools 2.27.1
5) samtools 1.3.1
6) CANOES, CODEX, CLAMMS, XHMM & CN-Learn

The complete list of softwares preinstalled in the image is provided in the [Dockerfile](https://github.com/girirajanlab/CN_Learn/blob/master/Dockerfile).

Docker must be installed using the following commands prior to running CN-Learn. Follow the steps provided in the [instructions](https://docs.docker.com/install/) page for the specific Linux distribution installed on the host machine. 

If docker installation is successful, the following command would return the current version of the docker installed on the host machine.

> **docker version**


# Citation

[A machine-learning approach for accurate detection of copy-number variants from exome sequencing](https://www.biorxiv.org/content/10.1101/460931v2)

Vijay Kumar, Matthew Jensen, Gopal Jayakar, Neil Kelkar, Santhosh Girirajan
bioRxiv 460931; doi: https://doi.org/10.1101/460931


# TLDR
>**`docker run --rm -v ${PROJ_DIR}:${PROJ_DIR} --user $(id -u):$(id -g) girirajanlab/cnlearn \`**

>**`python3 cn_learn.py [DATA_DIRECTORY] [TRAINING_DATA] [TEST_DATA] [CLASSIFIER_TYPE] [LOWER_SIZE_LIMIT] [UPPER_SIZE_LIMIT] [NUMBER_OF_TREES] [CALLER_COUNT] [CALLER_LIST]`**

# Input Files
CN-Learn takes a text file with the list of CNV predictions from multiple algorithms as its primary input and relies on five other important input files. This input files are expected to be placed in their corresponding directories without headers. Sample files for each of the following 6 files are readily provided with the project.


# How to run CN-Learn

### **`Step 1 | Clone CN-Learn repository from Github`** 
Run the following command to clone the required scripts to the local machine.
> **git clone --recursive https://github.com/girirajanlab/CN_Learn.git** 

> **cd CN_Learn**

### **`Step 2 | Place the required input files and update the file locations in config.params file`**
Place the three required input files and update the directory locations in the **config.params** file.

Once the files are available, run the following script to ensure the presence, quality and consistency of the input BAM files, exome capture targets and the reference genome.

> **bash prechecks.sh**

### **`Step 3 | Pull the required docker image`** 
Run the following command to pull the required docker image.

> **docker pull girirajanlab/cnlearn**

> **girirajanlab/cnlearn**

### **`Step 3 | Predict CNVs using CANOES`**
### **`Step 3A:`**
Run the following script to extract the read depth information required by CANOES to make CNV predictions.
> **bash canoes_extract.sh**

### **`Step 3B:`**
Run the following script to make CNV predictions using CANOES.
> **bash canoes_call_CNVs.sh**


### **`Step 4 | Predict CNVs using CLAMMS`**
### **`Step 4A:`**
Run the following script to extract the prerequisite data required by CLAMMS.
> **bash clamms_preprocess.sh**

### **`Step 4B:`**
Run the following script to extract read depth information required by CLAMMS to make CNV predictions.
> **bash clamms_extract.sh**

### **`Step 4C:`**
Run the following script to predict CNVs using CLAMMS.
> **bash clamms_postprocess.sh**


### **`Step 5 | Predict CNVs using CODEX`**
### **`Step 5A:`**
Run the following script to extract the read depth information required by CODEX to make CNV predictions.
> **bash codex_extract.sh**

### **`Step 5B:`**
Run the following script to generate a consolidated output file with the list of CNV calls.
> **bash codex_postprocess.sh**


### **`Step 6 | Predict CNVs using XHMM`**
### **`Step 6A:`**
Run the following script to extract the read depth information required by XHMM to make CNV predictions.
> **bash xhmm_extract.sh**

### **`Step 6B:`**
Run the following script to predict CNVs using XHMM.
> **bash xhmm_call_CNVs.sh**


### **`Step 7:`**
Run calculate_CNV_overlap.sh to measure the CNV overlap among all the callers used.

> **bash calculate_CNV_overlap.sh**


### **`Run either steps 8A & 8B together (Or) just step 8`**

### **`Step 8A:`**
Run generate_bp_coverage.sh to extract the basepair level coverage for each sample. Since this information can be extracted independently for each sample, make the necessary changes to this script to parallelize the process.

> **bash generate_bp_coverage.sh**

### **`Step 8B:`**
Run merge_overlapping_CNVs_readdepth.sh to resolve breakpoint conflicts of concordant CNVs.

> **bash merge_overlapping_CNVs_readdepth.sh**



### **`Step 8:`**
Run merge_overlapping_CNVs_endjoin.sh to resolve breakpoint conflicts of concordant CNVs.

> **bash merge_overlapping_CNVs_endjoin.sh**


### **`Step 9:`**
Run extract_gc_map_vals.sh to extract GC content and mappability scores for singletons and breakpoint-resolved CNVs

> **bash extract_gc_map_vals.sh**

### **`Step 10:`**
Run calc_valdata_overlap.sh to label the training data based on the overlap between CNVs in the training data and the “gold standard” validated CNVs. This script also reformats the CNVs in new samples (i.e., test data).

> **bash calc_valdata_overlap.sh**

### **`Step 11:`**
Run the python script cn_learn.py directly with the required parameters


> **python3 cn_learn.py [DATA_DIRECTORY] [TRAINING_DATA] [TEST_DATA] [CLASSIFIER_TYPE] [LOWER_SIZE_LIMIT] [UPPER_SIZE_LIMIT] [NUMBER_OF_TREES] [CALLER_COUNT] [CALLER_LIST]**

> or

Run cn_learn.sh to train CN-Learn and identify true CNVs in the test set. This script in turn invokes the python script cn_learn.py

> **bash cn_learn.sh**


# Output
Once the above listed steps finish successfully, the final output file named **CNV_list_with_predictions.csv** will be available in the **DATA** directory. The last two columns of the output file provides the probability of the CNVs being 'True' and the classification label **(1 = 'True'; 2 = 'False')** based on a cutoff threshold of **_0.5_**. 
The complete set of columns in the output file is listed below.

> CHR PRED_START PRED_END TYPE SAMPLE NUM_OVERLAPS RD_PROP GC PRED_SIZE MAP NUM_TARGETS SIZE_LABEL **`LIST OF CALLERS`** PRED_PROBS PRED_LABEL


# Copyright/License
    CN-Learn is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CN-Learn is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CN-Learn.  If not, see <https://www.gnu.org/licenses/>.

# Contact
For questions or comments, please contact Vijay Kumar (vxm915@psu.edu), Matthew Jensen (mpj5142@psu.edu) or Santhosh Girirajan (sxg47@psu.edu).
