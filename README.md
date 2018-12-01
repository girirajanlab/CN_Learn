# CN-Learn 
CN-Learn is a framework to integrate Copy Number Variant (CNV) predictions made by multiple algorithms using exome sequencing datasets. Unlike traditional methods that depend on a single measure of concordance, CN-Learn leverages the extent of concordance among multiple method as well as additional genomic contexts such as GC content, mappability and local read depth fluctuations. Using this exhaustive set of predictors extracted for CNVs from a “gold-standard” CNV call set, CN-Learn builds a ‘Random Forest’ classifier with hundreds of decision trees to estimate the probability of each CNV in a fresh set of samples being true. This strategy provides CN-Learn the ability to segregate true CNV calls from false positives in an informed fashion with extremely high precision.

While CN-Learn has been shown to perform best when built as a ‘Random Forest’ classifier, it can also be built as a ‘Logistic Regression’ or ‘Support Vector Machine’ classifier. CN-Learn can also be seamlessly extended to include newer set of CNV calling algorithms in the future by simply changing a single parameter supplied to CN-Learn. 

## Software Requirements
1.	python – Version 2.7
2.	R – Version 3.4.4
3.	Required python packages can be installed using `pip install python_packages_prereq.txt`
4.	Required R packages – sqldf (v0.4.11), reshape2 (v1.4.3).
5.	GATK (v3.5)
6.	bedtools (v2.0)
7.	samtools (v1.3.1)

# Citation
Kumar et al, "An integrative approach for accurate detection of copy-number variants from exome sequencing"

# TLDR
>**`python cn_learn.py [DATA_DIRECTORY] [TRAINING_DATA] [TEST_DATA] [CLASSIFIER_TYPE] [LOWER_SIZE_LIMIT] [UPPER_SIZE_LIMIT] [NUMBER_OF_TREES] [CALLER_COUNT] [CALLER_LIST]`**

### Directory and Data Organization
Once the ‘CN_Learn’ repository is cloned to a local folder, all the directories and subdirectories referenced in the scripts used across the project will be readily available. The ‘scripts’ directory holds all the scripts used in the project, while ‘data’ and ‘source’ are the two data directories. The ‘data’ directory acts as the working directory to hold the primary input file along with other intermediate and output files. All other reference/input files (bam files, gold-standard CNV list etc.) are placed in the ‘source’ directory.

### Script Organization
All the scripts required to run CN-Learn are provided in the ‘scripts’ directory. There are five bash scripts that preprocess data required for CN-Learn and one that executes CN-Learn. These bash scripts in turn invoke R and Python scripts to accomplish specific tasks. 

### config.params 
In order to make the workflow/scripts easy to read and execute, several commonly used parameters, directory names and file paths are maintained in a single parameter file named ‘config.params’. This file is located in the main project folder and sourced into each bash script using the ‘source’ command. 

NOTE: Since the config.params file is sourced into every script using hard-coded full path, the ‘source’ command in every bash script **MUST** be updated with the full local path.

# Input Files
CN-Learn takes a text file with the list of CNV predictions from multiple algorithms as its primary input and relies on five other important input files. This input files are expected to be placed in their corresponding directories without headers. Sample files for each of the following 6 files are readily provided with the project.

Place the following file in the **DATA** directory.

1. **consolidated_calls.bed** = List of CNV predictions made by multiple algorithms. 
>**`CHROMOSOME	START	END	CNV_TYPE	SAMPLE		CALLER`**

>**`CNV_TYPE = DUP or DEL`**

Place the following files in the **SOURCE** directory. 

2. **sample_list.txt**	= List of samples to be processed by CN-Learn
3. **sample_list_train.txt** = List of samples to be used to train CN-Learn
4. **sample_list_test.txt**	= List of samples in which true CNVs are to be identified by CN-Learn
5. **targets_auto_no_chr.bed** = List of exome capture probes in the format below.
>**`CHR	START	END`**
6. **validated_cnvs.txt** = List of “gold standard” CNVs in the format below.
>**`CHR	START	END	CNV_TYPE	SIZE	SAMPLE`**

# How to run CN-Learn
### **`Step 1:`** 
Make sure all the prerequisite softwares and packages are installed. Clone the ‘CN_Learn’ repository to a local directory.

### **`Step 2:`**
Make the following updates to the parameters listed in the ‘config.params’ file,
1.	Update the PROJ_DIR parameter with the full path of the local project directory. 
2.	Update the CALLER_LIST parameter to list all the CNV callers used in the project.
3.	Review and update the VALDATA_OV_THRESHOLD.
4.	Update the software tool and file locations.

Additional instructions are provided in the config.params file. 

Once the changes are complete, update the project directory path provided in the ‘source’ command in each bash script (`ls *.sh`) in the ‘scripts’ directory.

### **`Step 3:`**
Make sure that all the six input files as well as the bam files for each sample are present in their corresponding directories.


### **`Step 4:`**
Run calculate_CNV_overlap.sh to measure the CNV overlap among all the callers used.

> **bash calculate_CNV_overlap.sh**

<br/><br/>

### **`Run either steps 5A & 5B together (Or) just step 5`**

### **`Step 5A:`**
Run generate_bp_coverage.sh to extract the basepair level coverage for each sample. Since this information can be extracted independently for each sample, make the necessary changes to this script to parallelize the process.

> **bash generate_bp_coverage.sh**

### **`Step 5B:`**
Run merge_overlapping_CNVs_readdepth.sh to resolve breakpoint conflicts of concordant CNVs.

> **bash merge_overlapping_CNVs_readdepth.sh**

### **`Step 5:`**
Run merge_overlapping_CNVs_endjoin.sh to resolve breakpoint conflicts of concordant CNVs.

> **bash merge_overlapping_CNVs_endjoin.sh**

<br/><br/>

### **`Step 6:`**
Run extract_gc_map_vals.sh to extract GC content and mappability scores for singletons and breakpoint-resolved CNVs

> **bash extract_gc_map_vals.sh**

### **`Step 7:`**
Run calc_valdata_overlap.sh to label the training data based on the overlap between CNVs in the training data and the “gold standard” validated CNVs. This script also reformats the CNVs in new samples (i.e., test data).

> **bash calc_valdata_overlap.sh**

### **`Step 8:`**
Run the python script cn_learn.py directly with the required parameters


> **python cn_learn.py [DATA_DIRECTORY] [TRAINING_DATA] [TEST_DATA] [CLASSIFIER_TYPE] [LOWER_SIZE_LIMIT] [UPPER_SIZE_LIMIT] [NUMBER_OF_TREES] [CALLER_COUNT] [CALLER_LIST]**

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
