# CN-Learn 
CN-Learn is a framework to integrate Copy Number Variant (CNV) predictions made by multiple algorithms using exome sequencing datasets. CN-Learn can be modified and/or extended to include newer algorithms in the future. Using a small set of CNVs with cross-platform validations, CN-Learn learns to segregate true CNVs from false positives. 

## Software Requirements
1.	python 2.7
2.	Required python packages can be installed using `pip install python_packages_prereq.txt`
3.	GATK (v3.5)
4.	bedtools
5.	samtools (v1.3.1)

## TLDR
>**`python cn_learn.py [DATA_DIRECTORY] [TRAINING_DATA] [TEST_DATA] [CLASSIFIER_TYPE] [LOWER_SIZE_LIMIT] [UPPER_SIZE_LIMIT] [NUMBER_OF_TREES] [CALLER_COUNT] [CALLER_LIST]`**

## Script Organization
All the scripts required to run CN-Learn are provided in the scripts directory. There are five bash scripts that preprocess data required for CN-Learn and one that executes CN-Learn. These bash scripts in turn invoke R and Python scripts to accomplish specific tasks. The list of directories, files and global parameters are provided in the config.params file. This configuration file is sourced into each script to make the values accessible within each script.

## Directory and Data Organization
Once the ‘scripts’ directory is cloned to a local folder for a project, the first bash script is designed to automatically create all the directories and subdirectories. ‘data’ and ‘source’ are the two main directories used throughout the project, where the ‘data’ directory acts as the working directory to hold intermediate and output files while the input files are placed in the ‘source’ directory.

# Input Files
CN-Learn takes a text file with the list of CNV predictions from multiple algorithms as its primary input and relies on five other important input files. This input files are expected to be without headers.

1. **consolidated_calls.bed** = List of CNV predictions made by multiple algorithms. 
>**`CHROMOSOME	START	END	CNV_TYPE	SAMPLE		CALLER`**

>**`CNV_TYPE = DUP or DEL`**

Generate the following files required by CN-Learn without headers. Place them in the SOURCE directory. 
2. **sample_list.txt**	= List of samples to be processed by CN-Learn
3. **sample_list_train.txt** = List of samples to be used to train CN-Learn
4. **sample_list_test.txt**	= List of samples in which true CNVs are to be identified by CN-Learn
5. **targets_auto_no_chr.bed** = List of exome capture probes in the format below.
>**`CHR	START	END`**

6. **validated_cnvs.txt** = List of “gold standard” CNVs in the format below.
>**`CHR	START	END	CNV_TYPE	SIZE	SAMPLE`**


## How to run CN-Learn
**Step 1:** 
Make sure all the required software packages are installed. Clone the ‘scripts’ folder and download the config.params file into a local project directory.

**Step 2:**
Supply the local project directory location to the PROJ_DIR parameter in the config.params file. Make sure that the software tool locations, file locations and the global variables are initialized appropriately in the config.params file. Additional instructions are provided in the config.params file.

**Step 3:**
Run create_init_dirs.sh to generate the required directories and subdirectories inside the project directory.
```
bash create_init_dirs.sh
```

**Step 4:**
Run calculate_CNV_overlap.sh to measure the CNV overlap among all the the callers used.
```
bash calculate_CNV_overlap.sh
```

**Step 5:**
Run merge_overlapping_CNVs.sh to resolve breakpoint conflicts of concordant CNVs.
```
bash merge_overlapping_CNVs.sh
```

**Step 6:**
Run extract_gc_map_vals.sh to extract GC content and mappability scores for singletons and breakpoint-resolved CNVs
```
bash extract_gc_map_vals.sh
```
**Step 7:**
Run calc_valdata_overlap.sh to label the training data based on the overlap between CNVs in the training data and the “gold standard” validated CNVs. This script also reformats the CNVs in new samples (i.e., test data).
```
bash calc_valdata_overlap.sh
```
**Step 8:**
Run the python script cn_learn.py with the required input parameters to train and make predictions on CNVs in new samples/
```
python cn_learn.py [DATA_DIRECTORY] [TRAINING_DATA] [TEST_DATA] [CLASSIFIER_TYPE] [LOWER_SIZE_LIMIT] [UPPER_SIZE_LIMIT] 
                    [NUMBER_OF_TREES] [CALLER_COUNT] [CALLER_LIST]
```

> or

Run cn_learn.sh to train CN-Learn and identify true CNVs in the test set. This script in turn invokes the python script cn_learn.py.
```
bash cn_learn.sh
```


