# CN-Learn 
CN-Learn is a framework to integrate Copy Number Variant (CNV) predictions made by multiple algorithms using exome sequencing datasets. Unlike traditional integration methods that depend on a single measure of concordance, CN-Learn leverages the extent of concordance among the calls made by multiple methods, in addition to genomic contexts such as GC content, mappability and local read depth fluctuations. Using a wide range of predictors extracted for CNVs from a “gold-standard” CNV call set, CN-Learn builds a ‘Random Forest’ classifier with hundreds of decision trees to estimate the probability of each CNV in a fresh set of samples being true. This strategy provides CN-Learn the ability to segregate true CNV calls from false positives in an informed fashion with extremely high precision.

While CN-Learn has been shown to perform best when built as a ‘Random Forest’ classifier, it can also be built as a ‘Logistic Regression’ or ‘Support Vector Machine’ classifier. CN-Learn can also be seamlessly extended to include newer set of CNV calling algorithms in the future by simply changing a single parameter supplied to CN-Learn. 

# Citation

[A machine-learning approach for accurate detection of copy-number variants from exome sequencing](https://www.biorxiv.org/content/10.1101/460931v2)

Vijay Kumar, Matthew Jensen, Gopal Jayakar, Neil Kelkar, Santhosh Girirajan
bioRxiv 460931; doi: https://doi.org/10.1101/460931


# Software Requirements
Given the number of softwares/tools required to run the four individual CNV calling algorithms prior to running CN-Learn, every software/tool required to run CN-Learn end-to-end has been packaged into a Docker image. If the user chooses to make use of Docker to run CN-Learn, [Docker](https://www.docker.com/) must first be installed on the host machine. Please follow the steps provided in the [instructions](https://docs.docker.com/install/) page to install Docker for the specific Linux distribution installed on the host machine. If the installation is successful, the following command will return the current version of docker installed on the host machine.

> **docker version**

Once docker is available, the image can then be downloaded and made instantly available for use on the host machine using the following command.

>**docker pull girirajanlab/cnlearn**

Following are some of the software tools preinstalled on the docker image,
1) Python 3.7.3
2) R 3.4.4
3) GATK 3.5
4) bedtools 2.27.1
5) samtools 1.3.1
6) CANOES, CODEX, CLAMMS, XHMM & CN-Learn

The complete list of preinstalled softwares can be found in the [Dockerfile](https://github.com/girirajanlab/CN_Learn/blob/master/Dockerfile).


# Logistics
Running CN-Learn to identify CNVs involves the following tasks,

**1) Github:** Clone the CN_Learn github repo to a local host LINUX machine using the following command,

> **git clone --recursive https://github.com/girirajanlab/CN_Learn.git** 

**2) BAM files:** Place all the BAM files, along with their corresponding index files in a local directory. Ensure the following,

	a) All the bam files should be named <SAMPLE>.bam and the index file named <SAMPLE>.bam.bai, 
    where <SAMPLE> is the name of the sample without any special characters in them.
    
    b) Each bam file must have an index file associated with it.
    
    c) The directory with .bam and .bam.bai files should not have any other type of files in them.
    
**3) Reference genome:** Make sure that the version of reference genome to which the samples were mapped to, is available in a local directory, along with the index files. In addition to **<REFERENCE_GENOME>.fasta**, the following files must also be present in the same directory,
   
    a) <REFERENCE_GENOME>.fasta.fai
    
    b) <REFERENCE_GENOME>.dict

**4) Exome capture probes:** Name the file with the list of exome capture probes as **exome_capture_targets.bed** and place the file in the **/source/** directory inside the CN_Learn repository that was just cloned. 

**Important Note:** Make sure that the file is **tab separated** with the first three columns being Chromosome, Start Position and End Position.

**5) List of validated CNVs:** Place the file named **validated_cnvs.txt** in the source directory.

**6) config.params:** Update the following parameters in the config.params file in the CN_Learn directory that was just cloned;

	a) BAM_FILE_DIR     : Replace 'TBD' with the full path of the directory with all the BAM files.
    
    b) REF_GENOME       : Replace 'TBD' with the full path of the reference genome file. 
    
    c) SW_DIR           : This path is set to the directory inside the Docker image. If you are NOT 
                          using docker, please update this path to the location of the directory in 
                          the local file system.
    
    d) DOCKER_INDICATOR : This parameter is set to 'Y' by default. If you choose NOT to use Docker 
                          and prefer to use locally installed softwares, pelase update this parameter 
                          to 'N' prior to running rest of the steps. 
    
    
**6) Docker:** If you decide to use docker, download the image using the following command,
>**docker pull girirajanlab/cnlearn**

Run the following command and make sure that it lists the recently downloaded image,
> **docker images**

**8) Prechecks:** Once all the input files are available, run the following script to ensure the presence, quality and consistency of the input BAM files, exome capture targets and the reference genome.

> **bash prechecks.sh** 

**9)** Once the prechecks.sh executes successfully without errors, follow the steps below to generate CNVs.

## How to run CN-Learn?


### **`Step 1 | Predict CNVs using CANOES`**
### **`Step 1A:`**
Run the following script to extract the read depth information required by CANOES to make CNV predictions.
> **bash canoes_extract.sh**

### **`Step 1B:`**
Run the following script to make CNV predictions using CANOES.
> **bash canoes_call_CNVs.sh**


### **`Step 2 | Predict CNVs using CLAMMS`**
### **`Step 2A:`**
Run the following script to extract the prerequisite data required by CLAMMS.
> **bash clamms_preprocess.sh**

### **`Step 2B:`**
Run the following script to extract read depth information required by CLAMMS to make CNV predictions.
> **bash clamms_extract.sh**

### **`Step 2C:`**
Run the following script to predict CNVs using CLAMMS.
> **bash clamms_postprocess.sh**


### **`Step 3 | Predict CNVs using CODEX`**
### **`Step 3A:`**
Run the following script to extract the read depth information required by CODEX to make CNV predictions.
> **bash codex_extract.sh**

### **`Step 3B:`**
Run the following script to generate a consolidated output file with the list of CNV calls.
> **bash codex_postprocess.sh**


### **`Step 4 | Predict CNVs using XHMM`**
### **`Step 4A:`**
Run the following script to extract the read depth information required by XHMM to make CNV predictions.
> **bash xhmm_extract.sh**

### **`Step 4B:`**
Run the following script to predict CNVs using XHMM.
> **bash xhmm_call_CNVs.sh**


### **`Step 5 | Measure the overlap among callers`**
Run calculate_CNV_overlap.sh to measure the CNV overlap among all the callers used.

> **bash calculate_CNV_overlap.sh**


### **`Run either steps 6A & 6B together (Or) just step 6`**

### **`Step 6A | Extract basepair level coverage info`**
Run generate_bp_coverage.sh to extract the basepair level coverage for each sample. Since this information can be extracted independently for each sample, make the necessary changes to this script to parallelize the process.

> **bash generate_bp_coverage.sh**

### **`Step 6B | Resolve breakpoints`**
Run merge_overlapping_CNVs_readdepth.sh to resolve breakpoint conflicts of concordant CNVs.

> **bash merge_overlapping_CNVs_readdepth.sh**


### **`Step 6 | Resolve breakpoints`**
Run merge_overlapping_CNVs_endjoin.sh to resolve breakpoint conflicts of concordant CNVs.

> **bash merge_overlapping_CNVs_endjoin.sh**


### **`Step 7 | Extract GC content and mappability in breakpoint-resolved CNV regions`**
Run extract_gc_map_vals.sh to extract GC content and mappability scores for singletons and breakpoint-resolved CNVs

> **bash extract_gc_map_vals.sh**

### **`Step 8 | Label CNVs based on gold-standard validations`**
Run calc_valdata_overlap.sh to label the training data based on the overlap between CNVs in the training data and the “gold standard” validated CNVs. This script also reformats the CNVs in new samples (i.e., test data).

> **bash calc_valdata_overlap.sh**

### **`Step 9 | Classify CNVs`**
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
