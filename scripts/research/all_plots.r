#######################################################################################
# Author: Vijay Kumar                                                                 #
# Date  : 5/3/2018                                                                    #
# This is the main script that generates every summary plot required for the project. #
# The first section of the script imports a set of input files while the following    #
# sections use ggplot2 to generate the plots and save them as tiff files.             #
#                                                                                     #
# For readability and ease of use, the input files are listed below                   #
#    1) snp_overlap_data_w_callerinfo.txt - Original set of calls with MA SNP count   #
#    2) cons_overlap_w_caller_0.1.txt - Calls that overlap with atleast 1 SNP M.array #
#    3) cons_overlap_strat1_0.1.txt - Calls consolidated by the first strategy        #
#    4) cons_overlap_strat2_0.1.txt - Calls consolidated by the second strategy       #
#    5) cons_overlap_strat3_0.1.txt - Calls with a particular caller as validation    #
#    6) val_summ_by_sample.txt - List of samples with number of CNVs called by MA     #
#    7) RF_Metrics_SVIP.csv - Classifier metrics for all classifiers                  #
#    8) True_Pos_Calls_SVIP.csv - List of TP calls from first strategy                #
#    9) False_Pos_Calls_SVIP.csv - List of FP calls from third strategy               #
#   10) Samplelist_MA_SVIP.csv - List of samples in the test set (MA)                 #
#   11) Samplelist_CLAMMS_SVIP.csv - List of samples in the test set (CLAMMS)         #
#   12) Overlap_Stats_strat2 - Overlap stats with validated data at different props   #
#   13) Overlap_Stats_strat3 - Overlap stats with validated data at different props   #
#   14) 16p_dups_dels.txt - List of samples with 16p duplication or deletion          #
#   15) 16p_overlap.txt - Subset of calls on Chromosome 16                            #
#   16) 16p_overlap_only.txt - Subset of calls that overlapped with 16p events        #
#   17) 16p_overlap_cnlearn.txt - Subset of 16p overlap calls w/ first strategy       #
#######################################################################################
library(sqldf)
library(ggplot2)
library(reshape2)
library(gplots)
library(limma)
library(Rtsne)
library(ggfortify)
library(rgl)
library(corrplot)
library(Hmisc)

setwd("/Users/vijay/Documents/Projects/Exome_Caller/model")
list.files(getwd())

####################################################################
# Read the files that include CALLER information to generate metrics
####################################################################
original_calls_df <- read.table("snp_overlap_data_w_callerinfo.txt")
colnames(original_calls_df) <- c("CHR", "PRED_START", "PRED_END", "TYPE", "SAMPLE", "CALLER", "PRED_QUAL", "PRED_SIZE", "UNIQ_PRED_ID", 
                                 "CANOES_PROP", "CODEX_PROP", "CLAMMS_PROP", "XHMM_PROP", 
                                 "CANOES_OV", "CODEX_OV", "CLAMMS_OV", "XHMM_OV", "NUM_OVERLAPS", "NUM_TARGETS", "NUM_SNPS")
dim(original_calls_df)
retained_calls_df <- read.table("cons_overlap_w_caller_0.1.txt",  header = TRUE)
retained_calls_df$NUM_OVERLAPS <- retained_calls_df$NUM_OVERLAPS + 1

retained_calls_gt5kb_df <- sqldf("select * from retained_calls_df where PRED_SIZE > 5000")
retained_calls_gt50kb_df <- sqldf("select * from retained_calls_df where PRED_SIZE > 50000")

ma_val_count_df <- read.table("val_summ_by_sample.txt")
colnames(ma_val_count_df) <- c("SAMPLE", "COUNT")

######################################################################
# Read the files prior to and after merging the calls with concordance
######################################################################
all_calls_before_bpres_df <- read.table("cons_calls_w_caller_ov_prop.bed")
colnames(all_calls_before_bpres_df) <- c("CHR", "PRED_START", "PRED_END", "TYPE", "SAMPLE", "CALLER", "PRED_QUAL", "PRED_SIZE", "UNIQ_PRED_ID", 
                                 "CANOES_PROP", "CODEX_PROP", "CLAMMS_PROP", "XHMM_PROP", 
                                 "CANOES_OV", "CODEX_OV", "CLAMMS_OV", "XHMM_OV", "NUM_OVERLAPS")

all_calls_after_bpres_df <- read.table("final_preds_w_rd_stats.txt")
colnames(all_calls_after_bpres_df) <- c("CHR", "PRED_START", "PRED_END", "TYPE", "SAMPLE", "CANOES_PROP", "CODEX_PROP", "CLAMMS_PROP", "XHMM_PROP", 
                                      "CANOES_OV", "CODEX_OV", "CLAMMS_OV", "XHMM_OV", "NUM_OVERLAPS", "GROUP", "RD_PROP")

strat3_calls_after_bpres_df <- read.table("wo_clamms_final_preds_w_rd_stats.txt")
colnames(strat3_calls_after_bpres_df) <- c("CHR", "PRED_START", "PRED_END", "TYPE", "SAMPLE", "CANOES_PROP", "CODEX_PROP", "CLAMMS_PROP", "XHMM_PROP", 
                                        "CANOES_OV", "CODEX_OV", "CLAMMS_OV", "XHMM_OV", "NUM_OVERLAPS", "GROUP", "RD_PROP")
head(strat3_calls_after_bpres_df)
head(all_calls_before_bpres_df)
head(all_calls_after_bpres_df)

val_overlap_stats_strat2_df <- read.table("Overlap_Stats_strat2")
val_overlap_stats_strat3_df <- read.table("Overlap_Stats_strat3")

###############################################################
# Read the file with the list of exome target probe coordinates
###############################################################
target_list <- read.table("targets_auto_no_chr.bed")

#############################################################
# Read the files to generate metrics from consolidated calls
#############################################################
strat1_df <- read.table("cons_overlap_strat1_0.1.txt",  header = TRUE)
strat2_df <- read.table("cons_overlap_strat2_0.1.txt",  header = TRUE)
strat3_df <- read.table("cons_overlap_strat3_0.1.txt",  header = TRUE)

strat2_df$NUM_OVERLAPS <- strat2_df$NUM_OVERLAPS + 1

#############################################################
# Read the files to generate classifier performance measures
#############################################################
rf_metrics_df <- read.csv("RF_Metrics_SVIP.csv")
rf_metrics_size_df <- read.csv("RF_Metrics_Size.csv")


#############################################################
# Read the files to generate true positive distribution
#############################################################
tp_df <- read.csv("True_Pos_Calls_SVIP.csv")
fp_df <- read.csv("False_Pos_Calls_SVIP.csv")
fn_df <- read.csv("False_Neg_Calls_SVIP.csv")
tn_df <- read.csv("True_Neg_Calls_SVIP.csv")

tp_df$CLASS <- 'True Positive'
fp_df$CLASS <- 'False Positive'
fn_df$CLASS <- 'False Negative'
tn_df$CLASS <- 'True Negative'

tp_df$VAL_LABEL <- 'True CNV'
fn_df$VAL_LABEL <- 'True CNV'
tn_df$VAL_LABEL <- 'False CNV'
fp_df$VAL_LABEL <- 'False CNV'


tp_df$NUM_OVERLAPS <- tp_df$NUM_OVERLAPS + 1
fp_df$NUM_OVERLAPS <- fp_df$NUM_OVERLAPS + 1
tn_df$NUM_OVERLAPS <- tn_df$NUM_OVERLAPS + 1
fn_df$NUM_OVERLAPS <- fn_df$NUM_OVERLAPS + 1

tp_fp_df <- rbind(tp_df, fp_df)
tp_fp_fn_df <- rbind(tp_df, fp_df, fn_df)
tp_fp_tn_fn_df <- rbind(tp_df, fp_df, tn_df, fn_df)
tn_fn_df <- rbind(tn_df, fn_df)

tp_fp_tn_fn_df$VAL_LABEL <- ordered(tp_fp_tn_fn_df$VAL_LABEL, levels = c("True CNV", "False CNV"))
tp_fp_tn_fn_df$CLASS <- ordered(tp_fp_tn_fn_df$CLASS, levels = c("True Positive", "False Positive", "True Negative", "False Negative"))

head(tp_df)
####################################################################
# This portion generates the report on the number of 16P deletions #
# and duplications identified by each caller                       #
####################################################################
ma_sample_list_df <- as.data.frame(read.csv("Samplelist_MA_SVIP.csv",header = TRUE))
clamms_sample_list_df <- as.data.frame(read.csv("Samplelist_CLAMMS_SVIP.csv", header = TRUE))

chr16p_samples_df <- as.data.frame(read.table("16p_dups_dels.txt"))
chr16p_ind_callers_df <- as.data.frame(read.table("16p_overlap_only.txt"))
chr16p_cnlearn_df <- as.data.frame(read.table("16p_overlap_cnlearn.txt"))

colnames(chr16p_samples_df) <- c("CHR", "START", "END", "SAMPLE", "CNV_TYPE")
colnames(chr16p_ind_callers_df) <- c("CHR", "PRED_START", "PRED_END", "CNV_TYPE", "SAMPLE", "CALLER", "PRED_QUAL", "PRED_SIZE", "UNIQ_PRED_ID",  "CANOES_PROP", "CODEX_PROP", 
                                  "CLAMMS_PROP", "XHMM_PROP", "CANOES_OV", "CODEX_OV", "CLAMMS_OV", "XHMM_OV", "NUM_OVERLAPS",
                                  "NUM_TARGETS", "NUM_SNPS", "OV_STATUS")

colnames(chr16p_cnlearn_df) <- c("CHR", "PRED_START", "PRED_END", "CNV_TYPE", "SAMPLE", "CANOES_PROP", "CODEX_PROP", "CLAMMS_PROP", "XHMM_PROP", "NUM_OVERLAPS", "RD_PROP",
                                  "GC_CONTENT", "PRED_SIZE", "MAPPABILITY","NUM_TARGETS","NUM_SNPS", "OVERLAP_PROP", "SIZE_LABEL", "PRED_PROBS", "CLASSIFIER_NAME", "NUM_SPLITS", "STRATEGY", "TRAINING_PROP", "OV_STATUS")

##########################################################
# Read the files needed to generate supplementary tables #
##########################################################
setwd("/Users/vijay/Documents/Projects/Exome_Caller/model/")
sample_all_intvls <- read.table("14734.x11_DEL_preds_all_intvl_combos.txt")
colnames(sample_all_intvls) <- c("CHR", "PRED_START", "PRED_END", "TYPE", "SAMPLE", "CALLER_ST", "CALLER_END", "UNIQ_PRED_ST", "UNIQ_PRED_END", 
                                 "CANOES_PROP", "CODEX_PROP", "CLAMMS_PROP", "XHMM_PROP", 
                                 "CANOES_OV", "CODEX_OV", "CLAMMS_OV", "XHMM_OV", "NUM_OVERLAPS", "GROUP_NAME")

########################################
# Global variables for subsequent uses #
########################################
venn_colors <- c("red", "blue", "green", "orchid")

setwd("/Users/vijay/Documents/Projects/Exome_Caller/figures_and_tables")#######################################################################################
# Author: Vijay Kumar                                                                 #
# Date  : 5/3/2018                                                                    #
# This is the main script that generates every summary plot required for the project. #
# The first section of the script imports a set of input files while the following    #
# sections use ggplot2 to generate the plots and save them as tiff files.             #
#                                                                                     #
# For readability and ease of use, the input files are listed below                   #
#    1) snp_overlap_data_w_callerinfo.txt - Original set of calls with MA SNP count   #
#    2) cons_overlap_w_caller_0.1.txt - Calls that overlap with atleast 1 SNP M.array #
#    3) cons_overlap_strat1_0.1.txt - Calls consolidated by the first strategy        #
#    4) cons_overlap_strat2_0.1.txt - Calls consolidated by the second strategy       #
#    5) cons_overlap_strat3_0.1.txt - Calls with a particular caller as validation    #
#    6) val_summ_by_sample.txt - List of samples with number of CNVs called by MA     #
#    7) RF_Metrics_SVIP.csv - Classifier metrics for all classifiers                  #
#    8) True_Pos_Calls_SVIP.csv - List of TP calls from first strategy                #
#    9) False_Pos_Calls_SVIP.csv - List of FP calls from third strategy               #
#   10) Samplelist_MA_SVIP.csv - List of samples in the test set (MA)                 #
#   11) Samplelist_CLAMMS_SVIP.csv - List of samples in the test set (CLAMMS)         #
#   12) Overlap_Stats_strat2 - Overlap stats with validated data at different props   #
#   13) Overlap_Stats_strat3 - Overlap stats with validated data at different props   #
#   14) 16p_dups_dels.txt - List of samples with 16p duplication or deletion          #
#   15) 16p_overlap.txt - Subset of calls on Chromosome 16                            #
#   16) 16p_overlap_only.txt - Subset of calls that overlapped with 16p events        #
#   17) 16p_overlap_cnlearn.txt - Subset of 16p overlap calls w/ first strategy       #
#######################################################################################
library(sqldf)
library(ggplot2)
library(reshape2)
library(gplots)
library(limma)
library(Rtsne)
library(ggfortify)
library(rgl)
library(corrplot)
library(Hmisc)

setwd("/Users/vijay/Documents/Projects/Exome_Caller/model")
list.files(getwd())

####################################################################
# Read the files that include CALLER information to generate metrics
####################################################################
original_calls_df <- read.table("snp_overlap_data_w_callerinfo.txt")
colnames(original_calls_df) <- c("CHR", "PRED_START", "PRED_END", "TYPE", "SAMPLE", "CALLER", "PRED_QUAL", "PRED_SIZE", "UNIQ_PRED_ID", 
                                 "CANOES_PROP", "CODEX_PROP", "CLAMMS_PROP", "XHMM_PROP", 
                                 "CANOES_OV", "CODEX_OV", "CLAMMS_OV", "XHMM_OV", "NUM_OVERLAPS", "NUM_TARGETS", "NUM_SNPS")
dim(original_calls_df)
retained_calls_df <- read.table("cons_overlap_w_caller_0.1.txt",  header = TRUE)
retained_calls_df$NUM_OVERLAPS <- retained_calls_df$NUM_OVERLAPS + 1

retained_calls_gt5kb_df <- sqldf("select * from retained_calls_df where PRED_SIZE > 5000")
retained_calls_gt50kb_df <- sqldf("select * from retained_calls_df where PRED_SIZE > 50000")

ma_val_count_df <- read.table("val_summ_by_sample.txt")
colnames(ma_val_count_df) <- c("SAMPLE", "COUNT")

######################################################################
# Read the files prior to and after merging the calls with concordance
######################################################################
all_calls_before_bpres_df <- read.table("cons_calls_w_caller_ov_prop.bed")
colnames(all_calls_before_bpres_df) <- c("CHR", "PRED_START", "PRED_END", "TYPE", "SAMPLE", "CALLER", "PRED_QUAL", "PRED_SIZE", "UNIQ_PRED_ID", 
                                 "CANOES_PROP", "CODEX_PROP", "CLAMMS_PROP", "XHMM_PROP", 
                                 "CANOES_OV", "CODEX_OV", "CLAMMS_OV", "XHMM_OV", "NUM_OVERLAPS")

all_calls_after_bpres_df <- read.table("final_preds_w_rd_stats.txt")
colnames(all_calls_after_bpres_df) <- c("CHR", "PRED_START", "PRED_END", "TYPE", "SAMPLE", "CANOES_PROP", "CODEX_PROP", "CLAMMS_PROP", "XHMM_PROP", 
                                      "CANOES_OV", "CODEX_OV", "CLAMMS_OV", "XHMM_OV", "NUM_OVERLAPS", "GROUP", "RD_PROP")

strat3_calls_after_bpres_df <- read.table("wo_clamms_final_preds_w_rd_stats.txt")
colnames(strat3_calls_after_bpres_df) <- c("CHR", "PRED_START", "PRED_END", "TYPE", "SAMPLE", "CANOES_PROP", "CODEX_PROP", "CLAMMS_PROP", "XHMM_PROP", 
                                        "CANOES_OV", "CODEX_OV", "CLAMMS_OV", "XHMM_OV", "NUM_OVERLAPS", "GROUP", "RD_PROP")
head(strat3_calls_after_bpres_df)
head(all_calls_before_bpres_df)
head(all_calls_after_bpres_df)

val_overlap_stats_strat2_df <- read.table("Overlap_Stats_strat2")
val_overlap_stats_strat3_df <- read.table("Overlap_Stats_strat3")

###############################################################
# Read the file with the list of exome target probe coordinates
###############################################################
target_list <- read.table("targets_auto_no_chr.bed")

#############################################################
# Read the files to generate metrics from consolidated calls
#############################################################
strat1_df <- read.table("cons_overlap_strat1_0.1.txt",  header = TRUE)
strat2_df <- read.table("cons_overlap_strat2_0.1.txt",  header = TRUE)
strat3_df <- read.table("cons_overlap_strat3_0.1.txt",  header = TRUE)

strat2_df$NUM_OVERLAPS <- strat2_df$NUM_OVERLAPS + 1

#############################################################
# Read the files to generate classifier performance measures
#############################################################
rf_metrics_df <- read.csv("RF_Metrics_SVIP.csv")
rf_metrics_size_df <- read.csv("RF_Metrics_Size.csv")


#############################################################
# Read the files to generate true positive distribution
#############################################################
tp_df <- read.csv("True_Pos_Calls_SVIP.csv")
fp_df <- read.csv("False_Pos_Calls_SVIP.csv")
fn_df <- read.csv("False_Neg_Calls_SVIP.csv")
tn_df <- read.csv("True_Neg_Calls_SVIP.csv")

tp_df$CLASS <- 'True Positive'
fp_df$CLASS <- 'False Positive'
fn_df$CLASS <- 'False Negative'
tn_df$CLASS <- 'True Negative'

tp_df$VAL_LABEL <- 'True CNV'
fn_df$VAL_LABEL <- 'True CNV'
tn_df$VAL_LABEL <- 'False CNV'
fp_df$VAL_LABEL <- 'False CNV'


tp_df$NUM_OVERLAPS <- tp_df$NUM_OVERLAPS + 1
fp_df$NUM_OVERLAPS <- fp_df$NUM_OVERLAPS + 1
tn_df$NUM_OVERLAPS <- tn_df$NUM_OVERLAPS + 1
fn_df$NUM_OVERLAPS <- fn_df$NUM_OVERLAPS + 1

tp_fp_df <- rbind(tp_df, fp_df)
tp_fp_fn_df <- rbind(tp_df, fp_df, fn_df)
tp_fp_tn_fn_df <- rbind(tp_df, fp_df, tn_df, fn_df)
tn_fn_df <- rbind(tn_df, fn_df)

tp_fp_tn_fn_df$VAL_LABEL <- ordered(tp_fp_tn_fn_df$VAL_LABEL, levels = c("True CNV", "False CNV"))
tp_fp_tn_fn_df$CLASS <- ordered(tp_fp_tn_fn_df$CLASS, levels = c("True Positive", "False Positive", "True Negative", "False Negative"))

head(tp_df)
####################################################################
# This portion generates the report on the number of 16P deletions #
# and duplications identified by each caller                       #
####################################################################
ma_sample_list_df <- as.data.frame(read.csv("Samplelist_MA_SVIP.csv",header = TRUE))
clamms_sample_list_df <- as.data.frame(read.csv("Samplelist_CLAMMS_SVIP.csv", header = TRUE))

chr16p_samples_df <- as.data.frame(read.table("16p_dups_dels.txt"))
chr16p_ind_callers_df <- as.data.frame(read.table("16p_overlap_only.txt"))
chr16p_cnlearn_df <- as.data.frame(read.table("16p_overlap_cnlearn.txt"))

colnames(chr16p_samples_df) <- c("CHR", "START", "END", "SAMPLE", "CNV_TYPE")
colnames(chr16p_ind_callers_df) <- c("CHR", "PRED_START", "PRED_END", "CNV_TYPE", "SAMPLE", "CALLER", "PRED_QUAL", "PRED_SIZE", "UNIQ_PRED_ID",  "CANOES_PROP", "CODEX_PROP", 
                                  "CLAMMS_PROP", "XHMM_PROP", "CANOES_OV", "CODEX_OV", "CLAMMS_OV", "XHMM_OV", "NUM_OVERLAPS",
                                  "NUM_TARGETS", "NUM_SNPS", "OV_STATUS")

colnames(chr16p_cnlearn_df) <- c("CHR", "PRED_START", "PRED_END", "CNV_TYPE", "SAMPLE", "CANOES_PROP", "CODEX_PROP", "CLAMMS_PROP", "XHMM_PROP", "NUM_OVERLAPS", "RD_PROP",
                                  "GC_CONTENT", "PRED_SIZE", "MAPPABILITY","NUM_TARGETS","NUM_SNPS", "OVERLAP_PROP", "SIZE_LABEL", "PRED_PROBS", "CLASSIFIER_NAME", "NUM_SPLITS", "STRATEGY", "TRAINING_PROP", "OV_STATUS")

##########################################################
# Read the files needed to generate supplementary tables #
##########################################################
setwd("/Users/vijay/Documents/Projects/Exome_Caller/model/")
sample_all_intvls <- read.table("14734.x11_DEL_preds_all_intvl_combos.txt")
colnames(sample_all_intvls) <- c("CHR", "PRED_START", "PRED_END", "TYPE", "SAMPLE", "CALLER_ST", "CALLER_END", "UNIQ_PRED_ST", "UNIQ_PRED_END", 
                                 "CANOES_PROP", "CODEX_PROP", "CLAMMS_PROP", "XHMM_PROP", 
                                 "CANOES_OV", "CODEX_OV", "CLAMMS_OV", "XHMM_OV", "NUM_OVERLAPS", "GROUP_NAME")

########################################
# Global variables for subsequent uses #
########################################
venn_colors <- c("red", "blue", "green", "orchid")

setwd("/Users/vijay/Documents/Projects/Exome_Caller/figures_and_tables")
