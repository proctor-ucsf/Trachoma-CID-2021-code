CLEAN_YRS <- 2004:2019
#DISTRIBUTIONS <- c("renorm_beta","renorm_exp","renorm_gamma","renorm_pareto")
DISTRIBUTIONS <- c("renorm_exp","exp_norm_mix")
ZERO_THRESH <- 0.5
DISPLAY_WIDTH <- 5 # Width of histograms displayed
KEY_COL_NAMES <- c("survey_year","tf_prev")
BIN_WIDTH <- 0.5
SEED <- 1
NUM_RESAMP <- 500 # For bootsraps
COHORT_SIZE <- 4
EXP_THRESH <- 100

set.seed(1)

# library ("fitdistrplus")
library("gridExtra")
library("tictoc")
library("tidyverse")

setwd("/Users/sblumberg/Google Drive/Research/Trachoma/Trachoma prediction code/Global nowcast 091620/")
source("../forecast_fxns.R")
set.seed(SEED)

# To use rigel: ssh rigel / scp
# To initiate screen: screen - S name; ctrl-A ctrl-D to detach; screen -x name to resume
# Rscript name.R to run from script; To stop: Ctrl-c from terminal, or k from 'top'14
