
rm(list = ls(all = TRUE))
args <- commandArgs(TRUE)




source("get_spatial_cov.R")
source("spattempcl.R")
source("get_preds.R")
load("all_info_forpred.RData")



fn1 <- args[[1]]
rdatapath <- args[[2]]
outpath <- args[[3]]


#order is 
#1. Regional/matern
#2. Regional/exp
#3. Nnal/matern
#4. Nnal/exp

type <- NULL
if(k < 3) {
	type <- "reg"
}

# detrendmod
# detrendmodEXP
# detrendmodNNAL
# detrendmodNNALexp

# detrendmodPM
# detrendmodPMexp
# detrendmodPMNNAL
# detrendmodPMNNALexp

suffs <- c("", "EXP", "NNAL", "NNALexp")
suf1 <- suffs[k]
PARsufcons <- paste0("detrendmod", suf1)

if(k == 2) {suf1 <- tolower(suf1)}
PARsufPM <- paste0("detrendmodPM", suf1)



getpreds(fn1, PARsufcons, PARsufPM, rdatapath, outpath, type)