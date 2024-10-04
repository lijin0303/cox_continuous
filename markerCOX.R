#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
# args <- c("HPFS","alcohol","lineavgr","dt_crc","time_crc",
#           "NULL",3,
#           "list(interval,ageyrs)","HPFS_k3",0)

require(Hmisc)
require(optimx)
require(data.table)
require(rms)
require(ggplot2)
require(flextable)
require(dplyr)
require(MASS)
require(officer)

setwd("/udd/nhrli/CHANNING")
source("Code/Main.R")
User_data <- switch(args[1],
                    "NHS"=get(load("Data/Processed/NHS.RData")),
                    "HPFS"=get(load("Data/Processed/HPFS.RData")),
                    "both"=get(load("Data/Processed/both.RData")))
User_data <- datproc(User_data)
exposure_var <- args[2]
marker_var <- args[3]
case_var <- args[4]
time_var <- args[5]
covariates <- args[6]
Knots <- eval(parse(text=args[7]))
Groups <- args[8]
outdir <- args[9]
strt <- eval(parse(text=args[10]))

### Main Optimization Function
if (covariates=="NULL"){
  n <- Knots
  cov_strg <- c("Intercept",sapply(seq(Knots-1),function(v) paste("Z",v,sep="")))
}else{
    n <- Knots + length(strsplit(covariates,split=",")[[1]])
    cov_strg <-  c("Intercept",sapply(seq(Knots-1),function(v) paste("Z",v,sep="")),
                   strsplit(gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.","",
                                 x=covariates,perl=T),split=",")[[1]])}

opt <-optim(par = rep(strt,n),NegLoglik_spline,derv_spline,
            method ="BFGS",
            dat = User_data,
            main=exposure,
            cov=covariates,
            case_ind = case_var,
            marker = marker_var,
            time = time_var,
            by_var=Groups,
            k = n,
            hessian=T,
            control=list(ndeps=rep(1e-5,n), maxit=800, trace=T, REPORT = 50, reltol=1e-10))
save(opt,file=paste("./","optim_",outdir,".RData",sep=""))

### Extract significance information
par_intv_info <- uncertain_info(opt)
par_intv_info$Variable <- cov_strg
par_intv_info <- par_intv_info[,c("Variable",colnames(par_intv_info)[-dim(par_intv_info)[2]])]
fwrite(par_intv_info,file=paste(outdir,"/","Parameters_",outdir,".txt",sep=""))
