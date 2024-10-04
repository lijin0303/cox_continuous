## Utility function definition
caseinfo <- function(Dat){
  Dat <- Dat[order(time_case)]
  firstoccur <- Dat[!duplicated(time_case),.(id,time_case)]
  firstoccur[,index := match(id,Dat$id)]
  Dat <- merge(Dat,firstoccur[,.(time_case,index)],by="time_case")
  return(Dat)
}

NegLoglik_spline <-function(par,dat,main,cov,case_ind,marker,time,k,by_var="list(interval)"){
  #' This is a function to compute negative log partial likelihood for given dataset
  #' With Brewlow method implemented for handling ties
  #' @param par The estimated coefficient values (numeric vector)
  #' @param dat The given dataset (data.table) 
  #' @param exp The main exposure variable for modeling (quote)
  #' @param cov The covariates that would be controlled for in modeling (quote list)
  #' @param knots Number of knots for restricted cubic spline
  #' @return The negative log partial likelihood
  #' @export
  if(k<=2){knots <- 2}else{knots = k}
  nD <- nrow(dat)
  dat[,time_case:= get(time)]
  covF <- cov!="NULL" ### This is a flag for whether there is not an exposure-only model 
  if(covF){cov_par <- as.numeric(na.omit(par[knots+1:length(par)]))}
  if (knots==2){
    dat[,gxz := as.data.table(cbind(rep(1,nD),get(marker))%*%par[1:knots])]
  }else{
    dat[,gxz := as.data.table(cbind(rep(1,nD),rcspline.eval(get(marker),nk=knots,inclx=T))%*%par[1:knots])]
    }
  if(covF){
    part2 <- dat[get(case_ind)==1 & !is.na(gxz), sum(gxz*get(main) + as.matrix(.SD[,eval(parse(text=cov))])%*%cov_par)] 
  }else{
    part2 <- dat[get(case_ind)==1 & !is.na(gxz), sum(gxz*get(main))]}
  part1byp <- dat[,{
    datc <- caseinfo(Dat =.SD)
    caseindex <- which(datc[,get(case_ind)]==1 & !is.na(datc$gxz))
    if (length(caseindex)==0){logsum=0
      }else{
        N <- dim(datc)[1]
        if(covF){
          expsum <- sapply(caseindex, function(v) 
            sum(exp(datc$gxz[v]*datc[index[v]:N,get(main)]
                    +as.matrix(datc[index[v]:N,eval(parse(text=cov))])%*%cov_par))) 
        }else{expsum <- sapply(caseindex, function(v) sum(exp(datc$gxz[v]*datc[index[v]:N,get(main)])))}
    logsum <- sum(log(expsum))}
    .(logsum = logsum)},by=eval(parse(text=by_var))]
  part1 <- sum(part1byp$logsum)
  neg_loglike <- part1-part2
  return(as.numeric(neg_loglike))
}

derv_part1 <- function(v,DATA,main,N,cov,Cov_names,COV_PAR,K){
 covF <- unique(c(cov,Cov_names,COV_PAR))!="NULL"
 if (covF){
   denom <- exp(DATA$gxz[v]*DATA[index[v]:N,get(main)]
                +as.matrix(DATA[index[v]:N,eval(parse(text=cov))])%*%COV_PAR)
 }else{denom <- exp(DATA$gxz[v]*DATA[index[v]:N,get(main)])}
 ZIXL <- sapply(seq(K), function(m) sum(denom*as.numeric(DATA[v,paste("Z",m,sep=""),with=F])
                                            *DATA[index[v]:N,get(main)]))
 if (covF){
   ZIWL <- sapply(Cov_names, function(n) sum(denom*DATA[index[v]:N,n,with=F]))
   num_vec <- c(ZIXL,ZIWL)}else{num_vec <- c(ZIXL)}
 denom_sum <- sum(denom)
 num_denom <- num_vec/denom_sum
 if (covF){names(num_denom)<- c(sapply(seq(K), function(m) paste("Z",m,sep="")),Cov_names)}else{
  names(num_denom)<- sapply(seq(K), function(m) paste("Z",m,sep=""))}
 return(num_denom)
}

derv_spline <- function(par,dat,main,cov,case_ind,marker,time,k,by_var="list(interval)"){
  #' This is a function to compute the derivative in terms of all coefficients for the negative log likelihood
  #' With Brewlow method implemented for handling ties
  #' @param par The estimated coefficient values (numeric vector)
  #' @param dat The given dataset (data.table) 
  #' @param exp The main exposure variable for modeling (quote)
  #' @param cov The covariates that would be controlled for in modeling (quote list)
  #' @param knots Number of knots for restricted cubic spline
  #' @return The derivative vector
  #' @export
  if(k<=2){knots <- 2}else{knots = k}
  nD <- nrow(dat)
  dat[,time_case:= get(time)]
  covF <- cov!="NULL" ### This is a flag for whether there is not an exposure-only model 
  if(covF){cov_par <- as.numeric(na.omit(par[knots+1:length(par)]))}else{cov_par <- "NULL"}
  if (knots==2){
    dat[,gxz := as.data.table(cbind(rep(1,nD),get(marker))%*%par[1:knots])]
    dat[,sapply(seq(knots), function(m) paste("ZX",m,sep="")) 
        := as.data.table(cbind(rep(1,nD),get(marker))*get(main))]
    dat[,sapply(seq(knots), function(m) paste("Z",m,sep="")) 
        := as.data.table(cbind(rep(1,nD),get(marker)))]
  }else{
    dat[,gxz := as.data.table(cbind(rep(1,nD),rcspline.eval(get(marker),nk=knots,inclx=T))%*%par[1:knots])]
    dat[,sapply(seq(knots), function(m) paste("ZX",m,sep="")) 
        := as.data.table(cbind(rep(1,nD),rcspline.eval(get(marker),nk=knots,inclx=T))*get(main))]
    dat[,sapply(seq(knots), function(m) paste("Z",m,sep="")) 
        := as.data.table(cbind(rep(1,nD),rcspline.eval(get(marker),nk=knots,inclx=T)))]}
  if(covF){Cov_names <- strsplit(gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.","",
                                      x=cov,perl=T),split=",")[[1]]
  }else{Cov_names<- "NULL"}
  ### If we did not include covariates, objects cov_par & Cov_names all equal to NULL
  if(covF){
    part2VEC <- c(sapply(seq(knots), function(m) paste("ZX",m,sep="")),Cov_names)
    part1VEC <- c(sapply(seq(knots), function(m) paste("Z",m,sep="")),Cov_names)
  }else{
      part2VEC <- sapply(seq(knots), function(m) paste("ZX",m,sep=""))
      part1VEC <- sapply(seq(knots), function(m) paste("Z",m,sep=""))}
  part2 <- dat[get(case_ind)==1 & !is.na(gxz), sapply(.SD, sum, na.rm=TRUE), 
               .SDcols = names(dat) %in% part2VEC]
  part2 <- part2[order(factor(names(part2), levels = part2VEC))]
  part1 <- dat[,{
    datc <- caseinfo(.SD)
    caseindex <- which(datc[,get(case_ind)]==1 & !is.na(datc$gxz))
    if (length(caseindex)!=0){
      n <- nrow(datc)
      nom_denom <- rowSums(sapply(caseindex, function(m) derv_part1(m,DATA=datc,main=main,
                   N=n,cov=cov,Cov_names=Cov_names,COV_PAR = cov_par,K = knots))) 
    }else{
      nom_denom <-rep(0,length(part1VEC))
      if(covF){names(nom_denom) <- c(sapply(seq(knots), 
                                            function(m) paste("Z",m,sep="")),Cov_names)
      }else{names(nom_denom) <- sapply(seq(knots), function(m) paste("Z",m,sep=""))}}
    as.list(nom_denom)},by=eval(parse(text=by_var))]
  part1sum <- part1[, sapply(.SD, sum, na.rm=TRUE), .SDcols = names(part1) %in% part1VEC]
  part1sum <- part1sum[order(factor(names(part1sum), levels = part1VEC))]
  par_derv <- part1sum - part2
  return(par_derv)
}

uncertain_info <- function(opt){
  fisher_info<-ginv(opt$hessian) 
  prop_sigma<-sqrt(diag(fisher_info))
  prop_sigma<-diag(diag(prop_sigma))
  upper<-opt$par+1.96*prop_sigma
  lower<-opt$par-1.96*prop_sigma
  p_value <- pnorm(opt$par/prop_sigma,lower.tail = F)
  interval<-data.frame(value=opt$par,lower=lower,upper=upper,sig=p_value)
  return(interval)
}

hazard_plot <- function(Dat,optim_rlt,knots,FileName,both=F){
  if (both){
    marker_seq <- seq(20,100,0.5)
    PAR <- lapply(optim_rlt, function(v) v$value)
    hazard <- lapply(PAR, function(m) exp(as.data.table(cbind(rep(1,length(marker_seq)),
                                      rcspline.eval(marker_seq,nk=knots,inclx=T))%*%m[1:knots]))$V1)
    hazard <- unlist(hazard)
    ind <- c(rep("NHS",length(marker_seq)),rep("HPFS",length(marker_seq)))
    marker <-rep(marker_seq,2)
    df <- data.frame(hazard=hazard,ind=ind,marker=marker)
    ggplot() +
      geom_point(data=df, mapping=aes(x=marker, y=hazard,shape = factor(ind),colour=factor(ind)))+
      geom_hline(yintercept=1, linetype="dashed", color = "red", size=0.5)+
      xlab("Line-1 marker")+
      ylab("Corresponding Hazard")
    ggsave(FileName)
  }else{
    marker_range <- range(na.omit(Dat[,get(marker)]))
    marker_seq <- seq(marker_range[1],marker_range[2],0.5)
    PAR <- optim_rlt$value
    hazard <- exp(as.data.table(cbind(rep(1,length(marker_seq)),
                                      rcspline.eval(marker_seq,nk=knots,inclx=T))%*%PAR[1:knots]))
    qplot(marker_seq, hazard$V1, geom=c("point", "smooth"), 
          xlab="Line-1 marker", ylab="Corresponding Hazard")
    ggsave(FileName)
  }
}

datproc <- function(dat){
  ANA <- dat[,.(id,time_crc,dt_crc,lineavgr,
                alcohol,
                ageyrs,interval, # Stratification variables
                bmicont,bmicontmiss,
                mets,metsmiss,
                bmigpm,bmigpr,bmigp2,bmigp3, # Categorical physical activity
                ccafamhis,polyp,aspirin,mvitamin, # Binary covariates
                smkr,smkm, # Never, Ever, Missing smoking status
                actm,actr,act2,act3,act4, # Categorical physical activity measurement
                calcium,folate,# For quantiles
                rpmeats,sex)]
  
  ANA[,calq1:=ifelse(calcium<quantile(calcium)[2],1,0)][
    ,calq2:=ifelse(calcium<quantile(calcium)[3]&calcium>=quantile(calcium)[2],1,0)][
      ,calq3:=ifelse(calcium<quantile(calcium)[4]&calcium>=quantile(calcium)[3],1,0)][
        ,calq4:=ifelse(calcium>=quantile(calcium)[4],1,0)]
  
  ANA[,folq1:=ifelse(folate<quantile(folate)[2],1,0)][
    ,folq2:=ifelse(folate<quantile(folate)[3]&folate>=quantile(folate)[2],1,0)][
      ,folq3:=ifelse(folate<quantile(folate)[4]&folate>=quantile(folate)[3],1,0)][
        ,folq4:=ifelse(folate>=quantile(folate)[4],1,0)]
  
  ANA[,Eversmk:=ifelse(smkr==0&smkm==0,1,0)]
  
  ANA[,rmq1:=ifelse(rpmeats<quantile(rpmeats)[2],1,0)][
    ,rmq2:=ifelse(rpmeats<quantile(rpmeats)[3]&rpmeats>=quantile(rpmeats)[2],1,0)][
      ,rmq3:=ifelse(rpmeats<quantile(rpmeats)[4]&rpmeats>=quantile(rpmeats)[3],1,0)][
        ,rmq4:=ifelse(rpmeats>=quantile(rpmeats)[4],1,0)]
  
  ANA[,bmi:=bmicont*(1-bmicontmiss)]
  ANA[,pact:=mets*(1-metsmiss)]
  
  ANA <- ANA[,.(id,time_crc,dt_crc,lineavgr,
                alcohol,
                ageyrs,interval,sex, # Stratification variables
                bmi,bmicontmiss,
                pact,metsmiss,
                Neversmk = smkr,smkm,Eversmk, # Never, Ever, Missing smoking status
                ccafamhis,
                bmigpm,bmigpr,bmigp2,bmigp3, # Categorical physical activity
                polyp,aspirin,mvitamin, # Binary covariates
                actm,actr,act2,act3,act4, # Categorical physical activity measurement
                folq1,folq2,folq3,folq4,
                calq1,calq2,calq3,calq4,
                rmq1,rmq2,rmq3,rmq4)] 
  return(ANA)
}



