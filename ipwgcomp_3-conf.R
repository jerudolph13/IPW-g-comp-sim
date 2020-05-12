
###################################################################################################
#
# Purpose: Compare using simulation the performance of g-comp and IPTW to estimate same 
#          parameter in presence of:
#           -time-varying exposure
#           -3 time-varying confounders
#           -binary outcome
#
# Author: Jacqueline Rudolph (Credit to Young and Moodie for DGM)
#
# Last Update: 5 feb 2020
#
##################################################################################################

setwd("/home/jer/Documents/gcompipwsim")

packages <- c("survival", "nnet", "tidyverse", "data.table", "flexsurv", "parallel", 
              "doParallel", "geepack")
for (package in packages) {
  library(package, character.only=T)
}

expit <- function(x) {1/(1+exp(-x))}

ptm <- proc.time()


##################################################################################################

## This code generates data from a structural nested model 
## compatible with a marginal structural model.
## It is based on Jessica Young's algorithm, published here:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3635680/
#
## This code, which extend's Young's algorithm to multiple outcomes
## was written by Erica Moodie, published here:
# https://www.ncbi.nlm.nih.gov/pubmed/24272681
##
### Data generation
##
n <- 1000 # Number of subjects
N <- 10 #number of intervals per subject
K <- 1 # Number of causes of death

montecarlo <- 6000
nmethods <- 4 # Number of methods to estimate HR
nsim <- 50 # Number of simulations (increase this # after we are sure code works)
nboot <- 200 # Number of bootstrap resamples

# Prepare data set to hold simulation results
sim.res <- data.frame(
  method=c("Oracle PL", "Oracle Cox", "IPW PL", "IPW Cox", "G-comp PL", "G-comp Cox"),
  lhr=rep(NA, 6),
  sd=rep(NA, 6),
  coverage=rep(NA, 6),
  error=rep(NA, 6),
  stringsAsFactors = FALSE
)

simloop <- function(s, nboot, montecarlo, n=1000, K=1, N=10, nmethods=4){
  
  cat("Now running simulation",s,'\n')
  s <- 1
  sim <- s
  set.seed(sim)
  
  ## This is the matrix of parameters of interest, possibly different
  ## at each interval
  psi.mat <- matrix(0, nrow=K, ncol=N+1)
  
  ##Here are the effect sizes for the K=1 cause
  psi.mat[1, ] <- log(2)
  
  ##Here the (untreated) all-cause rate is set to lambda=0.01, with
  ##lambda/K per cause; muK=lambda is used in the algorithm.
  lambda <- 0.01
  gamma.vec <- rep(log(lambda/K))
  muK <- sum(exp(gamma.vec)) #So this bit is necessary to deal with that inequality in algorithm?
  A<-J<-M<-L<-ID<-Y<-Z<-Tv<-Int<-ALast<-LLast<-LFirst<-JLast<-JFirst<-MLast<-MFirst <- numeric()
  T0.vec<-T.vec<-Y.vec<-Z.vec <- rep(0, n)
  
  ##Here are the coefficients determining the
  ##mediation and treatment assignment mechanisms.
  bevec <- c(log(3/7), 2, log(0.5), log(1.5)) #Used to generate time-varying confounder L
  bevec2 <- c(log(4/7), 3, log(0.6), log(1.4)) #Used to generate time-varying confounder J
  bevec3 <- c(log(5/7), 4, log(0.7), log(1.3)) #Used to generate time-varying confounder M
  alvec <- c(log(2/7), 0.5, 0.5, log(4), 0.6, 0.6, 0.8, 0.8) #Used to generate exposure (Intercept, L, LLast, ALast, J, JLast, M, MLast)
  
  ##cval is used as in Young's algorithm to introduce the confounding
  cval <- 30
  
  ##Begin the data-generation loop
  simulation <- function (exposure) {
    
    for (i in 1:n) {
      ##Generate the counterfactual (untreated) survival time
      T0 <- rexp(1, lambda) #Generate T0 from an exponential dist with constant rate=lamba
      Ival <- as.numeric(T0 < cval)
      ##Begin the interval-by-interval simulation
      m <- 0
      mu.tot <- 0
      A.vec<-J.vec<-M.vec<-L.vec<-ALast.vec<-JLast.vec<-MLast.vec<-LLast.vec<-JFirst.vec<-MFirst.vec<-LFirst.vec <- rep(0, N+1)
      ##Implement Young's algorithm with multiple causes
      ##Generate the survival time, then the cause
      while (muK*T0 > mu.tot & m <= N) {
        if (m == 0) {
          ##First interval
          eta <- bevec[1] + bevec[2]*Ival + bevec[3]*0 + bevec[4]*0
          pval <- 1 / (1 + exp(-eta))
          L.vec[m+1] <- rbinom(1, 1, pval)
          eta <- bevec2[1] + bevec2[2]*Ival + bevec2[3]*0 + bevec2[4]*0
          pval <- 1 / (1 + exp(-eta))
          J.vec[m+1] <- rbinom(1, 1, pval)
          eta <- bevec3[1] + bevec3[2]*Ival + bevec3[3]*0 + bevec3[4]*0
          pval <- 1 / (1 + exp(-eta))
          M.vec[m+1] <- rbinom(1, 1, pval)
          
          eta <- alvec[1] + alvec[2]*L.vec[m+1] + alvec[3]*0 + alvec[4]*0 + alvec[5]*J.vec[m+1] + alvec[6]*0 + alvec[7]*M.vec[m+1] + alvec[8]*0
          pval <- 1 / (1 + exp(-eta))
          if (is.null(exposure)) {A.vec[m + 1] <- rbinom(1, 1, pval)}
          else {A.vec[m+1] <- exposure}
          
          ALast.vec[m+1] <- 0; LLast.vec[m+1] <- 0; JLast.vec[m+1] <- 0; MLast.vec[m+1] <- 0
          LFirst.vec <- rep(L.vec[m+1], N + 1)
          JFirst.vec <- rep(J.vec[m+1], N + 1)
          MFirst.vec <- rep(M.vec[m+1], N + 1)
          
        } else {
          ##Subsequent intervals
          eta <- bevec[1] + bevec[2]*Ival + bevec[3]*A.vec[m] + bevec[4]*L.vec[m]
          pval <- 1 / (1 + exp(-eta))
          L.vec[m+1] <- rbinom(1, 1, pval) 
          eta <- bevec2[1] + bevec2[2]*Ival + bevec2[3]*A.vec[m] + bevec2[4]*J.vec[m]
          pval <- 1 / (1 + exp(-eta))
          J.vec[m+1] <- rbinom(1, 1, pval) 
          eta <- bevec3[1] + bevec3[2]*Ival + bevec3[3]*A.vec[m] + bevec3[4]*M.vec[m]
          pval <- 1 / (1 + exp(-eta))
          M.vec[m+1] <- rbinom(1, 1, pval) 
          
          eta <- alvec[1] + alvec[2]*L.vec[m+1] + alvec[3]*L.vec[m] + alvec[4]*A.vec[m] +alvec[5]*J.vec[m+1] + alvec[6]*J.vec[m] 
          + alvec[7]*M.vec[m+1] + alvec[8]*M.vec[m]
          pval <- 1 / (1 + exp(-eta)) #A affected by L at this time point, last point and A at last time point
          if (is.null(exposure)) {A.vec[m+1] <- rbinom(1, 1, pval)}
          else {A.vec[m+1] <- exposure}
          ALast.vec[m+1] <- A.vec[m]; LLast.vec[m+1] <- L.vec[m]; JLast.vec[m+1] <- J.vec[m]; MLast.vec[m+1] <- M.vec[m]
        }
        
        muval <- sum(exp(gamma.vec + A.vec[m+1]*psi.mat[ , m+1]))
        
        ##Tval is computed for each interval, but is overwritten
        ##until the final interval
        Tval <- m + (muK * T0 - mu.tot) / muval
        mu.tot <- mu.tot + muval
        m <- m + 1
      }
      
      ##After exiting the loop, the survival time has been generated as Tval
      ##Now need to generate the failure type.
      if (m > N) {
        ##In the case of censoring at tenth interval, no failure.
        Tval <- m - 1
        Z.vec[i] <- 0
      } else {
        ##In the case of failure, use the ratio hazards to define the
        ##relevant multinomial distribution on the K causes.
        Z.vec[i] <- sample(c(1:K), 1, prob = exp(gamma.vec + A.vec[m]*psi.mat[ ,m])) # I don't really get this step
      }
      
      ##Store the outcomes
      T0.vec[i] <- T0
      T.vec[i] <- Tval
      Y.vec[i] <- m - 1
      ID <- c(ID, rep(i,m)) #Individual
      Int <- c(Int, c(1:m)) #Time point
      A <- c(A, A.vec[1:m]) #Time-updated treatment
      L <- c(L, L.vec[1:m]) #Time-updated covariate L
      J <- c(J, J.vec[1:m]) #Time-updated covariate J
      M <- c(M, M.vec[1:m]) #Time-updated covariate M
      ALast <- c(ALast, ALast.vec[1:m]) #Treatment at last t point
      LLast <- c(LLast, LLast.vec[1:m]) #Covariate L at last t point
      JLast <- c(JLast, JLast.vec[1:m]) #Covariate J at last t point
      MLast <- c(MLast, MLast.vec[1:m]) #Covariate M at last t point
      LFirst <- c(LFirst, LFirst.vec[1:m]) #Baseline covariate L value
      JFirst <- c(JFirst, JFirst.vec[1:m]) #Baseline covariate J value
      MFirst <- c(MFirst, MFirst.vec[1:m]) #Baseline covariate M value
      Z <- c(Z, rep(0,m - 1), Z.vec[i]) #Outcome: Z>0 indicates outcome of some type, determined by value)
      tv <- c(1:m); tv[m] <- Tval
      Tv <- c(Tv, tv) #If event occurs, exact time at which it occurred; o.w. equal to Int)
    }
    
    DeathsK.df <- data.frame(ID, Int, Tv, A, ALast, L, LLast, LFirst, J, JLast, JFirst, M, MLast, MFirst, Z)
    
    ##Trim off the intervals beyond the Nth (loop goes one too far)
    DeathsK.df <- DeathsK.df[DeathsK.df$Int <= N, ]
    DeathsK.df$Int0 <- DeathsK.df$Int - 1
    
    return(DeathsK.df)
  }
  

#Oracle: get counterfactuals for every simulated individual by setting exposure inside simulation
oracle1 <- simulation(exposure=1)
oracle0 <- simulation(exposure=0)
oracle <- data.table(rbind(oracle0,oracle1))

cat('\n')
cat("Fitting oracle",'\n')
oracle.pl<- glm(Z ~ A + as.factor(Int), data=oracle, family=binomial(link="logit"))
    sim.res$lhr[1] <- coef(oracle.pl)[2]
oracle.cox <- coxph(Surv(Int0, Tv, Z)  ~ A, data=oracle, ties="efron", timefix=FALSE)
    sim.res$lhr[2] <- oracle.cox$coef

#Now create simulation to be used in analysis steps
DeathsK.df <- simulation(exposure=NULL)


##################################################################################################
## IPTW with bootstrap resampling to get CI

# Set up data set to hold bootstrap results  
boot.res <- data.frame(
  method=c("IPW PL", "IPW Cox", "G-comp PL", "G-comp Cox"),
  lhr=rep(NA, 4),
  boot=rep(NA, 4),
  stringsAsFactors=FALSE
)
  
bootrep <- function(r) {
  cat("Now running simulation",s,"bootstrap",r,'\n')
  
  boot.res$boot <- rep(r, 4)
  set.seed(1000 * sim + r)
  firstobs <- DeathsK.df[DeathsK.df$Int == 1, ]
  samp <- table(firstobs[sample(1:nrow(firstobs),nrow(firstobs),replace=T), (names(DeathsK.df) == "ID")])

# The step below pulls in the simulated data for boot=0; otherwise grabs all records for the resampled observations
  boot <- NULL
  if(r==0){
    boot<-DeathsK.df
    boot$bid<-DeathsK.df$ID
  } else{
    for(zzz in 1:max(samp)){ # this counter goes from zero to select empirical data (no resample)
      cc <- DeathsK.df[DeathsK.df$ID %in% names(samp[samp %in% c(zzz:max(samp))]),]
      cc$bid<-paste0(cc$ID,zzz)
      boot <- rbind(boot, cc)
    }}
  
cat('\n')
cat("Fitting weight models, boot",r,'\n')    
#Denominator of weights
denominator <- rep(NA, nrow(boot))
logit <- predict(glm(A ~ ALast + L + LLast + J + JLast + M + MLast + as.factor(Int), family=binomial(link="logit"), data=boot))
denominator[boot$A == 1] <- expit(logit[boot$A == 1])
denominator[boot$A == 0] <- 1 - expit(logit[boot$A == 0])

#Numerator of weights
numerator <- rep(NA, nrow(boot))
logit <- predict(glm(A ~ ALast + as.factor(Int), family=binomial, data=boot))
numerator[boot$A == 1] <- expit(logit[boot$A == 1])
numerator[boot$A == 0] <- 1 - expit(logit[boot$A == 0])
wt <- unlist(tapply(numerator / denominator, boot$bid, cumprod))

#IP-weighted Pooled Logistic Model
cat('\n')
cat("Fitting IPW PL model, boot",r,'\n')  
iptw.pl <- glm(Z ~ A + as.factor(Int), data=boot, family=binomial(link="logit"), weights=wt)
  boot.res$lhr[1] <- coef(iptw.pl)[2]

#IP-weighted Cox model
cat('\n')
cat("Fitting IPW Cox model, boot",r,'\n')    
iptw.cox <- coxph(Surv(Int0, Tv, Z)  ~ A, data=boot, weights=wt, ties="efron", timefix=FALSE)
  boot.res$lhr[2] <- iptw.cox$coef
  

##################################################################################################
## G-computation
    ## Time-ordering: LLast, JLast, MLast, ALast, L, J, M, A, Z
    ## Bootstrap resampling to get CI
  
    #Model confounders
    cat('\n')
    cat("Fitting g-comp L model, boot",r,'\n')  
    mod.L <- glm(L ~ LLast + ALast + as.factor(Int), family=binomial(link="logit"), data=boot)
  
    cat('\n')
    cat("Fitting g-comp J model, boot",r,'\n')  
    mod.J <- glm(J ~ JLast + ALast + as.factor(Int), family=binomial(link="logit"), data=boot)    
  
    cat('\n')
    cat("Fitting g-comp L model, boot",r,'\n')  
    mod.M <- glm(M ~ MLast + ALast + as.factor(Int), family=binomial(link="logit"), data=boot)
    
    #Model exposure (if I ever want to generate NH)
    # cat('\n')
    # cat("Fitting g-comp A model, boot",r,'\n') 
    # mod.A <- glm(A ~ ALast + L + LLast + as.factor(Int), family=binomial(link="logit"), data=boot)
    
    #Model outcome
    cat('\n')
    cat("Fitting g-comp Y model, boot",r,'\n')
    mod.D <- glm(Z ~ A + ALast + L + LLast + J + JLast + M + MLast + as.factor(Int), family=binomial(link="logit"), data=boot)
    mod.D2<-flexsurvreg(Surv(Int0,Tv,Z) ~ A + ALast + L + LLast + J + JLast + M + MLast, data=boot, dist="exp")
    
   
#Take a MC Sample
    #Create monte carlo dataset
    #Select first obs for each person to obtain joint empirical distribution of baseline covariates
    set.seed(1000 * sim + 500 + r)
    MC0<-boot[boot$Int==1,(names(boot) %in% c("L", "J", "M", "A","rep"))]
    index<-sample(1:nrow(MC0),montecarlo,replace=T)
    MC<-MC0[index,]
    # head(MC)
    MC$id<-1:montecarlo
    
  cat('\n')
  cat("Running pgf function, boot",r,'\n')        
    
  #Predict follow-up based on g-formula using PGF function
  pgf<-function(ii, mc_data, length, exposure=NULL){
    
    pFunc<-function(mod,ndat){as.numeric(predict(mod,newdata=ndat,type="response")>runif(1))}
    
    d <- mc_data
    d<-d[d$id==ii,]
    lngth <- length
    Lp <- Jp <- Mp <- Ap <- Yp <- time <- numeric()
    time[1] <- j <- 1
    id <- d$id
    rep <- d$rep
    
    Lp[1] <- d$L; Jp[1] <- d$J; Mp[1] <- d$M
    if (is.null(exposure)) {
      Ap[1] <- d$A
    } else{
      Ap[1] <- exposure
    }
    
    dYp <- data.table(A=Ap[1], ALast=0, L=Lp[1], LLast=0, J=Jp[1], JLast=0, M=Mp[1], MLast=0)
    expSim<-function(dat){
      newD<-dat
      desX<-newD[,c("A", "ALast", "L","LLast", "J", "JLast", "M", "MLast")]
      y0<-rexp(1,exp(coef(mod.D2)[names(coef(mod.D2))=="rate"])*
                 exp(coef(mod.D2)[!names(coef(mod.D2))=="rate"]%*%t(desX)))
      return(y0)
    }
    t0<-expSim(dYp)
    if(t0<=1){
      Yp[1]<-1
      jj<-t0
    } else{
      Yp[1]<-0
      jj<-1
    }
    
    
    for (j in 2:lngth) {
      if (Yp[j-1]==0){
        ALast=Ap[j-1]; LLast=Lp[j-1]; JLast=Jp[j-1]; MLast=Mp[j-1]
        
        dLp <- data.table(LLast, ALast, Int=factor(j))
        Lp[j] <- pFunc(mod.L, dLp)
        dJp <- data.table(JLast, ALast, Int=factor(j))
        Jp[j] <- pFunc(mod.J, dJp)          
        dMp <- data.table(MLast, ALast, Int=factor(j))
        Mp[j] <- pFunc(mod.M, dMp)
        if (is.null(exposure)) {
          dAp <- data.table(L=Lp[j], LLast, J=Jp[j], JLast, M=Mp[j], MLast, ALast, Int=factor(j))
          Ap[j] <- pFunc(mod.A, dAp)
        } else{
          Ap[j] <- exposure
        }
        
        dYp <- data.table(A=Ap[j], ALast, L=Lp[j], LLast, J=Jp[j], JLast, M=Mp[j], MLast)
        t0<-expSim(dYp)
        if(t0<=0.001){
          Yp[j-1] <- 1 #If the time interval is short enough to cause an error, place event at end of previous time point
          Lp <- Lp[1:(j-1)]
          Jp <- Jp[1:(j-1)]
          Mp <- Mp[1:(j-1)]
          Ap <- Ap[1:(j-1)]
          id <- id[1:(j-1)]
          time <- time[1:(j-1)]
          jj <- jj[1:(j-1)]
          rep <- rep[1:(j-1)]
          break
        } else{
          if(t0>0.001 & t0<=1){
            Yp[j]<-1
            jj<-(j-1)+t0
          } else{ 
            Yp[j]<-0
            jj<-j
          }}
        
      } else {
        break
      }
      time[j] <- j
    }
    #print(ii)
    gdat <- data.table(id,time,jj,Ap,Lp,Jp,Mp,Yp,rep)
    gdat$last<-as.numeric(gdat$Yp!=0 | gdat$time==lngth)
    return(gdat)
  }
  
  cores <- detectCores()
  res0<-mclapply(1:montecarlo,function(x) {pgf(x,mc_data=MC,length=N,exposure=0)}, mc.cores=cores, mc.set.seed=FALSE)
  #res0<-lapply(1:montecarlo,function(x) {pgf(x,mc_data=MC,length=N,exposure=0)})
    res0<-do.call(rbind,res0)
    
  res1<-mclapply(1:montecarlo,function(x) {pgf(x,mc_data=MC,length=N,exposure=1)}, mc.cores=cores, mc.set.seed=FALSE)
  #res1<-lapply(1:montecarlo,function(x) {pgf(x,mc_data=MC,length=N,exposure=1)})
    res1<-do.call(rbind,res1)
    
    gcomp.dat <- data.table(rbind(res1, res0))
    gcomp.dat$Int0 <- gcomp.dat$time - 1
    gcomp.dat$Int <- ifelse(gcomp.dat$last==1, gcomp.dat$jj, gcomp.dat$time)
    
#Run outcome models
    #Pooled logistic
    gcomp.pl <- glm(Yp ~ Ap + as.factor(time), family=binomial(link="logit"), data=gcomp.dat)
      boot.res$lhr[3] <- coef(gcomp.pl)[2]
    
    #Cox model    
    gcomp.cox <- coxph(Surv(Int0, Int, Yp) ~ Ap, data=gcomp.dat, ties="efron", timefix=FALSE)
      boot.res$lhr[4] <- gcomp.cox$coef
    
    return(boot.res)
}

#cores <- detectCores()
#all.boot <- mclapply(0:nboot, function(tt) {bootrep(tt)}, mc.cores=cores, mc.set.seed=FALSE)
all.boot <- lapply(0:nboot, function(tt) {bootrep(tt)})
all.boot <- do.call(rbind, all.boot)

#For point estimate, pull out results where boot=0
boot0 <- filter(all.boot, boot == 0)
sim.res$lhr[3:6] <- boot0$lhr
sim.res$error <- abs(sim.res$lhr - log(2))
  

##################################################################################################
##Aggregate results

#Summarize over bootstraps
boot.summ <- all.boot %>% 
  group_by(method) %>% 
  summarize(b.avg.lhr = mean(lhr), b.sd.lhr = sd(lhr))
boot.summ <- boot.summ[order(boot.summ$method, decreasing=TRUE), ]
sim.res$coverage[3:6] <- (sim.res$lhr[3:6] - 1.96 * boot.summ$b.sd.lhr) <= log(2) & 
                          log(2) <= (sim.res$lhr[3:6] + 1.96 * boot.summ$b.sd.lhr)
sim.res$sd[3:6] <- boot.summ$b.sd.lhr
sim.res$sim <- s

return(sim.res)
}

#cores <- detectCores()
#all.res <- mclapply(1:nsim, function(x) {simloop(x,nboot,montecarlo)}, mc.cores=cores, mc.set.seed=FALSE)
all.res <- lapply(1:nsim, function(x) {simloop(x,nboot,montecarlo)})
all.res <- do.call(rbind, all.res)

#Summarize over simulations
res.summ <- all.res %>%
  group_by(method) %>%
  summarize(avg.lhr = mean(lhr), sd.lhr = sd(lhr), avg.se = mean(sd), avg.coverage = mean(coverage), avg.error = mean(error))
res.summ$bias <- res.summ$avg.lhr - log(2)
res.summ$mse <- res.summ$sd.lhr^2 + res.summ$bias^2

proc.time() - ptm

write.table(all.res, file="ipwgcomp_3-conf_50-results.txt", sep="\t", row.names=FALSE)
