
###################################################################################################
#
# Purpose: Compare using simulation the performance of g-comp and IPTW to estimate same 
#          parameter in presence of:
#           -time-varying exposure
#           -time-varying confounder
#           -binary outcome
#
# Author: Jacqueline Rudolph (Credit to Young and Moodie for DGM)
#
# Last Update: 24 feb 2020
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
nsim <- 200 # Number of simulations (increase this # after we are sure code works)
nboot <- 200 # Number of bootstrap resamples

# Prepare data set to hold simulation results

sim.res <- data.frame(
  method=c("Oracle PL", "Oracle Cox", "IPW PL RV", "IPW Cox RV", "IPW PL Boot", "IPW Cox Boot", 
           "G-comp PL", "G-comp Cox"),
  lhr=rep(NA, 8),
  sd=rep(NA, 8),
  coverage=rep(NA, 8),
  stringsAsFactors = FALSE
)

simloop <- function(s, nboot, montecarlo, n=1000, K=1, N=10, nmethods=4){
  
  cat("Now running simulation",s,'\n')
  
  sim <- s
  set.seed(sim)
  
## This is the matrix of parameters of interest, possibly different
## at each interval

psi.mat <- matrix(0, nrow=K, ncol=N+1)

##Here are the effect sizes for the K=3 causes

psi.mat[1, ] <- log(2)

##Here the (untreated) all-cause rate is set to lambda=0.01, with
##lambda/K per cause; muK=lambda is used in the algorithm.

lambda <- 0.05
gamma.vec <- rep(log(lambda/K))
muK <- sum(exp(gamma.vec)) #So this bit is necessary to deal with that inequality in algorithm?
A<-L<-ID<-Y<-Z<-Tv<-Int<-ALast<-LLast<-LFirst <- numeric()
T0.vec<-T.vec<-Y.vec<-Z.vec <- rep(0, n)

##Here are the coefficients determining the
##mediation and treatment assignment mechanisms.

bevec <- c(log(3/7), 2, log(0.5), log(1.5)) #Used to generate time-varying confounder
alvec <- c(log(2/7), 0.5, 0.5, log(4)) #Used to generate exposure

##cval is used as in Young's algorithm to introduce the confounding

cval <- 30

##Begin the data-generation loop
  #Do this 3 times so that I get oracle and simulation
  #First, where I get counterfactual for all at A=1
  #Second, where I get counterfactual for all at A=0
  #Third, where I simulate A value

simulation <- function (exposure) {

for (i in 1:n) {
  ##Generate the counterfactual (untreated) survival time
  T0 <- rexp(1, lambda) #Generate T0 from an exponential dist with constant rate=lamba
  Ival <- as.numeric(T0 < cval)
  ##Begin the interval-by-interval simulation
  m <- 0
  mu.tot <- 0
  A.vec<-L.vec<-ALast.vec<-LLast.vec<-LFirst.vec <- rep(0, N+1)
  ##Implement Young's algorithm with multiple causes
  ##Generate the survival time, then the cause
  while (muK*T0 > mu.tot & m <= N) {
    if (m == 0) {
      ##First interval
      eta <- bevec[1] + bevec[2] * Ival + bevec[3] * 0 + bevec[4] * 0
      pval <- 1 / (1 + exp(-eta))
      L.vec[m+1] <- rbinom(1, 1, pval)
      eta <- alvec[1] + alvec[2] * L.vec[m+1] + alvec[3] * 0 + alvec[4] * 0
      pval <- 1 / (1 + exp(-eta))
      if (is.null(exposure)) {A.vec[m + 1] <- rbinom(1, 1, pval)}
      else {A.vec[m + 1] <- exposure}

      ALast.vec[m + 1] <- 0; LLast.vec[m+1] <- 0 #The last values for A and L (but why create if it's not used?)
      LFirst.vec <- rep(L.vec[m + 1], N + 1) #The L at time=1
    } else {
      ##Subsequent intervals
      eta <- bevec[1] + bevec[2] * Ival + bevec[3] * A.vec[m] + bevec[4] * L.vec[m]
      pval <- 1 / (1 + exp(-eta))
      L.vec[m + 1] <- rbinom(1, 1, pval) #L affected by values of L and A at last time point
      eta <- alvec[1] + alvec[2] * L.vec[m + 1] + alvec[3] * L.vec[m] + alvec[4] * A.vec[m]
      pval <- 1 / (1 + exp(-eta)) #A affected by L at this time point, last point and A at last time point
      if (is.null(exposure)) {A.vec[m + 1] <- rbinom(1, 1, pval)}
      else {A.vec[m + 1] <- exposure}
      ALast.vec[m + 1] <- A.vec[m]; LLast.vec[m + 1] <- L.vec[m]
    }
    
    muval <- sum(exp(gamma.vec + A.vec[m + 1] * psi.mat[ , m + 1]))
    
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
    Z.vec[i] <- sample(c(1:K), 1, prob = exp(gamma.vec + A.vec[m] * psi.mat[ ,m])) # I don't really get this step
  }
  
  ##Store the outcomes
  T0.vec[i] <- T0
  T.vec[i] <- Tval
  Y.vec[i] <- m - 1
  ID <- c(ID, rep(i,m)) #Individual
  Int <- c(Int, c(1:m)) #Time point
  A <- c(A, A.vec[1:m]) #Time-updated treatment
  L <- c(L, L.vec[1:m]) #Time-updated covariate
  ALast <- c(ALast, ALast.vec[1:m]) #Treatment at last t point
  LLast <- c(LLast, LLast.vec[1:m]) #Covariate at last t point
  LFirst <- c(LFirst, LFirst.vec[1:m]) #Baseline covariate value
  Z <- c(Z, rep(0,m - 1), Z.vec[i]) #Outcome: Z>0 indicates outcome of some type, determined by value)
  tv <- c(1:m); tv[m] <- Tval
  Tv <- c(Tv, tv) #If event occurs, exact time at which it occurred; o.w. equal to Int)
}

DeathsK.df <- data.frame(ID, Int, Tv, A, ALast, L, LLast, LFirst, Z)

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
## Code for IPTW (Adapted from Moodie supp appendix)
    ## Get robust SE

cat('\n')
cat("Fitting weight models, no boot",'\n')

#Denominator of weights
denominator <- rep(NA, nrow(DeathsK.df))
logit <- predict(glm(A ~ ALast + L + LLast + as.factor(Int), family=binomial(link="logit"), data=DeathsK.df))
denominator[DeathsK.df$A == 1] <- expit(logit[DeathsK.df$A == 1])
denominator[DeathsK.df$A == 0] <- 1 - expit(logit[DeathsK.df$A == 0])

#Numerator of weights
numerator <- rep(NA, nrow(DeathsK.df))
logit <- predict(glm(A ~ ALast + as.factor(Int), family=binomial, data=DeathsK.df))
numerator[DeathsK.df$A == 1] <- expit(logit[DeathsK.df$A == 1])
numerator[DeathsK.df$A == 0] <- 1 - expit(logit[DeathsK.df$A == 0])
wt <- unlist(tapply(numerator / denominator, DeathsK.df$ID, cumprod))

#Save weight information for inspection
wts <- data.frame(wt, Int=DeathsK.df$Int, ID=DeathsK.df$ID)
wt.summ <- wts %>% 
  group_by(Int) %>% 
  summarize(avgwt=mean(wt), maxwt=max(wt))
wt.summ <- t(wt.summ[ , 2:3])
if (s==1) {write.table(wt.summ, file="wt_noboot.txt", sep="\t")}
else {write.table(wt.summ, file="wt_noboot.txt", sep="\t", col.names=FALSE, append=TRUE)}

#IP-weighted Pooled Logistic Model
cat('\n')
cat("Fitting IPW PL model, no boot",'\n')
iptw.pl <- geeglm(Z ~ A + as.factor(Int), data=DeathsK.df, family=binomial(link="logit"), weights=wt, id=ID)
  sim.res$lhr[3] <- coef(iptw.pl)[2]
  sim.res$sd[3] <- summary(iptw.pl)$coefficients[2,2]
  sim.res$coverage[3] <- (sim.res$lhr[3] - 1.96 * sim.res$sd[3]) <= log(2) &
                          log(2) <= (sim.res$lhr[3] + 1.96 * sim.res$sd[3])
  
#IP-weighted Cox model
cat('\n')
cat("Fitting IPW Cox model, no boot",'\n') 
iptw.cox <- coxph(Surv(Int0, Tv, Z) ~ (A) + cluster(ID), data=DeathsK.df, weights=wt, 
                  ties="efron", timefix=FALSE)
  sim.res$lhr[4] <- iptw.cox$coef
  sim.res$sd[4] <- summary(iptw.cox)$coefficients[4]
  sim.res$coverage[4] <- (sim.res$lhr[4] - 1.96 * sim.res$sd[4]) <= log(2) &
                          log(2) <= (sim.res$lhr[4] + 1.96 * sim.res$sd[4])


##################################################################################################
## Code for IPTW (Comes from Moodie supp appendix)
    ## Bootstrap resampling to get CI

# Set up data set to hold bootstrap results  
boot.res <- data.frame(
  method=c("IPW PL Boot", "IPW Cox Boot", "G-comp PL", "G-comp Cox"),
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
logit <- predict(glm(A ~ ALast + L + LLast + as.factor(Int), family=binomial(link="logit"), data=boot))
denominator[boot$A == 1] <- expit(logit[boot$A == 1])
denominator[boot$A == 0] <- 1 - expit(logit[boot$A == 0])

#Numerator of weights
numerator <- rep(NA, nrow(boot))
logit <- predict(glm(A ~ ALast + as.factor(Int), family=binomial, data=boot))
numerator[boot$A == 1] <- expit(logit[boot$A == 1])
numerator[boot$A == 0] <- 1 - expit(logit[boot$A == 0])
wt <- unlist(tapply(numerator / denominator, boot$bid, cumprod))

#Save weight information for inspection
wts <- data.frame(wt, Int=boot$Int)
wt.summ <- wts %>% 
  group_by(Int) %>% 
  summarize(avgwt=mean(wt), maxwt=max(wt))
wt.summ <- t(wt.summ[ , 2:3])
if (s==1) {write.table(wt.summ, file="wt_boot.txt", sep="\t")}
else {write.table(wt.summ, file="wt_boot.txt", sep="\t", col.names=FALSE, append=TRUE)}

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
    ## Time-ordering: LLast, ALast, L, A, Z
    ## Bootstrap resampling to get CI
  
    #Model confounder
    cat('\n')
    cat("Fitting g-comp L model, boot",r,'\n')  
    mod.L <- glm(L ~ LLast + ALast + as.factor(Int), family=binomial(link="logit"), data=boot)
    
    #Model exposure (if I ever want to generate NH)
    # cat('\n')
    # cat("Fitting g-comp A model, boot",r,'\n') 
    # mod.A <- glm(A ~ ALast + L + LLast + as.factor(Int), family=binomial(link="logit"), data=boot)
    
    #Model outcome
    cat('\n')
    cat("Fitting g-comp Y model, boot",r,'\n')
    mod.D <- glm(Z ~ L + LLast + ALast + A  + as.factor(Int), family=binomial(link="logit"), data=boot)
    mod.D2<-flexsurvreg(Surv(Int0,Tv,Z)~A+L+LLast+ALast,data=boot,dist="exp")
    
   
#Take a MC Sample
    #Create monte carlo dataset
    #Select first obs for each person to obtain joint empirical distribution of baseline covariates
    set.seed(1000 * sim + 500 + r)
    MC0<-boot[boot$Int==1,(names(boot) %in% c("L","A","rep"))]
    index<-sample(1:nrow(MC0),montecarlo,replace=T)
    MC<-MC0[index,]
    MC$id<-1:montecarlo
    
  cat('\n')
  cat("Running pgf function, boot",r,'\n')        
    
    #Predict follow-up based on g-formula using PGF function
    pgf<-function(ii, mc_data, length, exposure=NULL){
      
      pFunc<-function(mod,ndat){as.numeric(predict(mod,newdata=ndat,type="response")>runif(1))}
      d <- mc_data
      d<-d[d$id==ii,]
      lngth <- length
      Lp <- Ap <- Yp <- Yp2 <- time <- numeric()
      time[1] <- j <- 1
      id <- d$id
      rep <- d$rep
      
      Lp[1] <- d$L
      
      # if (is.null(exposure)) {
      #   Ap[1] <- d$A
      # } else {
      #   Ap[1] <- exposure
      # }
      
      Ap[1] <- exposure
      
      dYp <- data.table(L=Lp[1],LLast=0,ALast=0,A=Ap[1])
      expSim<-function(dat){
        newD<-dat
        desX<-newD[,c("A","L","LLast","ALast")]
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
      
      dYp2 <- data.table(L=Lp[1],LLast=0,ALast=0,A=Ap[1],Int=factor(j))
      Yp2[1] <- pFunc(mod.D,dYp2)
      
      for (j in 2:lngth) {
        if (Yp[j-1]==0){
          ALast=Ap[j-1];LLast=Lp[j-1]
          
          dLp <- data.table(LLast,ALast,Int=factor(j))
          Lp[j] <- pFunc(mod.L,dLp)
          
          # if (is.null(exposure)) {
          #   dAp <- data.table(L=Lp[j],LLast,ALast,Int=factor(j))
          #   Ap[j] <- pFunc(mod.A,dAp)
          # } else { 
          #   Ap[j] <- exposure
          # }

          Ap[j] <- exposure
                    
          dYp <- data.table(A=Ap[j],L=Lp[j],LLast,ALast)
          t0<-expSim(dYp)
          if(t0<=0.001){ #if the time interval is short enough to cause an error, place event at end of previous time point
            Yp[j-1] <- 1
            Lp <- Lp[1:(j-1)]
            Ap <- Ap[1:(j-1)]
            id <- id[1:(j-1)]
            time <- time[1:(j-1)]
            jj <- jj[1:(j-1)]
            rep <- rep[1:(j-1)]
            break
          } else{
          if(t0<=1){
            Yp[j]<-1
            jj<-(j-1)+t0
          } else{
            Yp[j]<-0
            jj<-j
          }}
          
          dYp2 <- data.table(L=Lp[j],LLast,ALast,A=Ap[j],Int=factor(j))
          Yp2[j]<-pFunc(mod.D,dYp2)
          
        } else {
          break
        }
        time[j] <- j
      }
      #print(ii)
      gdat <- data.table(id,time,jj,Ap,Lp,Yp,Yp2,rep)
      gdat$last<-as.numeric(gdat$Yp!=0 | gdat$time==lngth)
      return(gdat)
    }
    
    #cores <- detectCores()
    #res0 <- mclapply(1:montecarlo,function(x) {pgf(x,mc_data=MC,length=N,exposure=0)}, mc.cores=cores, mc.set.seed=FALSE)
    res0<-lapply(1:montecarlo,function(x) {pgf(x,mc_data=MC,length=N,exposure=0)})
    res0 <- do.call(rbind,res0)

    #res1 <- mclapply(1:montecarlo,function(x) {pgf(x,mc_data=MC,length=N,exposure=1)}, mc.cores=cores, mc.set.seed=FALSE)
    res1<-lapply(1:montecarlo,function(x) {pgf(x,mc_data=MC,length=N,exposure=1)})
    res1 <- do.call(rbind,res1)
    
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

cores <- detectCores()
all.boot <- mclapply(0:nboot, function(tt) {bootrep(tt)}, mc.cores=cores, mc.set.seed=FALSE)
#all.boot <- lapply(0:nboot, function(tt) {bootrep(tt)})
all.boot <- do.call(rbind, all.boot)

#For point estimate, pull out results where boot=0
boot0 <- filter(all.boot, boot == 0)
sim.res$lhr[5:8] <- boot0$lhr
  

##################################################################################################
##Aggregate results

#Summarize over bootstraps
boot.summ <- all.boot %>% 
  group_by(method) %>% 
  summarize(b.avg.lhr = mean(lhr), b.sd.lhr = sd(lhr))
boot.summ <- boot.summ[order(boot.summ$method, decreasing=TRUE), ]
sim.res$coverage[5:8] <- (sim.res$lhr[5:8] - 1.96 * boot.summ$b.sd.lhr) <= log(2) & 
                          log(2) <= (sim.res$lhr[5:8] + 1.96 * boot.summ$b.sd.lhr)
sim.res$sd[5:8] <- boot.summ$b.sd.lhr
sim.res$sim <- s

return(sim.res)
}

#cores <- detectCores()
#all.res <- mclapply(51:nsim, function(x) {simloop(x,nboot,montecarlo)}, mc.cores=cores, mc.set.seed=FALSE)
#all.res <- lapply(1:nsim, function(x) {simloop(x,nboot,montecarlo)})
all.res <- lapply(193, function(x) {simloop(x,nboot,montecarlo)})
all.res <- do.call(rbind, all.res)

#Summarize over simulations
res.summ <- all.res %>%
  group_by(method) %>%
  summarize(avg.lhr = mean(lhr), sd.lhr = sd(lhr), avg.coverage = mean(coverage))
res.summ$bias <- res.summ$avg.lhr - log(2)
res.summ$mse <- res.summ$sd.lhr^2 + res.summ$bias^2

proc.time() - ptm
 
write.table(all.res, file="ipwgcomp_common-outcome_1-conf_193.txt", sep="\t", row.names=FALSE)
