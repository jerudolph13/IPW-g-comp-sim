
###################################################################################################
#
# Purpose: Compare using simulation the performance of g-comp and IPTW to estimate same 
#          parameter in presence of:
#           -time-fixed exposure
#           -3 time-fixed confounders
#           -binary outcome
#
# Author: Jacqueline Rudolph
#
# Last Update: 29 apr 2020
#
##################################################################################################

setwd("/home/jackie/Documents/ipwgcomp/results")

packages <- c("survival", "nnet", "tidyverse", "data.table", "flexsurv", "parallel", 
              "doParallel", "geepack", "reda")
for (package in packages) {
  library(package, character.only=T)
}

expit <- function(x) {1/(1+exp(-x))}

ptm <- proc.time()


##################################################################################################

### Data generation
##
n <- 5000 # Number of subjects
K <- 1 # Number of causes of death

montecarlo <- 6000
nmethods <- 3 # Number of methods to estimate HR
nsim <- 200 # Number of simulations (increase this # after we are sure code works)
nboot <- 200 # Number of bootstrap resamples

# Prepare data set to hold simulation results

sim.res <- data.frame(
  method=c("Oracle", "IPW", "G-comp"),
  lhr=rep(NA, 3),
  sd=rep(NA, 3),
  coverage=rep(NA, 3),
  stringsAsFactors = FALSE
)

simloop <- function(s, nboot, montecarlo){
  
  cat("Now running simulation",s,'\n')
  
  sim <- s
  set.seed(sim)
  
simulation <- function (exposure) {
  
  ID <- c(1:n)
  L <- rbinom(n, 1, 0.3) 
  J <- rbinom(n, 1, 0.4)
  M <- rbinom(n, 1, 0.5) 
  
  if (is.null(exposure)){
    p_a <- 1/(1 + exp(-(-log(1/0.5 - 1) + 0.5*L - 0.5*0.3 + 0.6*J - 0.6*0.4 + 0.8*M - 0.8*0.5)))
    A <- rbinom(n, 1, p_a)
  }
  else {
    A <- exposure
  }
  
  p_t  <- (0.05)*exp(log(2)*A + 0.5*L + 0.6*J + 0.8*M) #In main analysis, lambda=0.01; for common outcome, lambda=0.05
  Tv <- rexp(n, p_t)
  Z <- as.numeric(Tv < 10)
  Tv <- ifelse(Tv > 10, 10, Tv)

DeathsK.df <- data.frame(ID, A, L, J, M, Tv, Z)
return(DeathsK.df)
}

#Oracle: get counterfactuals for every simulated individual by setting exposure inside simulation
oracle1 <- simulation(exposure=1)
oracle0 <- simulation(exposure=0)
oracle <- data.table(rbind(oracle0,oracle1))

cat('\n')
cat("Fitting oracle",'\n')
oracle.aft <- flexsurvreg(Surv(Tv, Z) ~ A, data=oracle, dist="exp")
sim.res$lhr[1] <- oracle.aft$coefficients[2]

#Now create simulation to be used in analysis steps
DeathsK.df <- simulation(exposure=NULL)


##################################################################################################
## Code for IPTW (Comes from Moodie supp appendix)
    ## Bootstrap resampling to get CI

# Set up data set to hold bootstrap results  
boot.res <- data.frame(
  method=c("IPW", "G-comp"),
  lhr=rep(NA, 2),
  boot=rep(NA, 2),
  stringsAsFactors=FALSE
)
  
bootrep <- function(r) {
  cat("Now running simulation",s,"bootstrap",r,'\n')
  
  boot.res$boot <- rep(r, 2)
  set.seed(1000 * sim + r)
  firstobs <- DeathsK.df
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
mod <- glm(A ~ L + J + M, family=binomial(link="logit"), data=boot)
logit <- predict(glm(A ~ L + J + M, family=binomial(link="logit"), data=boot))
denominator[boot$A == 1] <- expit(logit[boot$A == 1])
denominator[boot$A == 0] <- 1 - expit(logit[boot$A == 0])

#Numerator of weights
numerator <- rep(NA, nrow(boot))
logit <- predict(glm(A ~ 1, family=binomial, data=boot))
numerator[boot$A == 1] <- expit(logit[boot$A == 1])
numerator[boot$A == 0] <- 1 - expit(logit[boot$A == 0])
wt <- numerator / denominator 

#IP-weighted AFT Model
cat('\n')
cat("Fitting IPW model, boot",r,'\n')  
iptw.aft <- flexsurvreg(Surv(Tv, Z) ~ A, data=boot, dist="exp", weights=wt)
  boot.res$lhr[1] <- iptw.aft$coefficients[2]
  

##################################################################################################
## G-computation
    ## Time-ordering: L, J, M, A, Z
    ## Bootstrap resampling to get CI
    
    #Model exposure (if I ever want to generate NH)
    # cat('\n')
    # cat("Fitting g-comp A model, boot",r,'\n') 
    # mod.A <- glm(A ~ L + J + M, family=binomial(link="logit"), data=boot)
    
    #Model outcome
    cat('\n')
    cat("Fitting g-comp Y model, boot",r,'\n')
    mod.D <- flexsurvreg(Surv(Tv, Z) ~ A + L + J + M, data=boot, dist="exp")
    
    #Take a MC Sample
    #Select first obs for each person to obtain joint empirical distribution of baseline covariates
    set.seed(1000 * sim + 500 + r)
    MC0<-boot[,(names(boot) %in% c("L", "J", "M", "A", "rep"))]
    index<-sample(1:nrow(MC0), montecarlo, replace=T)
    MC<-MC0[index,]
    # head(MC)
    MC$id<-1:montecarlo
    
    #Predict follow-up based on g-formula using PGF function
    cat('\n')
    cat("Running pgf function, boot",r,'\n')    
    pgf<-function(mc_data, exposure){
      d <- mc_data
      id <- d$id
      rep <- d$rep
      
      L <- d$L
      J <- d$J
      M <- d$M
      
      if (is.null(exposure)) {
        A <- d$A
      } else {
        A <- exposure
      }
      
      p_t <- exp(coef(mod.D)[names(coef(mod.D))=="rate"])*
        exp(coef(mod.D)[names(coef(mod.D))=="A"]*A 
            + coef(mod.D)[names(coef(mod.D))=="L"]*L
            + coef(mod.D)[names(coef(mod.D))=="J"]*J
            + coef(mod.D)[names(coef(mod.D))=="M"]*M)
      Tv <- rexp(montecarlo, p_t)
      Z <- ifelse(Tv < 10, 1, 0)
      Tv <- ifelse(Tv > 10, 10, Tv)
      
      gdat <- data.table(id,A,L,J,M,Tv,Z,rep)
      return(gdat)
    }
    
    #cores <- detectCores()
    res0 <- pgf(mc_data=MC, exposure=0)
    res1 <- pgf(mc_data=MC, exposure=1)
      gcomp.dat <- data.table(rbind(res1, res0))

#Run outcome model
    #Pooled logistic
    gcomp.aft <- flexsurvreg(Surv(Tv, Z) ~ A, data=gcomp.dat, dist="exp")
      boot.res$lhr[2] <- gcomp.aft$coefficients[2]
    
    return(boot.res)
}

#cores <- detectCores()
all.boot <- lapply(0:nboot, function(tt) {bootrep(tt)})#, mc.cores=cores, mc.set.seed=FALSE)
all.boot <- do.call(rbind, all.boot)

#For point estimate, pull out results where boot=0
boot0 <- filter(all.boot, boot == 0)
sim.res$lhr[2:3] <- boot0$lhr
  

##################################################################################################
##Aggregate results

#Summarize over bootstraps
boot.summ <- all.boot %>% 
  group_by(method) %>% 
  summarize(b.avg.lhr = mean(lhr), b.sd.lhr = sd(lhr))
boot.summ <- boot.summ[order(boot.summ$method, decreasing=TRUE), ]
sim.res$coverage[2:3] <- (sim.res$lhr[2:3] - 1.96 * boot.summ$b.sd.lhr) <= log(2) & 
                          log(2) <= (sim.res$lhr[2:3] + 1.96 * boot.summ$b.sd.lhr)
sim.res$sd[2:3] <- boot.summ$b.sd.lhr
sim.res$sim <- s

return(sim.res)
}

cores <- detectCores()
all.res <- mclapply(1:nsim, function(x) {simloop(x,nboot,montecarlo)}, mc.cores=cores, mc.set.seed=FALSE)
#all.res <- lapply(1:nsim, function(x) {simloop(x,nboot,montecarlo)})
all.res <- do.call(rbind, all.res)

#Summarize over simulations
res.summ <- all.res %>%
  group_by(method) %>%
  summarize(avg.lhr=mean(lhr), sd.lhr=sd(lhr), avg.se=mean(sd), avg.coverage=mean(coverage))
res.summ$bias <- res.summ$avg.lhr - log(2)
res.summ$mse <- res.summ$sd.lhr^2 + res.summ$bias^2

write.table(all.res, file="timefix_3-conf_n-5000_common.txt", sep="\t")

proc.time() - ptm


