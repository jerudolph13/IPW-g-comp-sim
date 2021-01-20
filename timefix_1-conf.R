
###################################################################################################
#
# Purpose: Compare using simulation the performance of g-comp and IPTW to estimate same 
#          parameter in presence of:
#           -time-fixed exposure
#           -1 time-fixed confounder
#           -binary outcome
#
# Author: Jacqueline Rudolph and Ashley Naimi
#
# Last Update: 08 sept 2020
#
##################################################################################################

# Read in packages
#lib <- "~/R/x86_64-pc-linux-gnu-library/4.0"
packages <- c("tidyverse", "data.table", "survival", "parallel")
for (package in packages) {
  library(package, character.only=T)#, lib.loc=lib)
}

# Pull in command line arguments
#args <- commandArgs(trailingOnly=TRUE)

# Define parameters and functions
nsim <- 1000				#Number of simulations
nboot <- 200				#Number of bootstrap resamples
n <- 1000 #as.numeric(args[1])		#Sample size
lambda <- 0.01 #as.numeric(args[2])		#Baseline rate
montecarlo <- 2000 #as.numeric(args[3])	#Monte Carlo resample size (0 implies no MC)

expit <- function(x) {1/(1+exp(-x))}

# Prepare data set to hold simulation results
sim.res <- data.frame(
  method=c("Oracle", "IPW Robust", "IPW Boot", "G-computation"),
  lhr=rep(NA, 4),
  se=rep(NA, 4),
  coverage=rep(NA, 4),
  sim=rep(NA, 4),
  stringsAsFactors = FALSE
)


##################################################################################################
#Get the true, marginal log HR

truth_func <- function(simN){
  ID <- c(1:simN)
  L <- rbinom(simN, 1, 0.3) 

  p_t0  <- (lambda)*exp(log(2)*0 + 0.5*L) #In main analysis, lambda=0.01; for common outcome, lambda=0.05
  Tv0 <- rexp(simN, p_t0)
  Z0 <- as.numeric(Tv0 < 10)
  Tv0 <- ifelse(Tv0 > 10, 10, Tv0)
  
  p_t1  <- (lambda)*exp(log(2)*1 + 0.5*L) #In main analysis, lambda=0.01; for common outcome, lambda=0.05
  Tv1 <- rexp(simN, p_t1)
  Z1 <- as.numeric(Tv1 < 10)
  Tv1 <- ifelse(Tv1 > 10, 10, Tv1)
  
  DeathsK.df0 <- data.frame(ID, A=0, L, Tv=Tv0, Z=Z0)
  DeathsK.df1 <- data.frame(ID, A=1, L, Tv=Tv1, Z=Z1)
  
  #ORACLE
  truth_dat <- data.table(rbind(DeathsK.df0,DeathsK.df1))
  
  truth.aft <- survreg(Surv(Tv, Z) ~ A, data=truth_dat, dist="exp")
  truth <- -truth.aft$coefficients[2]
  
  return(truth)
}

set.seed(123)
true <- truth_func(1e6)


##################################################################################################
# Data generation

sim_func <- function(iter){
  set.seed(iter)
  sim.res$sim <- rep(iter, 4)

  ID <- c(1:n)
  L <- rbinom(n, 1, 0.3) 

  p_a <- expit(-log(1/0.5 - 1) + 0.5*L - 0.5*0.3)
  A <- rbinom(n, 1, p_a)
  
  p_t  <- (lambda)*exp(log(2)*A + 0.5*L) #In main analysis, lambda=0.01; for common outcome, lambda=0.05
  Tv <- rexp(n, p_t)
  Z <- as.numeric(Tv < 10)
  Tv <- ifelse(Tv > 10, 10, Tv)
  
  p_t0  <- (lambda)*exp(log(2)*0 + 0.5*L) #In main analysis, lambda=0.01; for common outcome, lambda=0.05
  Tv0 <- rexp(n, p_t0)
  Z0 <- as.numeric(Tv0 < 10)
  Tv0 <- ifelse(Tv0 > 10, 10, Tv0)
  
  p_t1  <- (lambda)*exp(log(2)*1 + 0.5*L) #In main analysis, lambda=0.01; for common outcome, lambda=0.05
  Tv1 <- rexp(n, p_t1)
  Z1 <- as.numeric(Tv1 < 10)
  Tv1 <- ifelse(Tv1 > 10, 10, Tv1)
  
  DeathsK.df <- data.frame(ID, A, L, Tv, Z)
  DeathsK.df0 <- data.frame(ID, A=0, L, Tv=Tv0, Z=Z0)
  DeathsK.df1 <- data.frame(ID, A=1, L, Tv=Tv1, Z=Z1)
  head(DeathsK.df)

##################################################################################################
# Point estimates
  
  #Oracle
  oracle <- data.table(rbind(DeathsK.df0,DeathsK.df1))
  oracle.aft <- survreg(Surv(Tv, Z) ~ A, data=oracle, dist="exp")
  sim.res$lhr[1] <- -oracle.aft$coefficients[2]
 
  #IPW with robust SE
  denominator <- rep(NA, n)
  ps <- glm(A ~ L, family=binomial(link="logit"), data=DeathsK.df)$fitted.values
  denominator <- A*ps + (1-A)*(1-ps)
  numerator <- rep(NA, n)
  ps <- glm(A ~ 1, family=binomial(link="logit"), data=DeathsK.df)$fitted.values
  numerator <- A*ps + (1-A)*(1-ps)

  wt <- numerator / denominator
 
  iptw.aft <- survreg(Surv(Tv, Z) ~ A, data=DeathsK.df, dist="exp", weights=wt, cluster=ID)
  sim.res$lhr[2] <- -iptw.aft$coefficients[2]
  sim.res$se[2] <- summary(iptw.aft)$table[2, 2] 
  sim.res$coverage[2] <- (sim.res$lhr[2] - 1.96*sim.res$se[2]) <= true &
	 		  true <= (sim.res$lhr[2] + 1.96*sim.res$se[2]) 

  #Bootstrap the IPW and g-computation
  boot.res <- data.frame(
    method = c("IPW Boot", "G-computation"),
    boot_num=rep(NA, 2),
    lhr=rep(NA, 2),
    stringsAsFactors = FALSE
  )
  
  boot_rep <- function(r) {
      set.seed(r+1)
      boot.res$boot_num=rep(r, 2)
      
      #Bootstrap resample
      index <- sample(1:nrow(DeathsK.df), nrow(DeathsK.df), replace=T)
      if (r==0) {
        boot <- DeathsK.df
      } else {
        boot <- DeathsK.df[index, ]
      }
      
      #IPW
      denominator <- rep(NA, n)
      ps <- glm(A ~ L, family=binomial(link="logit"), data=boot)$fitted.values
      denominator <- A*ps + (1-A)*(1-ps)
      numerator <- rep(NA, n)
      ps <- glm(A ~ 1, family=binomial(link="logit"), data=boot)$fitted.values
      numerator <- A*ps + (1-A)*(1-ps)
      
      wt <- numerator / denominator 
      
      iptw.aft <- survreg(Surv(Tv, Z) ~ A, data=boot, dist="exp", weights=wt)
      boot.res$lhr[1] <- -iptw.aft$coefficients[2]
      
      #G-computation
      mod.D <- survreg(Surv(Tv, Z) ~ A + L, data=boot, dist="exp")
      
      if (montecarlo==0) {      
      	MC <- boot[ , (names(boot) %in% c("L", "A", "rep"))]
      	MC$id <- 1:n
      } else {
      	MC0 <- boot[ ,(names(boot) %in% c("L", "A", "rep"))]
      	index <- sample(1:nrow(MC0), montecarlo, replace=T)
      	MC<-MC0[index, ]
      	MC$id<-1:montecarlo	
      }
      
      pgf<-function(mc_data, exposure){
        d <- mc_data
        id <- d$id
        rep <- d$rep
        
        L <- d$L

        if (is.null(exposure)) {
          A <- d$A
        } else {
          A <- exposure
        }
        
        p_t <- exp(-coef(mod.D)[names(coef(mod.D))=="(Intercept)"])*
          exp(-coef(mod.D)[names(coef(mod.D))=="A"]*A 
              + -coef(mod.D)[names(coef(mod.D))=="L"]*L)
        Tv <- rexp(length(d$id), p_t)
        Z <- ifelse(Tv < 10, 1, 0)
        Tv <- ifelse(Tv > 10, 10, Tv)
        
        gdat <- data.table(id,A,L,Tv,Z,rep)
        return(gdat)
      }
      res0 <- pgf(mc_data=MC, exposure=0)
      res1 <- pgf(mc_data=MC, exposure=1)
      gcomp.dat <- data.table(rbind(res1, res0))
      
      gcomp.aft <- survreg(Surv(Tv, Z) ~ A, data=gcomp.dat, dist="exp")
      boot.res$lhr[2] <- -gcomp.aft$coefficients[2]
      
      return(boot.res)
  }
  
  all.boot <- lapply(0:nboot, function(tt) {boot_rep(tt)})
  all.boot <- do.call(rbind, all.boot)
  
  #For point estimates, pull out results where boot=0
  boot0 <- filter(all.boot, boot_num == 0)
  sim.res$lhr[3:4] <- boot0$lhr
  all.boot <- filter(all.boot, boot_num>0)
  
  #Summarize over bootstraps
  boot.summ <- all.boot %>% 
    group_by(method) %>% 
    summarize(b.avg.lhr = mean(lhr), b.sd.lhr = sd(lhr)) %>% 
    arrange(desc(method))

  sim.res$coverage[3:4] <- (sim.res$lhr[3:4] - 1.96 * boot.summ$b.sd.lhr) <= true & 
    true <= (sim.res$lhr[3:4] + 1.96 * boot.summ$b.sd.lhr)
  sim.res$se[3:4] <- boot.summ$b.sd.lhr
  
  return(sim.res)
  
}

cores <- detectCores() - 2
all.res <- mclapply(1:nsim, function(ii) sim_func(ii), mc.cores=cores, mc.set.seed=FALSE)
all.res <- do.call(rbind,all.res)
all.res$truth <- rep(true, dim(all.res)[1])

filename <- paste("./results/timefix_1-conf_n-", n, "_mc-", montecarlo, sep="")
if (lambda==0.05) {
	filename <- paste(filename, "_common.txt", sep="")
} else {
	filename <- paste(filename, ".txt", sep="")
}
write.table(all.res, file=filename, sep="\t")

