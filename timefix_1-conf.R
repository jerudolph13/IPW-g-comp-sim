
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
# Last Update: 20 Jan 2021
#
##################################################################################################

# Read in packages
lib <- "~/R/x86_64-pc-linux-gnu-library/4.0"
packages <- c("tidyverse", "data.table", "survival", "parallel")
for (package in packages) {
  library(package, character.only=T, lib.loc=lib)
}

# Pull in command line arguments
args <- commandArgs(trailingOnly=TRUE)

# Define parameters and functions
nsim <- 1000			#Number of simulations
nboot <- 200			#Number of bootstrap resamples
n <- as.numeric(args[1])		#Sample size
lambda <- as.numeric(args[2])		#Baseline rate
montecarlo <- as.numeric(args[3])	#Monte Carlo resample size (0 implies no MC)

expit <- function(x) {1/(1+exp(-x))}

# Prepare data set to hold simulation results
sim.res <- data.frame(
  method=c("Oracle", "IPW", "G-computation"),
  r1=rep(NA, 3),
  r0=rep(NA, 3),
  rd=rep(NA, 3),
  se=rep(NA, 3),
  coverage=rep(NA, 3),
  sim=rep(NA, 3),
  stringsAsFactors = FALSE
)


##################################################################################################
#Get the true, marginal RD

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
  
  fit <- summary(survfit(Surv(Tv, Z) ~ A, data=truth_dat))
  truth.surv <- data.frame(time = fit$time, 
                           surv = fit$surv,
                           exposure = fit$strata)
  r0 <- 1- min(truth.surv$surv[truth.surv$exposure=="A=0"])
  r1 <- 1- min(truth.surv$surv[truth.surv$exposure=="A=1"])
  truth <- r1 - r0
  
  return(truth)
}

set.seed(123)
true <- truth_func(1e6)


##################################################################################################
# Data generation

sim_func <- function(iter){
  set.seed(iter)
  sim.res$sim <- rep(iter, 3)

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

  
##################################################################################################
# Point estimates
  
  #Oracle
  oracle <- data.table(rbind(DeathsK.df0,DeathsK.df1))
  
  fit <- summary(survfit(Surv(Tv, Z) ~ A, data=oracle))
  surv <- data.frame(time = fit$time, 
                     surv = fit$surv,
                     exposure = fit$strata)
  sim.res$r0[1] <- 1 - min(surv$surv[surv$exposure=="A=0"])
  sim.res$r1[1] <- 1 - min(surv$surv[surv$exposure=="A=1"])
  sim.res$rd[1] <- sim.res$r1[1] - sim.res$r0[1]
 
  #Bootstrap the IPW and g-computation
  boot.res <- data.frame(
    method = c("IPW", "G-computation"),
    boot_num=rep(NA, 2),
    r0=rep(NA, 2),
    r1=rep(NA, 2),
    rd=rep(NA, 2),
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
      
      fit <- summary(survfit(Surv(Tv, Z) ~ A, data=boot, weights=wt))
      surv <- data.frame(time = fit$time, 
                         surv = fit$surv,
                         exposure = fit$strata)
      boot.res$r0[1] <- 1 - min(surv$surv[surv$exposure=="A=0"])
      boot.res$r1[1] <- 1 - min(surv$surv[surv$exposure=="A=1"])
      boot.res$rd[1] <- boot.res$r1[1] - boot.res$r0[1]
      
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
              - coef(mod.D)[names(coef(mod.D))=="L"]*L)
        Tv <- rexp(length(d$id), p_t)
        Z <- ifelse(Tv < 10, 1, 0)
        Tv <- ifelse(Tv > 10, 10, Tv)
        
        gdat <- data.table(id,A,L,Tv,Z,rep)
        return(gdat)
      }
      res0 <- pgf(mc_data=MC, exposure=0)
      res1 <- pgf(mc_data=MC, exposure=1)
      gcomp.dat <- data.table(rbind(res1, res0))
      
      fit <- summary(survfit(Surv(Tv, Z) ~ A, data=gcomp.dat))
      surv <- data.frame(time = fit$time, 
                         surv = fit$surv,
                         exposure = fit$strata)
      boot.res$r0[2] <- 1 - min(surv$surv[surv$exposure=="A=0"])
      boot.res$r1[2] <- 1 - min(surv$surv[surv$exposure=="A=1"])
      boot.res$rd[2] <- boot.res$r1[2] - boot.res$r0[2]

      return(boot.res)
  }
  
  all.boot <- lapply(0:nboot, function(tt) {boot_rep(tt)})
  all.boot <- do.call(rbind, all.boot)
  
  #For point estimates, pull out results where boot=0
  boot0 <- filter(all.boot, boot_num == 0)
  sim.res$r0[2:3] <- boot0$r0
  sim.res$r1[2:3] <- boot0$r1
  sim.res$rd[2:3] <- boot0$rd
  all.boot <- filter(all.boot, boot_num>0)
  
  #Summarize over bootstraps
  boot.summ <- all.boot %>% 
    group_by(method) %>% 
    summarize(b.avg.rd = mean(rd), b.sd.rd = sd(rd)) %>% 
    arrange(desc(method))

  sim.res$coverage[2:3] <- (sim.res$rd[2:3] - 1.96 * boot.summ$b.sd.rd) <= true & 
    true <= (sim.res$rd[2:3] + 1.96 * boot.summ$b.sd.rd)
  sim.res$se[2:3] <- boot.summ$b.sd.rd
  
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

