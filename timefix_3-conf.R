
###################################################################################################
#
# Purpose: Compare using simulation the performance of g-comp and IPTW to estimate same 
#          parameter in presence of:
#           -time-fixed exposure
#           -3 time-fixed confounders
#           -binary outcome
#
# Author: Jacqueline Rudolph and Ashley Naimi
#
# Last Update: 25 Jan 2021
#
##################################################################################################

# Read in packages
lib <- "~/R/x86_64-pc-linux-gnu-library/4.0"
packages <- c("tidyverse", "data.table", "survival", "parallel", "ltmle", "tidyselect")
for (package in packages) {
  library(package, character.only=T, lib.loc=lib)
}

# Pull in command line arguments
args <- commandArgs(trailingOnly=TRUE)

# Define parameters and functions
nsim <- 1000			                #Number of simulations
nboot <- 200			                #Number of bootstrap resamples
n <- as.numeric(args[1])		      #Sample size
lambda <- as.numeric(args[2])		  #Baseline rate
montecarlo <- as.numeric(args[3])	#Monte Carlo resample size (0 implies no MC)

expit <- function(x) {1/(1+exp(-x))}

# Prepare data set to hold simulation results
sim.res <- data.frame(
  method=c("Oracle", "IPW", "G-computation", "ICE"),
  r1=rep(NA, 4),
  r0=rep(NA, 4),
  rd=rep(NA, 4),
  sim=rep(NA, 4),
  stringsAsFactors = FALSE
)


##################################################################################################
## Get the true, marginal RD

truth_func <- function(simN){
  ID <- c(1:simN)
  L <- rbinom(simN, 1, 0.3) 
  J <- rbinom(simN, 1, 0.4)
  M <- rbinom(simN, 1, 0.5) 
  
  p_t0  <- (lambda)*exp(log(2)*0 + 0.5*L + 0.6*J + 0.8*M) #In main analysis, lambda=0.01; for common outcome, lambda=0.05
  Tv0 <- rexp(simN, p_t0)
  Z0 <- as.numeric(Tv0 < 10)
  Tv0 <- ifelse(Tv0 > 10, 10, Tv0)
  
  p_t1  <- (lambda)*exp(log(2)*1 + 0.5*L + 0.6*J + 0.8*M) #In main analysis, lambda=0.01; for common outcome, lambda=0.05
  Tv1 <- rexp(simN, p_t1)
  Z1 <- as.numeric(Tv1 < 10)
  Tv1 <- ifelse(Tv1 > 10, 10, Tv1)
  
  DeathsK.df0 <- data.frame(ID, A=0, L, J, M, Tv=Tv0, Z=Z0)
  DeathsK.df1 <- data.frame(ID, A=1, L, J, M, Tv=Tv1, Z=Z1)
  
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
## Data generation

sim_func <- function(iter){
  set.seed(iter)
  sim.res$sim <- rep(iter, 4)

  ID <- c(1:n)
  L <- rbinom(n, 1, 0.3) 
  J <- rbinom(n, 1, 0.4)
  M <- rbinom(n, 1, 0.5) 
  
  p_a <- expit(-log(1/0.5 - 1) + 0.5*L - 0.5*0.3 + 0.6*J - 0.6*0.4 + 0.8*M - 0.8*0.5)
  A <- rbinom(n, 1, p_a)
  
  p_t  <- (lambda)*exp(log(2)*A + 0.5*L + 0.6*J + 0.8*M) #In main analysis, lambda=0.01; for common outcome, lambda=0.05
  Tv <- rexp(n, p_t)
  Z <- as.numeric(Tv < 10)
  Tv <- ifelse(Tv > 10, 10, Tv)
  
  p_t0  <- (lambda)*exp(log(2)*0 + 0.5*L + 0.6*J + 0.8*M) #In main analysis, lambda=0.01; for common outcome, lambda=0.05
  Tv0 <- rexp(n, p_t0)
  Z0 <- as.numeric(Tv0 < 10)
  Tv0 <- ifelse(Tv0 > 10, 10, Tv0)
  
  p_t1  <- (lambda)*exp(log(2)*1 + 0.5*L + 0.6*J + 0.8*M) #In main analysis, lambda=0.01; for common outcome, lambda=0.05
  Tv1 <- rexp(n, p_t1)
  Z1 <- as.numeric(Tv1 < 10)
  Tv1 <- ifelse(Tv1 > 10, 10, Tv1)
  
  DeathsK.df <- data.frame(ID, L, J, M, A, Tv, Z)
  DeathsK.df0 <- data.frame(ID, L, J, M, A=0, Tv=Tv0, Z=Z0)
  DeathsK.df1 <- data.frame(ID, L, J, M, A=1, Tv=Tv1, Z=Z1)
  
##################################################################################################
## Estimates
  
  # Oracle
  oracle <- data.table(rbind(DeathsK.df0,DeathsK.df1))
  
  fit <- summary(survfit(Surv(Tv, Z) ~ A, data=oracle))
  surv <- data.frame(time = fit$time, 
                     surv = fit$surv,
                     exposure = fit$strata)
  sim.res$r0[1] <- 1 - min(surv$surv[surv$exposure=="A=0"])
  sim.res$r1[1] <- 1 - min(surv$surv[surv$exposure=="A=1"])
  sim.res$rd[1] <- sim.res$r1[1] - sim.res$r0[1]
  
  # Bootstrap the IPW and g-computation
  boot.res <- data.frame(
    method = c("IPW", "G-computation", "ICE"),
    boot_num=rep(NA, 3),
    r0=rep(NA, 3),
    r1=rep(NA, 3),
    rd=rep(NA, 3),
    stringsAsFactors = FALSE
  )
  
  boot_rep <- function(r) {
    set.seed(r+1)
    boot.res$boot_num=rep(r, 3)
      
    # Bootstrap resample
    index <- sample(1:nrow(DeathsK.df), nrow(DeathsK.df), replace=T)
    if (r==0) {
      boot <- DeathsK.df %>% 
        rename(bid = ID)
    } else {
      boot <- DeathsK.df[index, ] %>% 
        mutate(bid = c(1:n)) %>% 
        select(-ID)
    }
      
    # IPW
    denominator <- rep(NA, n)
    ps <- glm(A ~ L + J + M, family=binomial(link="logit"), data=boot)$fitted.values
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
      
    # Classic g-computation
    mod.D <- survreg(Surv(Tv, Z) ~ A + L + J + M, data=boot, dist="exp")
      
    if (montecarlo==0) {      
      MC <- boot[ , (names(boot) %in% c("L", "J", "M", "A", "rep"))]
      MC$id <- 1:n
    } else {
      MC0 <- boot[ ,(names(boot) %in% c("L", "J", "M", "A", "rep"))]
      index <- sample(1:nrow(MC0), montecarlo, replace=T)
      MC<-MC0[index, ]
      MC$id<-1:montecarlo	
    }
      
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
        
         p_t <- exp(-coef(mod.D)[names(coef(mod.D))=="(Intercept)"])*
          exp(-coef(mod.D)[names(coef(mod.D))=="A"]*A 
              - coef(mod.D)[names(coef(mod.D))=="L"]*L
              - coef(mod.D)[names(coef(mod.D))=="J"]*J
              - coef(mod.D)[names(coef(mod.D))=="M"]*M)
        Tv <- rexp(length(d$id), p_t)
        Z <- ifelse(Tv < 10, 1, 0)
        Tv <- ifelse(Tv > 10, 10, Tv)
        
        gdat <- data.table(id,A,L,J,M,Tv,Z,rep)
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
    
    # G-computation via iterated conditional expectations (ICE)
    outcome <- boot %>% 
      mutate(time = ceiling(Tv)) %>% 
      select(bid, time, Z) %>% 
      pivot_wider(names_from=time, names_prefix="Z", values_from=Z, values_fill=0, names_sort=T)  
    
      # Once an event has occurred, all subsequent nodes must be 1
      for (i in 2:10) {
        outcome[ , i+1] <- outcome[ , i+1] + outcome[ , i]
      }
    
      # Merge back to baseline data   
      dat <- left_join(select(boot, -c(Tv, Z)), outcome, by="bid")
    
      # Set up call to ltmle package
      Anodes <- "A"
      Lnodes <- c("L", "J", "M")
      Ynodes <- vars_select(names(dat), starts_with("Z"))
    
      # Use ltmle to implement ICE g-comp
      res <- ltmle(data=select(dat, -bid), 
                   Anodes=Anodes, Lnodes=Lnodes, Ynodes=Ynodes,
                   abar=list(treatment=1, control=0),
                   survivalOutcome=T,
                   SL.library=NULL,
                   gcomp=T)
    
      summ <- summary(res)  
      boot.res$r1[3] <- summ$effect.measures$treatment$estimate
      boot.res$r0[3] <- summ$effect.measures$control$estimate
      boot.res$rd[3] <- summ$effect.measures$ATE$estimate
      
    return(boot.res)
  }
  
  all.boot <- lapply(0:nboot, function(tt) {boot_rep(tt)})
  all.boot <- do.call(rbind, all.boot)
  
  
##################################################################################################
## Aggregate results
  
  # For point estimates, pull out results where boot=0
  boot0 <- filter(all.boot, boot_num == 0)
  sim.res$r0[2:4] <- boot0$r0
  sim.res$r1[2:4] <- boot0$r1
  sim.res$rd[2:4] <- boot0$rd
  all.boot <- filter(all.boot, boot_num>0)
  
  # Summarize over bootstraps
  boot.summ <- all.boot %>% 
    group_by(method) %>% 
    summarize(se = sd(rd))
  
  sim.res <- left_join(sim.res, boot.summ, by="method") %>% 
    mutate(coverage = (rd - 1.96*se) <= true & true <= (rd + 1.96*se))
  
  return(sim.res)
}

cores <- detectCores() - 2
all.res <- mclapply(1:nsim, function(ii) sim_func(ii), mc.cores=cores, mc.set.seed=FALSE)
  #Use lapply when testing code for errors
  #all.res <- lapply(1:nsim, function(ii) sim_func(ii))
all.res <- do.call(rbind,all.res)
all.res$truth <- rep(true, dim(all.res)[1])

filename <- paste("./results/timefix_3-conf_n-", n, "_mc-", montecarlo, sep="")
if (lambda==0.05) {
  filename <- paste(filename, "_common.txt", sep="")
} else {
  filename <- paste(filename, ".txt", sep="")
}
write.table(all.res, file=filename, sep="\t")

