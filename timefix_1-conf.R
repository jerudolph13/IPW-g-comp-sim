
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
# Last Update: 29 Mar 2021
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
  
  truth_sim <- function(exposure) {
    Z1 <- rbinom(simN, 1, 0.3) 
    p_t  <- (lambda)*exp(log(2)*exposure + 0.5*Z1) 
    Tv <- rexp(simN, p_t)
    Y <- as.numeric(Tv < 5)
    Tv <- ifelse(Tv > 5, 5, Tv)
    
    dat <- data.frame(ID, X=exposure, Z1, Tv, Y)
    return(dat)
  }
  
  set.seed(123)
  dat1 <- truth_sim(exposure=1)
  
  set.seed(123)
  dat0 <- truth_sim(exposure=0)
  
  truth_dat <- data.table(rbind(dat0, dat1))
  
  r0 <- mean(truth_dat$Y[truth_dat$X==0])
  r1 <- mean(truth_dat$Y[truth_dat$X==1])
  truth <- r1 - r0
  
  return(truth)
}

true <- truth_func(1e6)

##################################################################################################
## Data generation

sim_func <- function(iter){
  sim.res$sim <- rep(iter, 4)
  ID <- c(1:n)
  
  sim <- function(exposure=NULL) {
    Z1 <- rbinom(n, 1, 0.3)
    
    if (is.null(exposure)){
      p_x <- expit(-log(1/0.5 - 1) - 0.5*0.3 + 0.5*Z1 )
      X <- rbinom(n, 1, p_x)
    } else {
      X <- exposure
    }
    
    p_t  <- (lambda)*exp(log(2)*X + 0.5*Z1)
    Tv <- rexp(n, p_t)
    Y <- as.numeric(Tv < 5)
    Tv <- ifelse(Tv > 5, 5, Tv)
    
    dat <- data.frame(ID, X, Z1, Tv, Y)
    return(dat)
  }

  set.seed(iter)
  DeathsK.df <- sim(exposure=NULL)
  
  set.seed(iter)
  DeathsK.df1 <- sim(exposure=1)
  
  set.seed(iter)
  DeathsK.df0 <- sim(exposure=0)

  
##################################################################################################
## Estimates
  
  # Oracle (no longer reported in the paper)
  oracle <- data.table(rbind(DeathsK.df0,DeathsK.df1))

  sim.res$r0[1] <- mean(oracle$Y[oracle$X==0])
  sim.res$r1[1] <- mean(oracle$Y[oracle$X==1])
  sim.res$rd[1] <- sim.res$r1[1] - sim.res$r0[1]
   
  # Bootstrap to get 95% CIs
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
    boot.res$boot_num <- rep(r, 3)
      
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
    ps <- glm(X ~ Z1, family=binomial(link="logit"), data=boot)$fitted.values
    denominator <- boot$X*ps + (1-boot$X)*(1-ps)
    numerator <- rep(NA, n)
    ps <- glm(X ~ 1, family=binomial(link="logit"), data=boot)$fitted.values
    numerator <- boot$X*ps + (1-boot$X)*(1-ps)
        
    wt <- numerator / denominator 
        
    boot.res$r0[1] <- weighted.mean(boot$Y[boot$X==0], w=wt[boot$X==0])
    boot.res$r1[1] <- weighted.mean(boot$Y[boot$X==1], w=wt[boot$X==1])
    boot.res$rd[1] <- boot.res$r1[1] - boot.res$r0[1]
        
    # Classic g-computation
    mod.D <- survreg(Surv(Tv, Y) ~ X + Z1, data=boot, dist="exp")
        
    if (montecarlo==0) {      
    	MC <- boot[ , (names(boot) %in% c("Z1", "X", "rep"))]
    	MC$id <- 1:n
    } else {
    	MC0 <- boot[ ,(names(boot) %in% c("Z1", "X", "rep"))]
    	index <- sample(1:nrow(MC0), montecarlo, replace=T)
    	MC <- MC0[index, ]
    	MC$id <- 1:montecarlo	
    }
        
    pgf <- function(mc_data, exposure){
          d <- mc_data
          id <- d$id
          rep <- d$rep
          
          Z1 <- d$Z1
  
          if (is.null(exposure)) {
            X <- d$X
          } else {
            X <- exposure
          }
          
          p_t <- exp(-coef(mod.D)[names(coef(mod.D))=="(Intercept)"])*
            exp(-coef(mod.D)[names(coef(mod.D))=="X"]*X 
                - coef(mod.D)[names(coef(mod.D))=="Z1"]*Z1)
          Tv <- rexp(length(d$id), p_t)
          Y <- ifelse(Tv < 5, 1, 0)
          Tv <- ifelse(Tv > 5, 5, Tv)
          
          gdat <- data.table(id,X,Z1,Tv,Y,rep)
          return(gdat)
    }
    
    set.seed(r+2)
    res0 <- pgf(mc_data=MC, exposure=0)
    set.seed(r+2)
    res1 <- pgf(mc_data=MC, exposure=1)
    gcomp.dat <- data.table(rbind(res1, res0))
        
    boot.res$r0[2] <- mean(gcomp.dat$Y[gcomp.dat$X==0])
    boot.res$r1[2] <- mean(gcomp.dat$Y[gcomp.dat$X==1])
    boot.res$rd[2] <- boot.res$r1[2] - boot.res$r0[2]
      
    # G-computation via iterated conditional expectations (ICE)
    outcome <- boot %>% 
      mutate(time = ceiling(Tv)) %>% 
      select(bid, time, Y) %>% 
      pivot_wider(names_from=time, names_prefix="Y", values_from=Y, values_fill=0, names_sort=T)
    
      # Deal with the possibility that no one might have an event at a particular t
      out <- data.frame(Y1=rep(0, nrow(outcome)),
                        Y2=rep(0, nrow(outcome)),
                        Y3=rep(0, nrow(outcome)),
                        Y4=rep(0, nrow(outcome)),
                        Y5=rep(0, nrow(outcome)))
      miss <- names(out)[!(names(out) %in% names(outcome))]
      out2 <- data.frame(out[ , miss])
      names(out2) <- names(out)[!(names(out) %in% names(outcome))]
      outcome <- bind_cols(outcome, out2) %>% 
        select(str_sort(c("bid", names(out)), numeric=TRUE))
      
      # Once an event has occurred, all subsequent nodes must be 1
      for (i in 2:5) {
        outcome[ , i+1] <- outcome[ , i+1] + outcome[ , i]
      }
      
      # Merge back to baseline data
      base <- boot[ , c("bid", "Z1", "X")]
      dat <- left_join(base, outcome, by="bid")
      
      # Set up call to ltmle package
      Anodes <- "X"
      Ynodes <- vars_select(names(dat), starts_with("Y"))
      
      # Use ltmle to implement ICE g-comp
      res <- ltmle(data=select(dat, -bid), 
                   Anodes=Anodes, Ynodes=Ynodes,
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

# Output results
filename <- paste("../results/timefix_1-conf_n-", n, "_mc-", montecarlo, sep="")
if (lambda==0.05) {
	filename <- paste(filename, "_common.txt", sep="")
} else {
	filename <- paste(filename, ".txt", sep="")
}
write.table(all.res, file=filename, sep="\t")

