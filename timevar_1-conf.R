
###################################################################################################
#
# Purpose: Compare using simulation the performance of g-comp and IPTW to estimate same 
#          parameter in presence of:
#           -time-varying exposure
#           -1 time-varying confounder
#           -binary outcome
#
# Authors: Jacqueline Rudolph, Ashley Naimi (Credit to Young and Moodie for DGM)
#
# Last Update: 29 Mar 2021
#
##################################################################################################

# Read in packages
lib <- "~/R/x86_64-pc-linux-gnu-library/4.0"
packages <- c("tidyverse", "data.table", "survival", "parallel", "flexsurv", "tidyselect", "ltmle")
for (package in packages) {
  library(package, character.only=T, lib.loc=lib)
}

# Pull in command line arguments
args <- commandArgs(trailingOnly=TRUE)

# Define parameters and functions
nsim <- 1000                            #Number of simulations
nboot <- 200                            #Number of bootstrap resamples
n <- as.numeric(args[1])                #Sample size
N <- 5                                  #Number of time points
k <- 1                                  #Number of outcomes
lambda <- as.numeric(args[2])           #Baseline rate
montecarlo <- as.numeric(args[3])       #Monte Carlo resample size (0 implies no MC)

expit <- function(x) {1/(1+exp(-x))}

# Read in the truth
truth <- read.table(file="../results/timevar_truth.txt", sep="\t") %>% 
  rename(lam=lambda) %>% 
  filter(n.conf==1 & lam==lambda)

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
## Data generation

## This code generates data from a structural nested model 
## compatible with a marginal structural model.
## It is based on Jessica Young's algorithm, published here:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3635680/
#
## This code, which extend's Young's algorithm to multiple outcomes
## was written by Erica Moodie, published here:
## https://www.ncbi.nlm.nih.gov/pubmed/24272681

simloop <- function(s, nboot, montecarlo){
  sim.res$sim <- s
  set.seed(s)
  
  # This is the matrix of parameters of interest, possibly different at each interval
  psi.mat <- matrix(0, nrow=k, ncol=N+1)
  
  # Here are the effect sizes for the K=1 cause
  psi.mat[1, ] <- log(2)
  
  # Here the (untreated) all-cause rate is set to lambda
  gamma.vec <- rep(log(lambda/k))
  muK <- sum(exp(gamma.vec))
  X<-Z1<-ID<-Y<-K<-Tv<-Int<-XLast<-Z1Last<-Z1First <- numeric()
  T0.vec<-T.vec<-Y.vec<-K.vec <- rep(0, n)
  
  # Here are the coefficients determining the mediation and treatment assignment mechanisms
  bevec <- c(log(3/7), 2, log(0.5), log(1.5)) # Used to generate time-varying confounder Z1
  alvec <- c(log(2/7), 0.5, 0.5, log(4)) # Used to generate exposure (Intercept, Z1, Z1Last, XLast)
  
  # cval is used to introduce the confounding
  cval <- 30
  
  # Begin the data-generation loop
  simulation <- function (exposure) {
      
    for (i in 1:n) {
      # Generate the counterfactual (untreated) survival time
      T0 <- rexp(1, lambda) # Generate T0 from an exponential dist with constant rate=lamba
      Ival <- as.numeric(T0 < cval)
      # Begin the interval-by-interval simulation
      m <- 0
      mu.tot <- 0
      X.vec<-Z1.vec<-XLast.vec<-Z1Last.vec<-Z1First.vec <- rep(0, N+1)
      # Implement Young's algorithm with multiple causes
      # Generate the survival time, then the cause
      while (muK*T0 > mu.tot & m <= N) {
        if (m == 0) {
          # First interval
          eta <- bevec[1] + bevec[2]*Ival + bevec[3]*0 + bevec[4]*0
          pval <- 1 / (1 + exp(-eta))
          Z1.vec[m+1] <- rbinom(1, 1, pval)
          
          eta <- alvec[1] + alvec[2]*Z1.vec[m+1] + alvec[3]*0 + alvec[4]*0
          pval <- 1 / (1 + exp(-eta))
          if (is.null(exposure)) {X.vec[m + 1] <- rbinom(1, 1, pval)}
          else {X.vec[m+1] <- exposure}
          
          XLast.vec[m+1] <- 0; Z1Last.vec[m+1] <- 0
          Z1First.vec <- rep(Z1.vec[m+1], N + 1)
        } else {
          # Subsequent intervals
          eta <- bevec[1] + bevec[2]*Ival + bevec[3]*X.vec[m] + bevec[4]*Z1.vec[m]
          pval <- 1 / (1 + exp(-eta))
          Z1.vec[m+1] <- rbinom(1, 1, pval) 
 
          eta <- alvec[1] + alvec[2]*Z1.vec[m+1] + alvec[3]*Z1.vec[m] + alvec[4]*X.vec[m] 
          pval <- 1 / (1 + exp(-eta))
          if (is.null(exposure)) {X.vec[m+1] <- rbinom(1, 1, pval)}
          else {X.vec[m+1] <- exposure}
          XLast.vec[m+1] <- X.vec[m]; Z1Last.vec[m+1] <- Z1.vec[m]
        }
        
        muval <- sum(exp(gamma.vec + X.vec[m+1]*psi.mat[ , m+1]))
        
        # Tval is computed for each interval, but is overwritten until the final interval
        Tval <- m + (muK * T0 - mu.tot) / muval
        mu.tot <- mu.tot + muval
        m <- m + 1
      }
      
      # After exiting the loop, the survival time has been generated as Tval
      # Now need to generate the failure type.
      if (m > N) {
        # In the case of censoring at tenth interval, no failure.
        Tval <- m - 1
        Y.vec[i] <- 0
      } else {
        # In the case of failure, use the ratio hazards to define the
        # relevant multinomial distribution on the k causes.
        Y.vec[i] <- sample(c(1:k), 1, prob = exp(gamma.vec + X.vec[m]*psi.mat[ ,m]))
      }
      
      # Store the outcomes
      T0.vec[i] <- T0
      T.vec[i] <- Tval
      K.vec[i] <- m - 1
      ID <- c(ID, rep(i,m)) # Individual
      Int <- c(Int, c(1:m)) # Time point
      X <- c(X, X.vec[1:m]) # Time-updated treatment
      Z1 <- c(Z1, Z1.vec[1:m]) # Time-updated covariate Z1
      XLast <- c(XLast, XLast.vec[1:m]) # Treatment at last t point
      Z1Last <- c(Z1Last, Z1Last.vec[1:m]) # Covariate Z1 at last t point
      Z1First <- c(Z1First, Z1First.vec[1:m]) # Baseline covariate Z1 value
      Y <- c(Y, rep(0,m - 1), Y.vec[i]) # Outcome: Y>0 indicates outcome of some type, determined by value)
      tv <- c(1:m); tv[m] <- Tval
      Tv <- c(Tv, tv) # If event occurs, exact time at which it occurred; o.w. equal to Int)
    }
    
    DeathsK.df <- data.frame(ID, Int, Tv, X, XLast, Z1, Z1Last, Z1First, Y)
    
    # Trim off the intervals beyond the Nth (loop goes one too far)
    DeathsK.df <- DeathsK.df[DeathsK.df$Int <= N, ]
    DeathsK.df$Int0 <- DeathsK.df$Int - 1
    
    return(DeathsK.df)
  }


  # Oracle: get counterfactuals for every simulated individual by setting exposure inside simulation
    # This gives us a measure of "optimal" performance given sample size and model
  set.seed(123+s)
  oracle1 <- simulation(exposure=1)
    oracle1$last <- as.numeric(!duplicated(oracle1$ID, fromLast=T))
  set.seed(123+s) 
  oracle0 <- simulation(exposure=0)
    oracle0$last <- as.numeric(!duplicated(oracle0$ID, fromLast=T))
  
  sim.res$r0[1] <- mean(oracle0$Y[oracle0$last==1])
  sim.res$r1[1] <- mean(oracle1$Y[oracle1$last==1])
  sim.res$rd[1] <- sim.res$r1[1] - sim.res$r0[1]
  
  # Now create simulation to be used in analysis steps
  set.seed(123+s)
  DeathsK.df <- simulation(exposure=NULL)


##################################################################################################
## Bootstrap resample to get CIs

  # Set up data set to hold bootstrap results
  boot.res <- data.frame(
    method=c("IPW", "G-computation", "ICE"),
    boot_num=rep(NA, 3),
    r1=rep(NA, 3),
    r0=rep(NA, 3),
    rd=rep(NA, 3),
    stringsAsFactors=FALSE
  )
    
  bootrep <- function(r) {
    boot.res$boot_num <- r
    set.seed(r+1)
    firstobs <- DeathsK.df[DeathsK.df$Int == 1, ]
    samp <- table(firstobs[sample(1:nrow(firstobs),nrow(firstobs),replace=T), (names(DeathsK.df) == "ID")])
  
  # The step below pulls in the simulated data for boot=0; otherwise grabs all records for the resampled observations
    boot <- NULL
    if(r==0){
      boot <- DeathsK.df %>% 
        rename(bid = ID)
    } else{
      for(zzz in 1:max(samp)){ 
        cc <- DeathsK.df[DeathsK.df$ID %in% names(samp[samp %in% c(zzz:max(samp))]),]
        cc$bid <- paste0(cc$ID, zzz)
        boot <- rbind(boot, cc)
      }
      boot <- select(boot, -ID)
    }
  
  
##################################################################################################
## IPW 
  
    # Denominator of weights
    ps <- glm(X ~ XLast + Z1 + Z1Last + as.factor(Int), family=binomial(link="logit"), data=boot)$fitted.values
    boot$denominator <- boot$X*ps + (1-boot$X)*(1-ps)
    
    # Numerator of weights
    ps <- glm(X ~ as.factor(Int), family=binomial(link="logit"), data=boot)$fitted.values
    boot$numerator <- boot$X*ps + (1-boot$X)*(1-ps)
    
    boot <- boot %>% 
      group_by(bid) %>%  
      mutate(wt=cumprod(numerator/denominator)) %>% 
      ungroup(bid) 

    # IP-weighted survival
    fit <- summary(survfit(Surv(Int0, Tv, Y)  ~ X, data=boot, weights=wt))
    surv <- data.frame(time = fit$time, 
                       surv = fit$surv,
                       exposure = fit$strata)
    boot.res$r0[1] <- 1 - min(surv$surv[surv$exposure=="X=0"])
    boot.res$r1[1] <- 1 - min(surv$surv[surv$exposure=="X=1"])
    boot.res$rd[1] <- boot.res$r1[1] - boot.res$r0[1]
  

##################################################################################################
## Classic g-computation
    # Time-ordering: Z1Last, XLast, Z1, X, Y

    # Model confounder
    mod.Z1 <- glm(Z1 ~ Z1Last + XLast + as.factor(Int), family=binomial(link="logit"), data=boot)

    # Model outcome (flexsurv used bc survreg doesn't support start/stop coding)
    mod.D <- flexsurvreg(Surv(Int0,Tv,Y) ~ X + XLast + Z1 + Z1Last, data=boot, dist="exp")
    
    # Take a Monte Carlo (MC) sample
    # Select first obs for each person to obtain joint empirical distribution of baseline covariates
    MC0 <- boot[boot$Int==1, (names(boot) %in% c("Z1", "X"))]
    index <- sample(1:nrow(MC0), montecarlo, replace=T)
    MC <- MC0[index, ]
    MC$id <- 1:montecarlo
    
    # Predict follow-up based on g-formula using PGF function
    pgf <- function (ii, mc_data, length, exposure=NULL) {
      
      pFunc <- function (mod,ndat) {
        as.numeric(predict(mod, newdata=ndat, type="response") > runif(1))
      }
      
      d <- mc_data
      d <- d[d$id==ii,]
      lngth <- length
      Z1p<-Xp<-Yp<-time <- numeric()
      time[1] <- j <- 1
      id <- d$id
  
      Z1p[1] <- d$Z1
      Xp[1] <- exposure

      # event status at first time point
      dYp <- data.table(X=Xp[1], XLast=0, Z1=Z1p[1], Z1Last=0)
      expSim <- function (dat) {
        newD <- dat
        desX <- newD[,c("X", "XLast", "Z1","Z1Last")]
        y0 <- rexp(1, exp(coef(mod.D)[names(coef(mod.D))=="rate"])*
                   exp(coef(mod.D)[!names(coef(mod.D))=="rate"]%*%t(desX)))
        return(y0)
      }
      t0 <- expSim(dYp)
      if (t0<=1) {
        Yp[1] <- 1
        jj <- t0
      } else {
        Yp[1] <- 0
        jj <- 1
      }
      
      # subsequent time points
      for (j in 2:lngth) {
        if (Yp[j-1]==0) {
          XLast <- Xp[j-1]; Z1Last <- Z1p[j-1]
          
          dZ1p <- data.table(Z1Last, XLast, Int=factor(j))
          Z1p[j] <- pFunc(mod.Z1, dZ1p)
  
          Xp[j] <- exposure
          
          dYp <- data.table(X=Xp[j], XLast, Z1=Z1p[j], Z1Last)
          t0 <- expSim(dYp)
          if (t0<=0.001) {
            Yp[j-1] <- 1 #If the time interval is short enough to cause an error, place event at end of previous time point
            Z1p <- Z1p[1:(j-1)]
            Xp <- Xp[1:(j-1)]
            id <- id[1:(j-1)]
            time <- time[1:(j-1)]
            jj <- jj[1:(j-1)]
            break
          } else {
            if (t0>0.001 & t0<=1) {
              Yp[j] <- 1
              jj <- (j-1) + t0
            } else { 
              Yp[j] <- 0
              jj <- j
            }}
          
        } else {
          break
        }
        time[j] <- j
      }
      gdat <- data.table(id, time, jj, Xp, Z1p, Yp)
      gdat$last <- as.numeric(gdat$Yp!=0 | gdat$time==lngth)
      return(gdat)
    }
  
    set.seed(r+2)
    res0 <- lapply(1:montecarlo,function(x) {pgf(x ,mc_data=MC, length=N, exposure=0)})
      res0 <- do.call(rbind, res0)
    set.seed(r+2)
    res1 <- lapply(1:montecarlo,function(x) {pgf(x, mc_data=MC, length=N, exposure=1)})
      res1 <- do.call(rbind, res1)
    
    # Estimate risk
    boot.res$r0[2] <- mean(res0$Yp[res0$last==1])
    boot.res$r1[2] <- mean(res1$Yp[res1$last==1])
    boot.res$rd[2] <- boot.res$r1[2] - boot.res$r0[2]

    
##################################################################################################
## G-computation via iterated conditional expectations (ICE)
    
    # Make confounder wide
    confounder <- boot %>% 
      select(bid, Int, Z1) %>% 
      pivot_wider(names_from=Int, names_prefix="Z1", values_from=Z1, values_fill=NA, names_sort=T) 
    
    # Make exposure wide
    exposure <- boot %>% 
      select(bid, Int, X) %>% 
      pivot_wider(names_from=Int, names_prefix="X", values_from=X, values_fill=NA, names_sort=T)  
    
    # Make outcome wide
    outcome <- boot %>% 
      select(bid, Int, Y) %>% 
      pivot_wider(names_from=Int, names_prefix="Y", values_from=Y, values_fill=0, names_sort=T)  
      # Once an event has occurred, all subsequent nodes must be 1
      for (i in 2:N) {
        outcome[ , i+1] <- outcome[ , i+1] + outcome[ , i]
      }
    
    # Interleave together in correct order (Z1, X, Y)
    wide <- select(boot, bid)
    for (i in 1:N){
      wide <- merge(wide, confounder[ , c(1, i+1)], by="bid")
      wide <- merge(wide, exposure[ , c(1, i+1)], by="bid")
      wide <- merge(wide, outcome[ , c(1, i+1)], by="bid")
    }
    
    # Set up call to ltmle package
    Anodes <- vars_select(names(wide), starts_with("X"))
    Lnodes <- vars_select(names(wide), starts_with("Z1"))
    Ynodes <- vars_select(names(wide), starts_with("Y"))
    
    Xformula <- rep(NA, N)
    for (i in 1:N){
      if (i==1) {
        #Formula for X(1)
        Xformula[i] <- paste(Anodes[i], "~", Lnodes[i], sep=" ")
      } else {
        #Formula for X(t), t>1
        Xformula[i] <- paste(Anodes[i], "~", Anodes[i-1], "+", Lnodes[i], "+", Lnodes[i-1], sep=" ")
      }
    }
    
    Yformula <- rep(NA, N)
    for (i in 1:N){
      if (i==1) {
        #Formula for Y(1)
        names(Yformula)[i] <- Ynodes[i]
        Yformula[i] <- paste("Q.kplus1 ~", Anodes[i], "+", Lnodes[i], sep=" ")
      } else {
        #Formula for Y(t), t>1
        names(Yformula)[i] <- Ynodes[i]
        Yformula[i] <- paste("Q.kplus1 ~", Anodes[i], "+", Anodes[i-1], "+", Lnodes[i], "+", Lnodes[i-1], sep=" ")
      }
    }

    # Use ltmle to implement ICE g-comp
    res <- ltmle(data=select(wide, -bid), 
                 Anodes=Anodes, Lnodes=Lnodes, Ynodes=Ynodes,
                 Qform=Yformula,
                 gform=Xformula,
                 abar=list(treatment=rep(1, N), control=rep(0, N)),
                 survivalOutcome=T,
                 SL.library=NULL,
                 gcomp=T)
    
    summ <- summary(res)  
    boot.res$r1[3] <- summ$effect.measures$treatment$estimate
    boot.res$r0[3] <- summ$effect.measures$control$estimate
    boot.res$rd[3] <- summ$effect.measures$ATE$estimate
    
    return(boot.res)
  }

  all.boot <- lapply(0:nboot, function(tt) {bootrep(tt)})
  all.boot <- do.call(rbind, all.boot)


##################################################################################################
## Aggregate results
  
  # For point estimate, pull out results where boot=0
  boot0 <- filter(all.boot, boot_num==0)
  sim.res$r0[2:4] <- boot0$r0
  sim.res$r1[2:4] <- boot0$r1
  sim.res$rd[2:4] <- boot0$rd
  all.boot <- filter(all.boot, boot_num>0)

  # Summarize over bootstraps
  boot.summ <- all.boot %>% 
    group_by(method) %>% 
    summarize(se = sd(rd))

  sim.res <- left_join(sim.res, boot.summ, by="method") %>% 
    mutate(coverage = (rd - 1.96*se) <= truth$rd & truth$rd <= (rd + 1.96*se))

  return(sim.res)
}

cores <- detectCores() - 2
all.res <- mclapply(1:nsim, function(x) {simloop(x, nboot, montecarlo)}, mc.cores=cores, mc.set.seed=FALSE)
  # Use lapply to test the code if there are errors
  #all.res <- lapply(1:nsim, function(x) {simloop(x, nboot, montecarlo)})
all.res <- do.call(rbind, all.res)
all.res$truth <- truth$rd

filename <- paste("../results/timevar_1-conf_n-", n, "_mc-", montecarlo, sep="")
if (lambda==0.05) {
  filename <- paste(filename, "_common.txt", sep="")
} else {
  filename <- paste(filename, ".txt", sep="")
}
write.table(all.res, file=filename, sep="\t")

