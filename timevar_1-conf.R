
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
# Last Update: 25 Jan 2021
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
N <- 10                                 #Number of time points
K <- 1                                  #Number of outcomes
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
  psi.mat <- matrix(0, nrow=K, ncol=N+1)
  
  # Here are the effect sizes for the K=1 cause
  psi.mat[1, ] <- log(2)
  
  # Here the (untreated) all-cause rate is set to lambda
  gamma.vec <- rep(log(lambda/K))
  muK <- sum(exp(gamma.vec))
  A<-L<-ID<-Y<-Z<-Tv<-Int<-ALast<-LLast<-LFirst <- numeric()
  T0.vec<-T.vec<-Y.vec<-Z.vec <- rep(0, n)
  
  # Here are the coefficients determining the mediation and treatment assignment mechanisms
  bevec <- c(log(3/7), 2, log(0.5), log(1.5)) # Used to generate time-varying confounder L
  alvec <- c(log(2/7), 0.5, 0.5, log(4)) # Used to generate exposure (Intercept, L, LLast, ALast)
  
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
      A.vec<-L.vec<-ALast.vec<-LLast.vec<-LFirst.vec <- rep(0, N+1)
      # Implement Young's algorithm with multiple causes
      # Generate the survival time, then the cause
      while (muK*T0 > mu.tot & m <= N) {
        if (m == 0) {
          # First interval
          eta <- bevec[1] + bevec[2]*Ival + bevec[3]*0 + bevec[4]*0
          pval <- 1 / (1 + exp(-eta))
          L.vec[m+1] <- rbinom(1, 1, pval)
          
          eta <- alvec[1] + alvec[2]*L.vec[m+1] + alvec[3]*0 + alvec[4]*0
          pval <- 1 / (1 + exp(-eta))
          if (is.null(exposure)) {A.vec[m + 1] <- rbinom(1, 1, pval)}
          else {A.vec[m+1] <- exposure}
          
          ALast.vec[m+1] <- 0; LLast.vec[m+1] <- 0
          LFirst.vec <- rep(L.vec[m+1], N + 1)
        } else {
          # Subsequent intervals
          eta <- bevec[1] + bevec[2]*Ival + bevec[3]*A.vec[m] + bevec[4]*L.vec[m]
          pval <- 1 / (1 + exp(-eta))
          L.vec[m+1] <- rbinom(1, 1, pval) 
 
          eta <- alvec[1] + alvec[2]*L.vec[m+1] + alvec[3]*L.vec[m] + alvec[4]*A.vec[m] 
          pval <- 1 / (1 + exp(-eta))
          if (is.null(exposure)) {A.vec[m+1] <- rbinom(1, 1, pval)}
          else {A.vec[m+1] <- exposure}
          ALast.vec[m+1] <- A.vec[m]; LLast.vec[m+1] <- L.vec[m]
        }
        
        muval <- sum(exp(gamma.vec + A.vec[m+1]*psi.mat[ , m+1]))
        
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
        Z.vec[i] <- 0
      } else {
        # In the case of failure, use the ratio hazards to define the
        # relevant multinomial distribution on the K causes.
        Z.vec[i] <- sample(c(1:K), 1, prob = exp(gamma.vec + A.vec[m]*psi.mat[ ,m]))
      }
      
      # Store the outcomes
      T0.vec[i] <- T0
      T.vec[i] <- Tval
      Y.vec[i] <- m - 1
      ID <- c(ID, rep(i,m)) # Individual
      Int <- c(Int, c(1:m)) # Time point
      A <- c(A, A.vec[1:m]) # Time-updated treatment
      L <- c(L, L.vec[1:m]) # Time-updated covariate L
      ALast <- c(ALast, ALast.vec[1:m]) # Treatment at last t point
      LLast <- c(LLast, LLast.vec[1:m]) # Covariate L at last t point
      LFirst <- c(LFirst, LFirst.vec[1:m]) # Baseline covariate L value
      Z <- c(Z, rep(0,m - 1), Z.vec[i]) # Outcome: Z>0 indicates outcome of some type, determined by value)
      tv <- c(1:m); tv[m] <- Tval
      Tv <- c(Tv, tv) # If event occurs, exact time at which it occurred; o.w. equal to Int)
    }
    
    DeathsK.df <- data.frame(ID, Int, Tv, A, ALast, L, LLast, LFirst, Z)
    
    # Trim off the intervals beyond the Nth (loop goes one too far)
    DeathsK.df <- DeathsK.df[DeathsK.df$Int <= N, ]
    DeathsK.df$Int0 <- DeathsK.df$Int - 1
    
    return(DeathsK.df)
  }


  # Oracle: get counterfactuals for every simulated individual by setting exposure inside simulation
    # This gives us a measure of "optimal" performance given sample size and model
  set.seed(123+s)
  oracle1 <- simulation(exposure=1)
  set.seed(123+s) #Same seed to ensure the same
  oracle0 <- simulation(exposure=0)
  
  oracle <- data.table(rbind(oracle0,oracle1))
  
  fit <- summary(survfit(Surv(Int0, Tv, Z)  ~ A, data=oracle))
  surv <- data.frame(time = fit$time, 
                     surv = fit$surv,
                     exposure = fit$strata)
  sim.res$r0[1] <- 1 - min(surv$surv[surv$exposure=="A=0"])
  sim.res$r1[1] <- 1 - min(surv$surv[surv$exposure=="A=1"])
  sim.res$rd[1] <- sim.res$r1[1] - sim.res$r0[1]
  
  # Now create simulation to be used in analysis steps
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
    ps <- glm(A ~ ALast + L + LLast + as.factor(Int), family=binomial(link="logit"), data=boot)$fitted.values
    boot$denominator <- boot$A*ps + (1-boot$A)*(1-ps)
    
    # Numerator of weights
    ps <- glm(A ~ as.factor(Int), family=binomial(link="logit"), data=boot)$fitted.values
    boot$numerator <- boot$A*ps + (1-boot$A)*(1-ps)
    
    boot <- boot %>% 
      group_by(bid) %>%  
      mutate(wt=cumprod(numerator/denominator)) %>% 
      ungroup(bid) 

    # IP-weighted survival
    fit <- summary(survfit(Surv(Int0, Tv, Z)  ~ A, data=boot, weights=wt))
    surv <- data.frame(time = fit$time, 
                       surv = fit$surv,
                       exposure = fit$strata)
    boot.res$r0[1] <- 1 - min(surv$surv[surv$exposure=="A=0"])
    boot.res$r1[1] <- 1 - min(surv$surv[surv$exposure=="A=1"])
    boot.res$rd[1] <- boot.res$r1[1] - boot.res$r0[1]
  

##################################################################################################
## Classic g-computation
    # Time-ordering: LLast, ALast, L, A, Z

    # Model confounder
    mod.L <- glm(L ~ LLast + ALast + as.factor(Int), family=binomial(link="logit"), data=boot)

    # Model exposure (if you want natural course)
    # mod.A <- glm(A ~ ALast + L + LLast + as.factor(Int), family=binomial(link="logit"), data=boot)
    
    # Model outcome (flexsurv used bc survreg doesn't support start/stop coding)
    mod.D <- flexsurvreg(Surv(Int0,Tv,Z) ~ A + ALast + L + LLast, data=boot, dist="exp")
    
    # Take a Monte Carlo (MC) sample
    # Select first obs for each person to obtain joint empirical distribution of baseline covariates
    MC0 <- boot[boot$Int==1, (names(boot) %in% c("L", "A"))]
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
      Lp<-Ap<-Yp<-time <- numeric()
      time[1] <- j <- 1
      id <- d$id
  
      Lp[1] <- d$L
      if (is.null(exposure)) {
        Ap[1] <- d$A
      } else{
        Ap[1] <- exposure
      }
      
      # event status at first time point
      dYp <- data.table(A=Ap[1], ALast=0, L=Lp[1], LLast=0)
      expSim <- function (dat) {
        newD <- dat
        desX <- newD[,c("A", "ALast", "L","LLast")]
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
          ALast <- Ap[j-1]; LLast <- Lp[j-1]
          
          dLp <- data.table(LLast, ALast, Int=factor(j))
          Lp[j] <- pFunc(mod.L, dLp)
  
          if (is.null(exposure)) {
            dAp <- data.table(L=Lp[j], LLast, ALast, Int=factor(j))
            Ap[j] <- pFunc(mod.A, dAp)
          } else{
            Ap[j] <- exposure
          }
          
          dYp <- data.table(A=Ap[j], ALast, L=Lp[j], LLast)
          t0 <- expSim(dYp)
          if (t0<=0.001) {
            Yp[j-1] <- 1 #If the time interval is short enough to cause an error, place event at end of previous time point
            Lp <- Lp[1:(j-1)]
            Ap <- Ap[1:(j-1)]
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
      gdat <- data.table(id, time, jj, Ap, Lp, Yp)
      gdat$last <- as.numeric(gdat$Yp!=0 | gdat$time==lngth)
      return(gdat)
    }
  
    res0 <- lapply(1:montecarlo,function(x) {pgf(x ,mc_data=MC, length=N, exposure=0)})
      res0 <- do.call(rbind,res0)
    
    res1 <- lapply(1:montecarlo,function(x) {pgf(x, mc_data=MC, length=N, exposure=1)})
      res1 <- do.call(rbind, res1)
    
    gcomp.dat <- data.table(rbind(res1, res0))
    gcomp.dat$Int0 <- gcomp.dat$time - 1
    gcomp.dat$Int <- ifelse(gcomp.dat$last==1, gcomp.dat$jj, gcomp.dat$time)
    
    # Run outcome model
    fit <- summary(survfit(Surv(Int0, Int, Yp) ~ Ap, data=gcomp.dat))
    surv <- data.frame(time = fit$time, 
                       surv = fit$surv,
                       exposure = fit$strata)
    boot.res$r0[2] <- 1 - min(surv$surv[surv$exposure=="Ap=0"])
    boot.res$r1[2] <- 1 - min(surv$surv[surv$exposure=="Ap=1"])
    boot.res$rd[2] <- boot.res$r1[2] - boot.res$r0[2]

    
##################################################################################################
## G-computation via iterated conditional expectations (ICE)
    
    # Make confounder wide
    confounder <- boot %>% 
      select(bid, Int, L) %>% 
      pivot_wider(names_from=Int, names_prefix="L", values_from=L, values_fill=NA, names_sort=T) 
    
    # Make exposure wide
    exposure <- boot %>% 
      select(bid, Int, A) %>% 
      pivot_wider(names_from=Int, names_prefix="A", values_from=A, values_fill=NA, names_sort=T)  
    
    # Make outcome wide
    outcome <- boot %>% 
      select(bid, Int, Z) %>% 
      pivot_wider(names_from=Int, names_prefix="Z", values_from=Z, values_fill=0, names_sort=T)  
      # Once an event has occurred, all subsequent nodes must be 1
      for (i in 2:N) {
        outcome[ , i+1] <- outcome[ , i+1] + outcome[ , i]
      }
    
    # Interleave together in correct order (L, A, Z)
    wide <- select(boot, bid)
    for (i in 1:N){
      wide <- merge(wide, confounder[ , c(1, i+1)], by="bid")
      wide <- merge(wide, exposure[ , c(1, i+1)], by="bid")
      wide <- merge(wide, outcome[ , c(1, i+1)], by="bid")
    }
    
    # Set up call to ltmle package
    Anodes <- vars_select(names(wide), starts_with("A"))
    Lnodes <- vars_select(names(wide), starts_with("L"))
    Ynodes <- vars_select(names(wide), starts_with("Z"))
    
    Aformula <- rep(NA, N)
    for (i in 1:N){
      if (i==1) {
        #Formula for A(1)
        Aformula[i] <- paste(Anodes[i], "~", Lnodes[i], sep=" ")
      } else {
        #Formula for A(t), t>1
        Aformula[i] <- paste(Anodes[i], "~", Anodes[i-1], "+", Lnodes[i], "+", Lnodes[i-1], sep=" ")
      }
    }
    
    Yformula <- rep(NA, N)
    for (i in 1:N){
      if (i==1) {
        #Formula for Z(1)
        names(Yformula)[i] <- Ynodes[i]
        Yformula[i] <- paste("Q.kplus1 ~", Anodes[i], "+", Lnodes[i], sep=" ")
      } else {
        #Formula for Z(t), t>1
        names(Yformula)[i] <- Ynodes[i]
        Yformula[i] <- paste("Q.kplus1 ~", Anodes[i], "+", Anodes[i-1], "+", Lnodes[i], "+", Lnodes[i-1], sep=" ")
      }
    }

    # Use ltmle to implement ICE g-comp
    res <- ltmle(data=select(wide, -bid), 
                 Anodes=Anodes, Lnodes=Lnodes, Ynodes=Ynodes,
                 Qform=Yformula,
                 gform=Aformula,
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

