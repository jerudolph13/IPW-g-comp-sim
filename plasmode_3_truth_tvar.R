
###############################################################################################
#
# Project: Time-varying plasmode simulation from EAGeR data (true RD)
#
# Author: Jacqueline Rudolph
#
# Last Update: 20 Dec 2021
#
###############################################################################################


lib <- "~/R/x86_64-pc-linux-gnu-library/4.0"
packages <- c("dplyr", "magrittr", "readr", "survival", "tidyselect")
for (package in packages) {
  library(package, character.only=T, lib.loc=lib)
}

# Define parameters and functions
n <- 1e6		        # Sample size
k <- 1              # Number of outcomes
N <- 10             # Number of time points
beta.lagx <- log(3) # Effect of past exposure
beta.lagc <- log(3) # Effect of past aspirin compliance
beta.lagn <- log(3) # Effect of past nausea

expit <- function(x) {1/(1+exp(-x))}


# Read in data ------------------------------------------------------------

base <- read_csv(file="../data/eager_base_limited.csv")

tvar.param <- read_csv(file="../results/tvar_param.csv")
  coef.x <- tvar.param$coef.x
  coef.c <- tvar.param$coef.c
  coef.n <- tvar.param$coef.n
  coef.y <- tvar.param$coef.y
  

# Resample data -----------------------------------------------------------

set.seed(123)
  
index <- sample(1:nrow(base), n, replace=T)
resample <- base[index, ]
resample$id <- 1:n


# Moodie DGM --------------------------------------------------------------

# Here the (untreated) rate of the outcome is set to lambda
lambda <- 0.05
gamma.vec <- log(lambda)
muK <- exp(gamma.vec)

# Baseline variables
id<-smoke<-white<-age<-BMI <- numeric()
# Time-varying covariates (Z1=aspirin; Z2=nausea)
X<-XLast<-Z1<-Z1Last<-Z2<-Z2Last <- numeric()
# Outcomes
K<-Y<-Tv<-Int <- numeric()

T0.vec<-T.vec<-K.vec<-Y.vec <- rep(0, n)

# Introduce confounding by creating indicator Ival, which affects Y and time-varying confounders
cval <- 30
beta.ival <- 2

# Begin the data-generation loop
simulation <- function (exposure) {
  
  for (i in 1:n) {
    
    # Generate the counterfactual (untreated) survival time
    T0 <- rexp(1, lambda) # Generate T0 from an exponential dist with constant rate=lambda
    Ival <- as.numeric(T0 < cval)
    
    # Begin the interval-by-interval simulation
    m <- 0
    mu.tot<-mu.tot_d <- 0
    X.vec<-XLast.vec<-Z1.vec<-Z1Last.vec<-Z2.vec<-Z2Last.vec <- rep(0, N+1)
    smoke.vec<-age.vec<-white.vec<-BMI.vec <- rep(0, N+1)
    
    # Implement Young's algorithm with multiple causes
    # Generate the survival time, then the cause
    while (muK*T0 > mu.tot & m <= N) {
      if (m == 0) {
        # First interval
        white.vec <- rep(resample$white[i], N+1)
        smoke.vec <- rep(resample$smoke[i], N+1)
        age.vec <- rep(resample$age[i], N+1)
        BMI.vec <- rep(resample$BMI[i], N+1)
        XLast.vec[m+1] <- 0
        Z1.vec[m+1] <- resample$compliance[i]
        Z1Last.vec[m+1] <- 0
        Z2.vec[m+1] <- resample$nausea[i]
        Z2Last.vec[m+1] <- 0
        
        # Generate exposure
        pval <- expit(coef.x[1] + beta.lagx*XLast.vec[m+1] + 
                        coef.x[2]*white.vec[m+1] + coef.x[3]*smoke.vec[m+1] +
                        coef.x[4]*age.vec[m+1] + coef.x[5]*log(BMI.vec[m+1]) +
                        coef.x[6]*Z1.vec[m+1] + coef.x[6]*Z1Last.vec[m+1] +
                        coef.x[7]*Z2.vec[m+1] + coef.x[7]*Z2Last.vec[m+1])
        
        if (is.null(exposure)) {X.vec[m+1] <- rbinom(1, 1, pval)} 
        else {X.vec[m+1] <- exposure}
        
      } else {
        # Subsequent intervals
        XLast.vec[m+1] <- X.vec[m]
        Z1Last.vec[m+1] <- Z1.vec[m]
        Z2Last.vec[m+1] <- Z2.vec[m]
        
        # Generate aspirin use
        pval <- expit(coef.c[1] + beta.lagc*Z1Last.vec[m+1] + coef.c[2]*XLast.vec[m+1] + beta.ival*Ival +
                        coef.c[3]*white.vec[m+1] + coef.c[4]*smoke.vec[m+1] +
                        coef.c[5]*age.vec[m+1] + coef.c[6]*log(BMI.vec[m+1]))
        Z1.vec[m+1] <- rbinom(1, 1, pval)
        
        # Generate nausea
        pval <- expit(coef.n[1] + beta.lagn*Z2Last.vec[m+1] + coef.n[2]*XLast.vec[m+1] + beta.ival*Ival +
                        coef.n[3]*white.vec[m+1] + coef.n[4]*smoke.vec[m+1] +
                        coef.n[5]*age.vec[m+1] + coef.n[6]*log(BMI.vec[m+1]))
        Z2.vec[m+1] <- rbinom(1, 1, pval)
        
        # Generate exposure
        pval <- expit(coef.x[1] + beta.lagx*XLast.vec[m+1] + 
                        coef.x[2]*white.vec[m+1] + coef.x[3]*smoke.vec[m+1] +
                        coef.x[4]*age.vec[m+1] + coef.x[5]*log(BMI.vec[m+1]) +
                        coef.x[6]*Z1.vec[m+1] + coef.x[6]*Z1Last.vec[m+1] +
                        coef.x[7]*Z2.vec[m+1] + coef.x[7]*Z2Last.vec[m+1])
        if (is.null(exposure)) {X.vec[m+1] <- rbinom(1, 1, pval)} 
        else {X.vec[m+1] <- exposure}
      }
      
      # Include baseline confounders in this equation
      muval <- sum(exp(gamma.vec + X.vec[m+1]*coef.y[2] - white.vec[m+1]*coef.y[3] - smoke.vec[m+1]*coef.y[4] - 
                         age.vec[m+1]*coef.y[5] - log(BMI.vec[m+1])*coef.y[6]))
      
      # Tval is computed for each interval, but is overwritten until the final interval
      Tval <- m + (muK*T0 - mu.tot)/muval
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
      # Otherwise outcome occurred
      Y.vec[i] <- 1
    }
    
    # Store the outcomes
    T0.vec[i] <- T0
    T.vec[i] <- Tval
    id <- c(id, rep(i,m)) # Individual
    Int <- c(Int, c(1:m)) # Time point
    X <- c(X, X.vec[1:m]) # Time-updated treatment
    XLast <- c(XLast, XLast.vec[1:m]) # Treatment at last t point
    white <- c(white, white.vec[1:m]) # Baseline covariate white
    smoke <- c(smoke, smoke.vec[1:m]) # Baseline covariate smoke
    age <- c(age, age.vec[1:m]) # Baseline covariate age
    BMI <- c(BMI, BMI.vec[1:m]) # Baseline covariate BMI
    Z1 <- c(Z1, Z1.vec[1:m]) # Time-updated covariate aspirin use
    Z1Last <- c(Z1Last, Z1Last.vec[1:m]) # Covariate aspirin at last t point
    Z2 <- c(Z2, Z2.vec[1:m]) # Time-updated covariate nausea
    Z2Last <- c(Z2Last, Z2Last.vec[1:m]) #Covariate nausea at last t point
    Y <- c(Y, rep(0, m-1), Y.vec[i]) # Outcome: Y>0 indicates outcome
    tv <- c(1:m); tv[m] <- Tval
    Tv <- c(Tv, tv) # If event occurs, exact time at which it occurred; o.w. equal to Int)
  }
  
  DeathsK.df <- data.frame(id, Int, Tv, X, XLast, white, smoke, age, BMI, Z1, Z1Last,
                           Z2, Z2Last, Y)
  
  # Trim off the intervals beyond the Nth (loop goes one too far)
  DeathsK.df <- DeathsK.df[DeathsK.df$Int <= N, ]
  DeathsK.df$Int0 <- DeathsK.df$Int - 1
  
  return(DeathsK.df)
}


# Generate simulation (X=1) -----------------------------------------------

# Reset seed for consistency with X=0 scenario
set.seed(123)

# Set exposure
sim.dat1 <- simulation(exposure=1) %>% 
  group_by(id) %>% 
  mutate(last = as.numeric(!duplicated(id, fromLast=T))) %>% 
  filter(last==1) %>% 
  select(-last)
  

# Generate simulation (X=0) -----------------------------------------------

# Reset seed for consistency with X=1 scenario
set.seed(123)

# Set exposure
sim.dat0 <- simulation(exposure=0) %>% 
  group_by(id) %>% 
  mutate(last = as.numeric(!duplicated(id, fromLast=T))) %>% 
  filter(last==1) %>% 
  select(-last)


# Estimate risk -----------------------------------------------------------

sim.dat <- bind_rows(sim.dat1, sim.dat0)

fit <- summary(survfit(Surv(Tv, (Y==1)) ~ X, data=sim.dat))
surv <- data.frame(time = fit$time, 
                   surv = fit$surv,
                   exposure = fit$strata)

r0 <- 1 - min(surv$surv[surv$exposure=="X=0"])
r1 <- 1 - min(surv$surv[surv$exposure=="X=1"])
rd <- r1 - r0


# Output results ----------------------------------------------------------

res <- tibble(method="truth", r0=r0, r1=r1, rd=rd)
write_csv(res, file="../results/tvar_truth.csv")
