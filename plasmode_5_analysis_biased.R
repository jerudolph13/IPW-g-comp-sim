
###############################################################################################
#
# Project: Time-varying plasmode simulation from EAGeR data (bias due to t-var confounding)
#
# Author: Jacqueline Rudolph
#
# Last Update: 31 Jan 2022
#
###############################################################################################


lib <- "~/R/x86_64-pc-linux-gnu-library/4.0"
packages <- c("dplyr", "magrittr", "readr", "broom", "tidyr", "data.table", "tidyselect", "survival", 
              "flexsurv", "ltmle", "parallel")
for (package in packages) {
  library(package, character.only=T, lib.loc=lib)
}

# Define parameters and functions
start_sim <- 1		      
end_sim <- 100
nboot <- 500	      # Number of bootstrap resamples
n <- 1226		        # Sample size
n_mc <- 1226        # Size of Monte Carlo resample
k <- 1              # Number of outcomes
N <- 10             # Number of time points
beta.lagx <- log(3) # Effect of past exposure
beta.lagc <- log(3) # Effect of past aspirin compliance
beta.lagn <- log(3) # Effect of past nausea

expit <- function(x) {1/(1+exp(-x))}

# Prepare data set to hold simulation results
sim.res <- data.frame(
  method=c("Crude", "IPW", "MC1", "MC2", "ICE"),
  r1=rep(NA, 5),
  r0=rep(NA, 5),
  rd=rep(NA, 5),
  sim=rep(NA, 5),
  stringsAsFactors = FALSE
)


# Read in data ------------------------------------------------------------

base <- read_csv(file="../data/eager_base_limited.csv")

tvar.param <- read_csv(file="../results/tvar_param.csv")
  coef.x <- tvar.param$coef.x
  coef.c <- tvar.param$coef.c
  coef.n <- tvar.param$coef.n
  coef.y <- tvar.param$coef.y
  coef.d <- tvar.param$coef.d


# Resample data -----------------------------------------------------------

sim_rep <- function(iter) {
  sim.res$sim <- rep(iter, 5)
  set.seed(iter)

index <- sample(1:nrow(base), n, replace=T)
resample <- base[index, ]
resample$id <- 1:n


# Moodie DGM --------------------------------------------------------------

# Here the (untreated) rate of the outcome is set to lambda
lambda <- 0.05
gamma.vec <- log(lambda)
muK <- exp(gamma.vec)

# Here the (untreated) rate of censoring is set to lambda_d
lambda_d <- 0.1
gamma.vec_d <- log(lambda_d)
muK_d <- exp(gamma.vec_d)

# Baseline variables
id<-smoke<-white<-age<-BMI <- numeric()
# Time-varying covariates (Z1=aspirin; Z2=nausea)
X<-XLast<-Z1<-Z1Last<-Z2<-Z2Last <- numeric()
# Outcomes
K<-Y<-D<-Tv<-Int <- numeric()

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
    
    # Generate the (untreated) time to censoring
    T0_d <- rexp(1, lambda_d)
    
    # Begin the interval-by-interval simulation
    m <- 0
    mu.tot<-mu.tot_d <- 0
    X.vec<-XLast.vec<-Z1.vec<-Z1Last.vec<-Z2.vec<-Z2Last.vec<-D.vec <- rep(0, N+1)
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
                      coef.x[4]*age.vec[m+1] + coef.x[5]*BMI.vec[m+1] +
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
                        coef.c[5]*age.vec[m+1] + coef.c[6]*BMI.vec[m+1])
        Z1.vec[m+1] <- rbinom(1, 1, pval)
        
        # Generate nausea
        pval <- expit(coef.n[1] + beta.lagn*Z2Last.vec[m+1] + coef.n[2]*XLast.vec[m+1] + beta.ival*Ival +
                        coef.n[3]*white.vec[m+1] + coef.n[4]*smoke.vec[m+1] +
                        coef.n[5]*age.vec[m+1] + coef.n[6]*BMI.vec[m+1])
        Z2.vec[m+1] <- rbinom(1, 1, pval)
        
        # Generate exposure
        pval <- expit(coef.x[1] + beta.lagx*XLast.vec[m+1] + 
                        coef.x[2]*white.vec[m+1] + coef.x[3]*smoke.vec[m+1] +
                        coef.x[4]*age.vec[m+1] + coef.x[5]*BMI.vec[m+1] +
                        coef.x[6]*Z1.vec[m+1] + coef.x[6]*Z1Last.vec[m+1] +
                        coef.x[7]*Z2.vec[m+1] + coef.x[7]*Z2Last.vec[m+1])
        if (is.null(exposure)) {X.vec[m+1] <- rbinom(1, 1, pval)} 
        else {X.vec[m+1] <- exposure}
      }
      
      # Include baseline confounders in this equation
      muval <- sum(exp(gamma.vec + X.vec[m+1]*coef.y[2] - white.vec[m+1]*coef.y[3] - smoke.vec[m+1]*coef.y[4] - 
                         age.vec[m+1]*coef.y[5] - BMI.vec[m+1]*coef.y[6]))
      
      muval_d <- sum(exp(gamma.vec_d - white.vec[m+1]*coef.d[2] - smoke.vec[m+1]*coef.d[3] - 
                         age.vec[m+1]*coef.d[4] - BMI.vec[m+1]*coef.d[5]))
      
      # Tval is computed for each interval, but is overwritten until the final interval
      Tval <- m + (muK*T0 - mu.tot)/muval
      Tval_d <- m + (muK_d*T0_d - mu.tot_d)/muval_d
      
      mu.tot <- mu.tot + muval
      mu.tot_d <- mu.tot_d + muval_d
      
      D.vec[m+1] <- ifelse(muK_d*T0_d > mu.tot_d, 0, 1)
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
    D <- c(D, D.vec[1:m]) # Time-updated censoring
    Y <- c(Y, rep(0, m-1), Y.vec[i]) # Outcome: Y>0 indicates outcome
    tv <- c(1:m); tv[m] <- Tval
    Tv <- c(Tv, tv) # If event occurs, exact time at which it occurred; o.w. equal to Int)
  }
  
  DeathsK.df <- data.frame(id, Int, Tv, X, XLast, white, smoke, age, BMI, Z1, Z1Last,
                           Z2, Z2Last, Y, D)
  
  # Trim off the intervals beyond the Nth (loop goes one too far)
  DeathsK.df <- DeathsK.df[DeathsK.df$Int <= N, ]
  DeathsK.df$Int0 <- DeathsK.df$Int - 1
  
  return(DeathsK.df)
}

sim.dat <- simulation(exposure=NULL)

sim.dat_censored <- sim.dat %>% 
  group_by(id) %>% 
  mutate(cum.D = cumsum(D)) %>% 
  ungroup() %>% 
  filter(cum.D<=1) %>% 
  mutate(D = ifelse(Y==1, 0, D)) %>% 
  select(-cum.D)


# Bootstrap ---------------------------------------------------------------

boot.res <- data.frame(
  method = c("Crude", "IPW", "MC1", "MC2", "ICE"),
  boot_num=rep(NA, 5),
  r0=rep(NA, 5),
  r1=rep(NA, 5),
  rd=rep(NA, 5),
  stringsAsFactors = FALSE
)

boot_rep <- function(r) {
  set.seed(r+1)
  boot.res$boot_num <- rep(r, 5)

# Sample with replacement
firstobs <- sim.dat_censored[sim.dat_censored$Int == 1, ]
samp <- table(firstobs[sample(1:nrow(firstobs),nrow(firstobs),replace=T), (names(sim.dat_censored) == "id")])

boot <- NULL
if (r==0) {
  boot <- sim.dat_censored %>% 
    rename(bid = id)
} else {
  for(zzz in 1:max(samp)){ 
    cc <- sim.dat_censored[sim.dat_censored$id %in% names(samp[samp %in% c(zzz:max(samp))]),]
    cc$bid <- paste0(cc$id, zzz)
    boot <- rbind(boot, cc)
  }
  boot <- select(boot, -id)
}
  

# Crude -------------------------------------------------------------------

fit <- summary(survfit(Surv(Int0, Int, Y)  ~ X, data=boot))
surv <- data.frame(time = fit$time, 
                   surv = fit$surv,
                   exposure = fit$strata)
boot.res$r0[1] <- 1 - min(surv$surv[surv$exposure=="X=0"])
boot.res$r1[1] <- 1 - min(surv$surv[surv$exposure=="X=1"])
boot.res$rd[1] <- boot.res$r1[1] - boot.res$r0[1]


# IPW ---------------------------------------------------------------------

# Exposure weights
ps <- glm(X ~ XLast + white + smoke + age + BMI + as.factor(Int), 
          family=binomial(link="logit"), data=boot)$fitted.values
denominator <- boot$X*ps + (1-boot$X)*(1-ps)

ps <- glm(X ~ as.factor(Int), family=binomial(link="logit"), data=boot)$fitted.values
numerator <- boot$X*ps + (1-boot$X)*(1-ps)
wt_x <- numerator/denominator

# Censoring weights
den_d <- glm((D==0) ~ X + white + smoke + age + BMI + as.factor(Int), 
          family=binomial(link="logit"), data=boot)$fitted.values
num_d <- glm((D==0) ~ X + as.factor(Int), family=binomial(link="logit"), data=boot)$fitted.values
wt_d <- ifelse(boot$D==1, 0, num_d/den_d)

boot$wt <- wt_x*wt_d
boot <- boot %>% 
  group_by(bid) %>%  
  mutate(cum.wt=cumprod(wt)) %>% 
  ungroup(bid) 

# IP-weighted survival
fit <- summary(survfit(Surv(Int0, Int, Y)  ~ X, data=boot, weights=cum.wt))
surv <- data.frame(time = fit$time, 
                   surv = fit$surv,
                   exposure = fit$strata)
boot.res$r0[2] <- 1 - min(surv$surv[surv$exposure=="X=0"])
boot.res$r1[2] <- 1 - min(surv$surv[surv$exposure=="X=1"])
boot.res$rd[2] <- boot.res$r1[2] - boot.res$r0[2]


# MC g-computation --------------------------------------------------------

# Model outcome (flexsurv used bc survreg doesn't support start/stop coding)
mod.Y <- flexsurvreg(Surv(Int0, Int, Y) ~ X + XLast + white + smoke + age + BMI, 
                     data=boot, dist="exp")

# Take a Monte Carlo (MC) sample
# Select first obs for each person to obtain joint empirical distribution of baseline covariates
MC0 <- boot %>% filter(Int==1) %>% select(X, Z1, Z2, white, smoke, age, BMI)
index <- sample(1:nrow(MC0), n_mc, replace=T)
MC <- MC0[index, ]
MC$id <- 1:n_mc

# Predict follow-up based on g-formula using PGF function
  # Intervention on exposure
  # Intervention to remove censoring
pgf <- function(ii, mc_data, length, exposure=NULL){

  pFunc <- function(mod,ndat){as.numeric(predict(mod, newdata=ndat, type="response") > runif(1))}
  
  expSim <- function (dat) {
    newD <- dat
    desX <- newD[, c("X", "XLast", "white", "smoke", "age", "BMI")]
    p_y <- exp(coef(mod.Y)[names(coef(mod.Y))=="rate"])*
      exp(coef(mod.Y)[!names(coef(mod.Y))=="rate"]%*%t(desX))
    y0 <- rexp(1, p_y)
    return(list(p_y, y0))
  }
  
  d <- mc_data
  d <- d[d$id==ii, ]
  lngth <- length
  white<-smoke<-age<-BMI <- numeric()
  Xp<-Yp<-p_y<-time <- numeric()
  time[1] <- j <- 1
  id <- d$id
  
  white[1] <- d$white; smoke[1] <- d$smoke; age[1] <- d$age; BMI[1] <- d$BMI
  Xp[1] <- exposure
  
  # Event status at first time point
  dYp <- data.table(X=Xp[1], XLast=0, white=white[1], smoke=smoke[1], age=age[1], BMI=BMI[1])

  pred <- expSim(dYp)
  p_y[1] <- pred[[1]]
  t0 <- pred[[2]]
  if (t0<=1) {
    Yp[1] <- 1
  } else {
    Yp[1] <- 0
  }
  
  # subsequent time points
  for (j in 2:lngth) {

    white[j] <- d$white; smoke[j] <- d$smoke; age[j] <- d$age; BMI[j] <- d$BMI
    XLast <- Xp[j-1]
    
    Xp[j] <- exposure
      
    dYp <- data.table(X=Xp[j], XLast, white=white[j], smoke=smoke[j], age=age[j], BMI=BMI[j])
    pred <- expSim(dYp)
    p_y[j] <- pred[[1]]
    t0 <- pred[[2]]
    
    if (Yp[j-1]==1) {
      Yp[j] <- 1
    } else {
      if (t0<=0.001) {
        Yp[j-1] <- 1 #If the time interval is short enough to cause an error, place event at end of previous time point
        Yp[j] <- 1
      } else {
        if (t0>0.001 & t0<=1) {
          Yp[j] <- 1
        } else { 
          Yp[j] <- 0
        }
      }
    }
 
    time[j] <- j
  }
  gdat <- data.table(id, time, Xp, white, smoke, age, BMI, Yp, p_y)
  return(gdat)
}

set.seed(1)
res0 <- lapply(1:n_mc, function(x) {pgf(x, mc_data=MC, length=N, exposure=0)})
res0 <- do.call(rbind, res0)

set.seed(1)
res1 <- lapply(1:n_mc, function(x) {pgf(x, mc_data=MC, length=N, exposure=1)})
res1 <- do.call(rbind, res1)

# Estimate risk
  # Using outcome predicted in pgf function
  res0_trim <- res0 %>% 
    group_by(id) %>% 
    mutate(cumYp = cumsum(Yp)) %>% 
    filter(cumYp<=1) %>% 
    select(-cumYp)
  fit <- summary(survfit(Surv(time-1, time, Yp)  ~ 1, data=res0_trim))
  boot.res$r0[3] <- 1 - min(fit$surv)
  
  res1_trim <- res1 %>% 
    group_by(id) %>% 
    mutate(cumYp = cumsum(Yp)) %>% 
    filter(cumYp<=1) %>% 
    select(-cumYp)
  fit <- summary(survfit(Surv(time-1, time, Yp)  ~ 1, data=res1_trim))
  boot.res$r1[3] <- 1 - min(fit$surv)
  
  boot.res$rd[3] <- boot.res$r1[3] - boot.res$r0[3]

  # Using predicted hazards of event
  res0_summ <- res0 %>% 
    group_by(time) %>% 
    summarize(avg_p_y = mean(p_y)) %>% 
    ungroup() %>% 
    mutate(risk = 1 - cumprod(1 - avg_p_y))

  res1_summ <- res1 %>% 
    group_by(time) %>% 
    summarize(avg_p_y = mean(p_y)) %>% 
    ungroup() %>% 
    mutate(risk = 1 - cumprod(1 - avg_p_y))
  
  boot.res$r0[4] <- res0_summ$risk[N]
  boot.res$r1[4] <- res1_summ$risk[N]
  boot.res$rd[4] <- boot.res$r1[4] - boot.res$r0[4]
  
  
# ICE g-computation -------------------------------------------------------

long <- boot %>% expand(bid, Int)
# Use last record to get the time at which outcome or censoring occurred
last <- boot %>% mutate(last = as.numeric(!duplicated(bid, fromLast=T))) %>% 
  filter(last==1) %>% 
  select(bid, Int, white, smoke, age, BMI) %>% 
  rename(time=Int)
long2 <- left_join(long, last, by=c("bid"))
long3 <- left_join(long2, select(boot, bid, Int, X, Y, D), by=c("bid", "Int")) %>% 
  group_by(bid) %>% 
  # Once someone has event, Y must be 1 for rest of time points
  # Once someone is censored, Y must be NA for rest of time points
  mutate(Y = ifelse(Int>time, 0, Y),
         Y = cumsum(Y),
         Y = ifelse(Y==0 & is.na(D), NA, Y),
         D = BinaryToCensoring(is.censored=D)) %>% 
  ungroup()

# Make exposure wide
exposure <- long3 %>% 
  select(bid, Int, X) %>% 
  pivot_wider(names_from=Int, names_prefix="X", values_from=X)  

# Make censoring wide
censor <- long3 %>% 
  select(bid, Int, D) %>% 
  pivot_wider(names_from=Int, names_prefix="D", values_from=D)

# Make outcome wide
outcome <- long3 %>% 
  select(bid, Int, Y) %>% 
  pivot_wider(names_from=Int, names_prefix="Y", values_from=Y)  

# Interleave together in correct order (baseline, Z1, Z2, X, D, Y)
wide <- long3 %>% filter(Int==1) %>% select(bid, white, smoke, age, BMI)
for (i in 1:N){
  wide <- merge(wide, exposure[ , c(1, i+1)], by="bid")
  wide <- merge(wide, censor[ , c(1, i+1)], by="bid")
  wide <- merge(wide, outcome[ , c(1, i+1)], by="bid")
}

# Set up call to ltmle package
Anodes <- vars_select(names(wide), starts_with("X"))
Ynodes <- vars_select(names(wide), starts_with("Y"))
Cnodes <- vars_select(names(wide), starts_with("D"))

# Use ltmle to implement ICE g-comp
res <- ltmle(data=select(wide, -bid), 
             Anodes=Anodes, Ynodes=Ynodes, Cnodes=Cnodes,
             abar=list(treatment=rep(1, N), control=rep(0, N)),
             survivalOutcome=T,
             SL.library=NULL,
             gcomp=T)

summ <- summary(res)  
boot.res$r0[5] <- summ$effect.measures$control$estimate
boot.res$r1[5] <- summ$effect.measures$treatment$estimate
boot.res$rd[5] <- summ$effect.measures$ATE$estimate

return(boot.res)
}


# Aggregate results -------------------------------------------------------

all.boot <- lapply(0:nboot, function(tt) {boot_rep(tt)})
all.boot <- do.call(rbind, all.boot)

# For point estimates, pull out results where boot=0
boot0 <- filter(all.boot, boot_num == 0)
sim.res$r0 <- boot0$r0
sim.res$r1 <- boot0$r1
sim.res$rd <- boot0$rd
all.boot <- filter(all.boot, boot_num>0)

# Summarize over bootstraps
boot.summ <- all.boot %>%
  group_by(method) %>%
  summarize(se = sd(rd))

sim.res <- left_join(sim.res, boot.summ, by="method")

return(sim.res)
}

cores <- detectCores() - 2
all.res <- mclapply(1:nsim, function(ii) sim_rep(ii), mc.cores=cores, mc.set.seed=FALSE)
#Use lapply when testing code for errors
#all.res <- lapply(1:nsim, function(ii) sim_rep(ii))
all.res <- do.call(rbind,all.res)


# Output results ----------------------------------------------------------

filename <- paste("../results/tvar_cens_bias", end_sim/100, ".csv", sep="")
write_csv(all.res, file=filename)

