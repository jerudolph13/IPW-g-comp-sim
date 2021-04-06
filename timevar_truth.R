
###################################################################################################
#
# Purpose: Obtain true RDs for the time-varying analyses
#
# Authors: Jacqueline Rudolph, Ashley Naimi (Credit to Young and Moodie for DGM)
#
# Last Update: 02 Apr 2021
#
##################################################################################################

# Read in packages
lib <- "~/R/x86_64-pc-linux-gnu-library/4.0"
packages <- c("tidyverse", "data.table", "survival", "parallel")
for (package in packages) {
  library(package, character.only=T, lib.loc=lib)
}

# Define parameters and functions
n <- 1e6          #Sample size
N <- 5            #Number of time points
K <- 1            #Number of outcomes

expit <- function(x) {1/(1+exp(-x))}

# Prepare data set to hold the truth
truth <- data.frame(n.conf = NA,
                    lambda = NA,
                    rd = NA)

# Loop for parallelization
loop <- function(i) {
  
# 1 confounder ------------------------------------------------------------
  
  conf1 <- function(lambda){
    psi.mat <- matrix(0, nrow=K, ncol=N+1)
    psi.mat[1, ] <- log(2)
    
    gamma.vec <- rep(log(lambda/K))
    muK <- sum(exp(gamma.vec))
    A<-L<-ID<-Y<-Z<-Tv<-Int<-ALast<-LLast<-LFirst <- numeric()
    T0.vec<-T.vec<-Y.vec<-Z.vec <- rep(0, n)
    
    bevec <- c(log(3/7), 2, log(0.5), log(1.5)) 
    alvec <- c(log(2/7), 0.5, 0.5, log(4))
    cval <- 30
    
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
            ##Subsequent intervals
            #L affected by past L and past A
            eta <- bevec[1] + bevec[2]*Ival + bevec[3]*A.vec[m] + bevec[4]*L.vec[m]
            pval <- 1 / (1 + exp(-eta))
            L.vec[m+1] <- rbinom(1, 1, pval) 
            
            #A affected by L at this time point, last point and A at last time point
            eta <- alvec[1] + alvec[2]*L.vec[m+1] + alvec[3]*L.vec[m] + alvec[4]*A.vec[m] 
            pval <- 1 / (1 + exp(-eta))
            if (is.null(exposure)) {A.vec[m+1] <- rbinom(1, 1, pval)}
            else {A.vec[m+1] <- exposure}
            ALast.vec[m+1] <- A.vec[m]; LLast.vec[m+1] <- L.vec[m]
          }
          
          muval <- sum(exp(gamma.vec + A.vec[m+1]*psi.mat[ , m+1]))
          
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
          Z.vec[i] <- sample(c(1:K), 1, prob = exp(gamma.vec + A.vec[m]*psi.mat[ ,m]))
        }
        
        ##Store the outcomes
        T0.vec[i] <- T0
        T.vec[i] <- Tval
        Y.vec[i] <- m - 1
        ID <- c(ID, rep(i,m)) #Individual
        Int <- c(Int, c(1:m)) #Time point
        A <- c(A, A.vec[1:m]) #Time-updated treatment
        L <- c(L, L.vec[1:m]) #Time-updated covariate L
        ALast <- c(ALast, ALast.vec[1:m]) #Treatment at last t point
        LLast <- c(LLast, LLast.vec[1:m]) #Covariate L at last t point
        LFirst <- c(LFirst, LFirst.vec[1:m]) #Baseline covariate L value
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
    
    set.seed(123)
    dat1 <- simulation(exposure=1)
    dat1$last <- as.numeric(!duplicated(dat1$ID, fromLast=T))
    
    set.seed(123)
    dat0 <- simulation(exposure=0)
    dat0$last <- as.numeric(!duplicated(dat0$ID, fromLast=T))
    
    r0 <- mean(dat0$Z[dat0$last==1])
    r1 <- mean(dat1$Z[dat1$last==1])
    rd <- r1 - r0
    
    return(rd)
  }
  
  if (i==1) {
    truth$n.conf <- 1
    truth$lambda <- 0.01
    truth$rd <- conf1(lambda=0.01)
  }
  if (i==2) {
    truth$n.conf <- 1
    truth$lambda <- 0.05
    truth$rd <- conf1(lambda=0.05)
  }
  
# 3 confounders -----------------------------------------------------------
  
  conf3 <- function(lambda){
    psi.mat <- matrix(0, nrow=K, ncol=N+1)
    psi.mat[1, ] <- log(2)
    
    gamma.vec <- rep(log(lambda/K))
    muK <- sum(exp(gamma.vec))
    A<-J<-M<-L<-ID<-Y<-Z<-Tv<-Int<-ALast<-LLast<-LFirst<-JLast<-JFirst<-MLast<-MFirst <- numeric()
    T0.vec<-T.vec<-Y.vec<-Z.vec <- rep(0, n)
    
    bevec <- c(log(3/7), 2, log(0.5), log(1.5)) 
    bevec2 <- c(log(4/7), 2, log(0.6), log(1.4)) 
    bevec3 <- c(log(5/7), 2, log(0.7), log(1.3)) 
    alvec <- c(log(2/7), 0.5, 0.5, log(4), 0.6, 0.6, 0.8, 0.8) 
    cval <- 30
    
    simulation <- function (exposure) {
      
      for (i in 1:n) {
        T0 <- rexp(1, lambda) #Generate T0 from an exponential dist with constant rate=lamba
        Ival <- as.numeric(T0 < cval)
        ##Begin the interval-by-interval simulation
        m <- 0
        mu.tot <- 0
        A.vec<-J.vec<-M.vec<-L.vec<-ALast.vec<-JLast.vec<-MLast.vec<-LLast.vec<-JFirst.vec<-MFirst.vec<-LFirst.vec <- rep(0, N+1)
        ##Implement Young's algorithm with multiple causes
        ##Generate the survival time, then the cause
        while (muK*T0 > mu.tot & m <= N) {
          if (m == 0) {
            ##First interval
            eta <- bevec[1] + bevec[2]*Ival + bevec[3]*0 + bevec[4]*0
            pval <- 1 / (1 + exp(-eta))
            L.vec[m+1] <- rbinom(1, 1, pval)
            eta <- bevec2[1] + bevec2[2]*Ival + bevec2[3]*0 + bevec2[4]*0
            pval <- 1 / (1 + exp(-eta))
            J.vec[m+1] <- rbinom(1, 1, pval)
            eta <- bevec3[1] + bevec3[2]*Ival + bevec3[3]*0 + bevec3[4]*0
            pval <- 1 / (1 + exp(-eta))
            M.vec[m+1] <- rbinom(1, 1, pval)
            
            eta <- alvec[1] + alvec[2]*L.vec[m+1] + alvec[3]*0 + alvec[4]*0 + alvec[5]*J.vec[m+1] + alvec[6]*0 + alvec[7]*M.vec[m+1] + alvec[8]*0
            pval <- 1 / (1 + exp(-eta))
            if (is.null(exposure)) {A.vec[m + 1] <- rbinom(1, 1, pval)}
            else {A.vec[m+1] <- exposure}
            
            ALast.vec[m+1] <- 0; LLast.vec[m+1] <- 0; JLast.vec[m+1] <- 0; MLast.vec[m+1] <- 0
            LFirst.vec <- rep(L.vec[m+1], N + 1)
            JFirst.vec <- rep(J.vec[m+1], N + 1)
            MFirst.vec <- rep(M.vec[m+1], N + 1)
            
          } else {
            ##Subsequent intervals
            eta <- bevec[1] + bevec[2]*Ival + bevec[3]*A.vec[m] + bevec[4]*L.vec[m]
            pval <- 1 / (1 + exp(-eta))
            L.vec[m+1] <- rbinom(1, 1, pval) 
            eta <- bevec2[1] + bevec2[2]*Ival + bevec2[3]*A.vec[m] + bevec2[4]*J.vec[m]
            pval <- 1 / (1 + exp(-eta))
            J.vec[m+1] <- rbinom(1, 1, pval) 
            eta <- bevec3[1] + bevec3[2]*Ival + bevec3[3]*A.vec[m] + bevec3[4]*M.vec[m]
            pval <- 1 / (1 + exp(-eta))
            M.vec[m+1] <- rbinom(1, 1, pval) 
            
            eta <- alvec[1] + alvec[2]*L.vec[m+1] + alvec[3]*L.vec[m] + alvec[4]*A.vec[m] +alvec[5]*J.vec[m+1] + alvec[6]*J.vec[m] 
            + alvec[7]*M.vec[m+1] + alvec[8]*M.vec[m]
            pval <- 1 / (1 + exp(-eta)) 
            if (is.null(exposure)) {A.vec[m+1] <- rbinom(1, 1, pval)}
            else {A.vec[m+1] <- exposure}
            ALast.vec[m+1] <- A.vec[m]; LLast.vec[m+1] <- L.vec[m]; JLast.vec[m+1] <- J.vec[m]; MLast.vec[m+1] <- M.vec[m]
          }
          
          muval <- sum(exp(gamma.vec + A.vec[m+1]*psi.mat[ , m+1]))
          
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
          Z.vec[i] <- sample(c(1:K), 1, prob = exp(gamma.vec + A.vec[m]*psi.mat[ ,m]))
        }
        
        ##Store the outcomes
        T0.vec[i] <- T0
        T.vec[i] <- Tval
        Y.vec[i] <- m - 1
        ID <- c(ID, rep(i,m)) #Individual
        Int <- c(Int, c(1:m)) #Time point
        A <- c(A, A.vec[1:m]) #Time-updated treatment
        L <- c(L, L.vec[1:m]) #Time-updated covariate L
        J <- c(J, J.vec[1:m]) #Time-updated covariate J
        M <- c(M, M.vec[1:m]) #Time-updated covariate M
        ALast <- c(ALast, ALast.vec[1:m]) #Treatment at last t point
        LLast <- c(LLast, LLast.vec[1:m]) #Covariate L at last t point
        JLast <- c(JLast, JLast.vec[1:m]) #Covariate J at last t point
        MLast <- c(MLast, MLast.vec[1:m]) #Covariate M at last t point
        LFirst <- c(LFirst, LFirst.vec[1:m]) #Baseline covariate L value
        JFirst <- c(JFirst, JFirst.vec[1:m]) #Baseline covariate J value
        MFirst <- c(MFirst, MFirst.vec[1:m]) #Baseline covariate M value
        Z <- c(Z, rep(0,m - 1), Z.vec[i]) #Outcome: Z>0 indicates outcome of some type, determined by value)
        tv <- c(1:m); tv[m] <- Tval
        Tv <- c(Tv, tv) #If event occurs, exact time at which it occurred; o.w. equal to Int)
      }
      
      DeathsK.df <- data.frame(ID, Int, Tv, A, ALast, L, LLast, LFirst, J, JLast, JFirst, M, MLast, MFirst, Z)
      
      ##Trim off the intervals beyond the Nth (loop goes one too far)
      DeathsK.df <- DeathsK.df[DeathsK.df$Int <= N, ]
      DeathsK.df$Int0 <- DeathsK.df$Int - 1
      
      return(DeathsK.df)
    }
    
    set.seed(123)
    dat1 <- simulation(exposure=1)
    dat1$last <- as.numeric(!duplicated(dat1$ID, fromLast=T))
    
    set.seed(123)
    dat0 <- simulation(exposure=0)
    dat0$last <- as.numeric(!duplicated(dat0$ID, fromLast=T))
    
    r0 <- mean(dat0$Z[dat0$last==1])
    r1 <- mean(dat1$Z[dat1$last==1])
    rd <- r1 - r0
    
    return(rd)
  }
  
  if (i==3) {
    truth$n.conf <- 3
    truth$lambda <- 0.01
    truth$rd <- conf3(lambda=0.01)
  }
  if (i==4) {
    truth$n.conf <- 3
    truth$lambda <- 0.05
    truth$rd <- conf3(lambda=0.05)
  }
  return(truth)
}

cores <- detectCores() - 1
truth <- mclapply(1:4, function(ii) loop(ii), mc.cores=cores, mc.set.seed=FALSE)
  truth <- do.call(rbind, truth)

write.table(truth, file="../results/timevar_truth.txt", sep="\t")

