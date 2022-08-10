
###############################################################################################
#
# Project: Models to inform parameters for EAGeR plasmode simulation
#
# Author: Jacqueline Rudolph
#
# Last Update: 20 Dec 2021
#
###############################################################################################


packages <- c("tidyverse", "broom", "survival")
for (package in packages) {
  library(package, character.only=T)
}

# Read in data ------------------------------------------------------------

eager_base <- read_csv(file="../data/eager_base.csv") %>% 
  mutate(exposure = as.numeric(exercise>1))


# Time-fixed parameters ---------------------------------------------------

# Exposure affected by confounders
mod.x <- glm(exposure ~ white + smoke + age + BMI, family=binomial(link="logit"), data=eager_base)
coef.x <- tidy(mod.x)$estimate*2 # Increase magnitude of confounding
coef.x[1] <- coef.x[1]/2         # Leave reference probability the same
coef.x[5] <- coef.x[5]*10        # Make BMI effect even larger

# Outcome affected by exposure and confounders
mod.y <- survreg(Surv(week, conception) ~ exposure + white + smoke + age + BMI, 
                 dist="exponential", data=eager_base)
coef.y <- tidy(mod.y)$estimate*2 # Increase magnitude of confounding
coef.y[1] <- coef.y[1]/2         # Leave reference rate the same
coef.y[2] <- log(2)               # True effect of X on Y

# Censoring affected by confounders
mod.d <- survreg(Surv(week, drop) ~ white + smoke + age + BMI, dist="exponential", data=eager_base)
coef.d <- tidy(mod.d)$estimate     # Reverse direction of selection bias (?)
#coef.d[1] <- coef.d[1]*-1         # Leave reference rate the same

# Store parameters
tfix.param <- data.frame(coef.x = c(coef.x, NA),
                         coef.y = coef.y,
                         coef.d = c(coef.d, NA))


# Time-varying parameters -------------------------------------------------

# Exposure (moderate to high exercise per week)
mod.x <- glm(exposure ~ white + smoke + age + BMI + compliance + nausea, 
             family=binomial(link="logit"), data=eager_base)
coef.x <- tidy(mod.x)$estimate*2 # Increase magnitude of confounding
coef.x[1] <- coef.x[1]/2         # Leave reference probability the same
coef.x[5] <- coef.x[5]*10        # Make BMI effect even larger

# Compliance affected by exposure and baseline confounders
mod.c <- glm(compliance ~ exposure + white + smoke + age + BMI,
             family=binomial(link="logit"), data=eager_base)
coef.c <- tidy(mod.c)$estimate

# Nausea affected by exposure and baseline confounders
mod.n <- glm(nausea ~ exposure + white + smoke + age + BMI,
             family=binomial(link="logit"), data=eager_base)
coef.n <- tidy(mod.n)$estimate

# Outcome affected by exposure and confounders
mod.y <- survreg(Surv(week, conception) ~ exposure + white + smoke + age + BMI + 
                   compliance + nausea, dist="exponential", data=eager_base)
coef.y <- tidy(mod.y)$estimate*2 # Increase magnitude of confounding
coef.y[1] <- coef.y[1]/2         # Leave reference rate the same
coef.y[2] <- log(2)              # True effect of X on Y

# Censoring informed by exposure and confounders
mod.d <- survreg(Surv(week, drop) ~ white + smoke + age + BMI + compliance + nausea, 
                 dist="exponential", data=eager_base)
coef.d <- tidy(mod.d)$estimate    # Reverse direction of selection bias (?)
# coef.d[1] <- coef.d[1]*-1         # Leave reference rate the same
# coef.d[6] <- coef.d[6]/-2         # Make compliance effect smaller

# Store parameters
tvar.param <- data.frame(coef.x = c(coef.x, NA),
                         coef.c = c(coef.c, rep(NA, 2)),
                         coef.n = c(coef.n, rep(NA, 2)),
                         coef.y = coef.y,
                         coef.d = c(coef.d, NA))


# Output results ----------------------------------------------------------

write_csv(tfix.param, file="../results/tfix_param.csv")
write_csv(tvar.param, file="../results/tvar_param.csv")

