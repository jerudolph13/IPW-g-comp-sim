
###############################################################################################
#
# Project: EAGeR data preparation
#
# Author: Jacqueline Rudolph
#
# Last Update: 06 Dec 2021
#
###############################################################################################


packages <- c("tidyverse", "broom", "survival", "tidyselect", "ltmle")
for (package in packages) {
  update.packages(package, ask=F)
  library(package, character.only=T)
}


# Read in data ------------------------------------------------------------

eager <- read.csv(file="../data/2020_11_06-eager_weekly.csv", header=TRUE) %>% 
  #Remove ID 417199 due to missing information
  filter(id != 417199) %>% 
  mutate(
    #baseline variables
    eligibility = as.numeric(eligibility == "new"),
    
    #t-varying compliance exposure
    compliance = as.numeric(studymed_b_imputed >= 5/7), 
    
    #t-varying confounders
    bleed = as.numeric(bleed_d_imputed > 0), 
    nausea = as.numeric(nausea_d_imputed > 0),
    
    #flag last record
    last = as.numeric(!duplicated(id, fromLast=T))) %>% 
  select(-c(bleed_d_imputed, nausea_d_imputed, studymed_b_imputed, conception))


# Manipulate data ---------------------------------------------------------

# Grab baseline data
eager_base <- filter(eager, week==1) %>% 
  select(-week)

# Define outcomes
eager_last <- eager %>% 
  filter(last==1) %>% 
  mutate(c_week = week - GA) %>% 
  select(id, c_week)
eager <- left_join(eager, eager_last, by="id") %>% 
  mutate(conception = ifelse(GA==99, 0, as.numeric((week - c_week) >= 0)))

# Define drop out
eager$drop <- ifelse(eager$last==1 & eager$outcome=="withdrawal", 1, 
                     ifelse(eager$last==1 & eager$conception==0 & eager$outcome %in% 
                              c("live birth", "pregnancy loss"), 1, 0))

# Subset to records through event or censoring
eager <- filter(eager, week <= 10)

# Remove extra records after conception
eager <- eager %>% 
  group_by(id) %>% 
  mutate(cum_concep = cumsum(conception)) %>% 
  filter(cum_concep <= 1) %>% 
  ungroup(id) %>% 
  select(-c(cum_concep, last))

# Grab last record
eager_last2 <- eager %>% 
  filter(as.numeric(!duplicated(id, fromLast=T))==1) %>% 
  select(id, week, conception, drop)

# Merge outcomes to baseline data
eager_base <- left_join(eager_base, eager_last2, by="id")


# Output data -------------------------------------------------------------

write_csv(eager_base, file="../data/eager_base.csv")

# Limited data set
base <- select(eager_base, age, white, smoke, BMI, compliance, nausea)

write_csv(base, file="../data/eager_base_limited.csv")