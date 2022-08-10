
###################################################################################################################
#
# Purpose: Summarize and visualize the results from the time-varying IPW/g-comp simulation
#
# Author: Jacqueline Rudolph
#
# Last Update: 21 Feb 2022
#
###################################################################################################################


library("tidyverse")

truth <- read_csv("../results/tvar_truth.csv")
  true_rd <- truth$rd


# Helper functions --------------------------------------------------------

read.res <- function(i) {
  
  dat <- read_csv(file=paste0("../results/tvar_cens_n-",n,"_mc-",mc,"_",i,".csv")) %>% 
    filter(!grepl('Error', sim)) %>% 
    mutate(rd = as.numeric(rd),
           se = as.numeric(se),
           sim = as.numeric(sim))
  res <- dat %>% 
    mutate(truth = true_rd,
           coverage = ((rd - 1.96*se) <= truth) & (truth <= (rd + 1.96*se)))
  
  return(res)
}
  
read.res2 <- function(i) {
  
  dat <- read_csv(file=paste0("../results/tvar_cens_bias", i, ".csv")) %>% 
    filter(!grepl('Error', sim)) %>% 
    mutate(rd = as.numeric(rd),
           se = as.numeric(se),
           sim = as.numeric(sim))
  res <- dat %>% 
    mutate(truth = true_rd,
           coverage = ((rd - 1.96*se) <= truth) & (truth <= (rd + 1.96*se)))
    
  return(res)
}

read.res3 <- function(i) {
  
  dat <- read_csv(file=paste0("../results/tvar_cens_cont", i, ".csv")) %>% 
    filter(!grepl('Error', sim)) %>% 
    mutate(rd = as.numeric(rd),
           se_rd = as.numeric(se_rd),
           loghr = as.numeric(loghr),
           se_hr = as.numeric(se_hr),
           sim = as.numeric(sim))
  res <- dat %>% 
    mutate(truth = true_rd,
           cov_rd = ((rd - 1.96*se_rd) <= truth) & (truth <= (rd + 1.96*se_rd)),
           cov_hr = ((loghr - 1.96*se_hr) <= log(2)) & (log(2) <= (loghr + 1.96*se_hr)))
  
  return(res)
}

summ.sim <- function (data) {
  
  summ <- data %>% 
    group_by(method) %>% 
    summarise(avg_rd=mean(rd), sd_rd=sd(rd), avg_se=mean(se), cov=mean(coverage), n=n()) %>% 
    mutate(bias = avg_rd - true_rd, ser = avg_se/sd_rd) %>% 
    select(method, avg_rd, bias, sd_rd, avg_se, ser, cov, n)
  
  return(summ)
}


# Get results -------------------------------------------------------------

n <- 2500
mc <- 5000

all.res <- lapply(1:5, function(tt) {read.res(tt)})
all.res <- do.call(rbind, all.res) 

summ <- summ.sim(all.res)
view(summ)


# Supplementary analyses --------------------------------------------------

# Don't control for time-varying confounders
all.res <- lapply(1:5, function(tt) {read.res2(tt)})
all.res <- do.call(rbind, all.res) 

summ2 <- summ.sim(all.res)
view(summ2)

# Use continuous time
all.res <- lapply(1:5, function(tt) {read.res3(tt)})
all.res <- do.call(rbind, all.res) 

summ3 <- all.res %>% 
  group_by(method) %>% 
  summarise(avg_rd=mean(rd), sd_rd=sd(rd), avg_se_rd=mean(se_rd), cov_rd=mean(cov_rd), 
            avg_hr=mean(loghr), sd_hr=sd(loghr), avg_se_hr=mean(se_hr), cov_hr=mean(cov_hr), n=n()) %>% 
  mutate(bias = avg_rd - true_rd, ser = avg_se_rd/sd_rd)
view(summ3)


# Visualize results -------------------------------------------------------

thm <- theme_classic() +
  theme(
    #Format axes
    axis.title = element_text(family="Helvetica", size=14, color="black"),
    axis.text = element_text(family="Helvetica", size=14, color="black"),
    axis.line = element_line(size=0.75),
    axis.ticks = element_line(size=0.75),
    
    #Format legend
    legend.text = element_text(family="Helvetica", size=14, color="black",
                               margin=margin(t=0.25,b=0.25, unit="lines")),
    legend.title = element_text(family="Helvetica", size=14, color="black"),
    legend.title.align = 0.5,
    legend.position = c(0.85, 0.9),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.direction="vertical",
    
    #Add space around plot
    plot.margin = unit(c(2, 2, 2, 2), "lines")
  )

# Comparing results of 3 approaches
fig <- all.res %>% 
  rename(Estimator=method, Estimate=rd) %>% 
  filter(Estimator!="MC2" & Estimator!="Crude")
fig$Estimator <- recode(fig$Estimator, 
                        "MC1"="MC g-computation",
                        "ICE"="ICE g-computation")
fig$Estimator <- factor(fig$Estimator, 
                        levels=c("IPW", "MC g-computation", "ICE g-computation"))

# Freqpoly histogram
jpeg("../figures/res_n-1226_mc-1226.jpg", height=6, width=7, units="in", res=300)
ggplot(data=fig, aes(x=Estimate, linetype=Estimator)) + thm +
  labs(x="\nRisk difference", y="Number of simulations\n") +
  geom_vline(xintercept=true_rd, linetype="solid", color="gray", size=0.75) +
  geom_freqpoly(color="black", binwidth=0.005, size=0.75) +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0))
dev.off()

