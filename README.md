# IPW-g-comp-sim

This code accompanies an In Press paper at the *American Journal of Epidemiology*, in which we compare using a plasmode simulation the performance of inverse probability weighting, Monte Carlo g-computation, and iterated conditional expectations g-computation when estimating the average treatment effect of a time-varying exposure on a survival outcome.

The main programs do the following:
  - plasmode_1_data.R -- Set up the observed data from the Effects of Aspirin in Gestation and Reproduction (EAGeR) trial.
  - plasmode_2_models.R -- Model the observed EAGeR data to obtain parameters for plasmode simulation
  - plasmode_3_truth_tvar.R -- Determine the true risk difference for the data generating mechanism
  - plasmode_4_analysis_tvar_cens.R -- Generate plasmode simulation and carry out analysis
  - plasmode_5_analysis_biased.R -- Supplementary analysis examining bias when only time-fixed confounding is accounted for
  - plasmode_5_analaysis_cont.R -- Supplementary analysis where time is treated as continuous
  - plasmode_6_results.R -- Read in, organize, and visualize results
 
 R programs included in the archive folder were used in earlier iterations of the paper, where we a non-plasmode simulation.
