# nb_bkmr
This repository contains R code for implementing the simulation study presented in Mutiso et al. (2024) Bayesian kernel machine regression for
count data: modelling the association between social vulnerability and COVID-19 deaths in South Carolina. JRSS-C

Comments: 

1) ordered_adjmat_sc.txt contains adjacency matrix used in the simulation

2) The 'SVI2018_US_COUNTY.csv' used in the simulation study is the CDC's SVI data for all US counties downloaded from the 
cdc website: https://www.atsdr.cdc.gov/placeandhealth/svi/data_documentation_download.html

3) The script bkmr_sim_data_nb2.R generates the simulation study data

4) The script bkmr_sim_initialize_nb2.R initializes different parameters

5) The script bkmr_sim_mcmc_nb2.R performs the mcmc estimation

6) The script bkmr_sim_runcode_nb2.R calls the scripts in 4, 5, and 6 to run the mcmc. 

7) The script bkmr_sim_pred_nb2.R makes predictions for new exposure profiles

8) The script bkmr_sim_heatmap_nb2.R calls the script in 8 to create heatmaps for new exposure profiles

9) The script bkmr_sim_pred_nb2_rr.R makes predictions for new exposure profiles used to create relative importances

10) The script bkmr_sim_relative_importance_nb2.R calls the script in 10 to create relative importances of various variables

11) Some notation in the programs deviates from the notation used in pap
