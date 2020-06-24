# - - - - - - - 
# set a working directory with the tutorial:
#  **/Tutorial_scripts/COVID-19/MultinomialApproximation
# change ** with your directory
setwd("~/Dropbox/PhD/Research/HMM-SEIR/NeurIPS/Tutorial_scripts/COVID-19/MultinomialApproximation")
dropbox_path = ""

# This file loads the dataset, the parameters and all the fuctions that are needed for the analysis
source("R/HMM-covid.R")

# - - - - - - - 
# Grid search
nn = 3000; dt = 1; 
# q_local_grid = seq(0.001,0.005,0.001); q_grid = seq(0.4,0.9,0.1)

# Run a grid search over the emission density parameters or use our results 
# likelihood_surface = MLE_check_2D(q_local_grid, q_grid, nn, dt, filename = 1)


# - - - - - - - 
# Run 100 particle filters
rep_plot <- 100;  cut_off=0; filename="1"

# Choose the best parameters on the grid
q_local = 0.00175 #likelihood_surface$q_local[which(likelihood_surface$lik==max(likelihood_surface$lik))]
q = 0.8 #likelihood_surface$q[which(likelihood_surface$lik==max(likelihood_surface$lik))]

#if you want to rerun the simulation uncomment
# run_fits(rep_plot, nn, cut_off, dt, q, q_local, filename)

# Plot the results
plot_outputs(rep_plot, cut_off, filename="1")

