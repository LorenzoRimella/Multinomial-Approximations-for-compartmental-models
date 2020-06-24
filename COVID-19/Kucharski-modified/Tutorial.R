# - - - - - - - 
# set a working directory with the tutorial:
#  **/Tutorial_scripts/COVID-19/Kucharski-modified
# change ** with your directory
setwd("~/Tutorial_scripts/COVID-19/Kucharski-modified")
dropbox_path = ""

# This file loads the dataset, the parameters and all the fuctions that are needed for the analysis
source("R/outputs_main.R")

# - - - - - - - 
# Run 100 particle filters
rep_plot <- 100; nn <- 3e3; cut_off=0; dt=1; filename="1"

# uncomment to Run bootstrap SMC 
# run_fits(rep_plot, nn=3e3, dt=1, filename="1")

# Output plots
plot_outputs(rep_plot, filename="1") # Figure 2

# plot_dispersion(filename="1") # Figure 3


# Run models --------------------------------------------------------------




