# frap_par_estimation
Sample Matlab code for estimating parameters from FRAP microscopy data

This cose accompanies the article "Ciocanel, Maria-Veronica, et al. "Analysis of active transport by fluorescence recovery after photobleaching." Biophysical journal 112.8 (2017): 1714-1725." This sample code uses an average dataset of fluorescence recovery and assumes models protein dynamics in the bleach spot using a 2-state diffusion-transport partial differential equation model. 

The code parameter_sweep_2state.m carries out an initial sweep through relevant parameters, which can be run in parallel. This code calls find_bt.m which sets up an initial value for a parameter defining the initial bleach spot model. It also calls fcn_2state.m which simulates a sample FRAP curve given known parameter.

Once some candidate parameter sets are chosen from the sweep, the code find_par-estimates_2state.m performs parameter estimation using the Matlab multisearch function, which can be run in parallel. Besides the functions above, this code also calls function parameter_2state.m which evaluates the residual between simulated FRAP curves and data.
