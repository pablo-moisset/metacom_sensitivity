# metacom_sensitivity

1- System requirements

PC with Linux, or OS X and MATLAB. The code is simple and should work on pretty much every version of MATLAB but we have only tested it using 2019a on Debian Linux 12 and OS X v10.15. gzip is required

2- Instalation

Download all the provided files in any user directory. 

3- Running the program

The main programs are:
run_sims.sh: The parameters are described in detail in main_cluster_v16.m
This runs simulations and leaves the time series in hlp.gz files

main_star_experiment.m
Runs a star experiment and creates timeseries files

pp_tseries.sh
Reads a  hlp.gz file with simulation time series and computes the response variables in a .mat file

main_regression_foodweb.m
main_regression_landscape.m 
Compute regressions and effects of predictors

main_akaike.m
Compute AIC for different linear models

main_create_SI_panels.m
main_hist_response.m
main_plot_diversity_vs_predictors.m
main_plotpanels.m
main_prelim_fig.m
main_scatter_hist.m
main_star_plots.m
Scripts to create graphics

4- Licence related issues

All program files were written by us and are in the public domain with the exception of csvimport.m, hlp_deserialize.m,  hlp_serialize.m, and randraw.m. These file was obtained from MATLAB central and are used for formatting the plots and generating random numbers. They have no code related to the numerical experiments thenselves.
