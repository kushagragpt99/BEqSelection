##################################################################
################### READ ME for code for #########################
################### Gupta, Vats, Chatterjee:   ###################
############## Bayesian Equation Selection paper #################
##################################################################

The final code used in the paper is in final_implementation, and the additional code can be found in additional_code.
There is a running script "final_implementation/plots/minimum.R" that produces all the 
Figures and the contents of the Tables. However, 
the script "minimum.R" loads objects produced by other codes
that take a while to run. Thus, the output files are provided as well in final_implementation/plots. 
These output files are copies of the output files inside each example folder in final_implementation.

####
"lorenz63/l63_jitter_1e5_sd_*" are produced by "lorenz63/l63_jitter.R". This may take about 11 hours to run.

"lorenz63/truth/l63_linch_T_20_5e4_cwise_1_spikes_init_theta_try" are produced by "lorenz63/truth/l63_spike_slab.R". This may take about 19 hours to run.
"lorenz63/truth/l63_linch_2e6_fin_var" are produced by "lorenz63/truth/lorenz_63_linchpin.R". This may take about 40 hours to run.

"lorenz63/interp/l63_linch_T_20_*e5_cwise_1_spikes_interp_diffuse_6_by_10_scale_try" are produced by "lorenz63/interp/l63_spike_interp.R". This may take about 29 hours to run.
"lorenz63/interp/l63_interp_5e6" are produced by "lorenz63/interp/l63_interp.R". This may take about 40 hours to run.

####
"lorenz96/truth/l96_1e5_MH" are produced by "lorenz96/truth/l96_RWMH.R". This may take about 1.4 hours to run.

"lorenz96/truth/l96_5e4_cwise_spikes_truth_diffuse_try" are produced by "lorenz96/truth/l96_spike_slab.R". This may take about 20 hours to run.
"lorenz96/truth/l96_linch_2e6_fin_var" are produced by "lorenz96/truth/lorenz_96_fin.R". This may take about 38 hours to run.

"lorenz96/interp/l96_*e5_cwise_spikes_interp_diffuse_init_theta_try" are produced by "lorenz96/interp/l96_spike_interp.R". This may take about 19 hours to run.
"lorenz96/interp/l96_linch_5e6_interp" are produced by "lorenz96/interp/l96_interp.R". This may take about 38 hours to run.

####
"OU/truth/ou_linch_spike_truth_T2_1e5_try" are produced by "OU/truth/spike_slab.R". This may take about 1 hour to run.
"OU/truth/OU_linch_1e6_T2_fin" are produced by "OU/truth/ou_linchpin.R". This may take about 8 hours to run.

"OU/interp/ou_linch_spike_interp_T2_1e5_try" are produced by "OU/interp/spike_interp.R". This may take about 1 hour to run.
"OU/interp/OU_linch_interp_1e6_T2_fin" are produced by "OU/interp/ou_interp_linch.R". This may take about 8 hours to run.

