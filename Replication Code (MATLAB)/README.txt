README FILE

Helper Functions:
        1.) DSGEautocovar - computes ACF for parameters theta at length h
        2.) DSGEforecasterror - computes Forecast errors up to length h
        3.) DSGEspectrum - computes spectrum of DSGE model at frequencies omega
        4.) ImpulseResponse - obtain DSGE impulse response function from i to j
        5.) kal_wrapper - wrapper for kalman filter that outputs time series
                of likelihood and filtered states
        6.) kalman_filter - filter based on a single date of observables
        7.) particle_wrapper - wrapper for particle filter that outputs time series objects
        8.) particle_filter - performs particle filtering


month2quarter - inputs csv files for output, labor, inflation and interest rates and converts them all into quarterly observations and correct units (logs, growth rates, etc…)
	      -pictures go into data/data_pics
    Plots:
	1.) lab_share_full - plots labor share from 1965:2014
	2.) cons_share_S_ND - plots consumption share of nominal GDP 1965:2014

All Plots go into folder Matlab/figuresv2

Codes:
	1.) empiricalanalogues.m - plots ACF, spectrum and IRF based on real data
	2.) irfmatching.m        - performs IRF matching to estimate model
	3.) liklihood.m          - kalman and particle filtering
	4.) methodofmoments.m    - Moment estimation of model
	5.) modelimplications.m  - model implied ACF, spectrum and IRF
Folders:
	1.) Metropolis
	2.) SMC
