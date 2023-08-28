# Shallow water simple model (SWM)

# Developing the simple model based on the SWE from Würsch and Craig 2013 (WC13)

Contains files necessary to run SWM with a large ensemble. Allows for saving across multiple netCDF files.
Initiatialise, then begin with data assimilation (with the 'DA_init_vis_dt__' function) and then free run(with the 'model_initss__' function). Jupyter notebook (Example_2) provides (small) example of how to run code. 'Run-script.py' provides a script which can be used to run the model on SLURM.

Details of this first version: 
- base-version
- 1D 
- saves to netcdf files 
- use parameters which are already set
- 'ss__' attached to function-ends indicates that data is saved every tf timesteps.
- There will be errors if mindex is set to anything but 0.

Files for this first version: 
- DA_2019.py: Data Assimilation functions
- Example_2.ipynb: Jupyter notebook example
- constants.py: constants needed for functions
- msw_model_expanded.py: functions to run ensemble
- run_script.py: script to run ensemble

This version of the model was used for the simulation ran in the paper titled "Convergence of forecast distributions in a 100,000‐member idealised convective‐scale ensemble" by Tempest, Craig and Brehmer 2023. 

# Extending the model

An extended version of the model was developed to allow for weak and strong forcing convective weather regimes. See the PDF "Convergence_of_forecast_distributions_in_weak_and_strong_forcing_convective_weather_regimes__to_submit_" for details on extensions to model. All files with "_weakregime" and "_strongregime" at end of name are for this extended version, the weak and strong referring to the type of forcing regime. Runs in similar manner as the original version. Files uploaded for weak forcing: 

Analysing scripts:
- basic_analysis_weak_forcing.py: script to create basic parameters to analyse data. Eg convective timescale. 
- basic_analysiss_weakforcing.sh: slurm script to run above script
- Calculate_spread_mean_abserror_weakforcing.py: calculate spread properties of ensemble.
- Calculate_spread_mean_abserrors_weakforcing.sh: slurm script to run script above
- make_bootstrap_neighbourhood_2_weakforcing.py: create bootstraps with neighbourhood distributions
- make_bootstrap_neighbourhood_2s_weakforcing.sh: script to run script above
- make_bootstrap_weakforcing.py: create boostrap from 1 grid points distribution
- make_bootstraps_weakforcing.sh: script to run script above
- make_convergence_eescomparison_plot_weakforcing.py: create convergence measure plots to compare weak and strong forcing
- make_convergence_eescomparison_plots_weakforcing.sh: script to run script above
- make_distribution_data_weakforcing.py
- make_distribution_datas_weakforcing.sh
- make_neighbourhood_distribution_data_weakforcing.py
- make_neighbourhood_distribution_datas_weakforcing.sh

Simulation scripts: 
- constants_weakforcing.py: constants for model
- DA_2019_weakforcing.py: Data Assimilation functions
- msw_model_expanded_4_weakforcing.py: functions to run simulation
- run_script_weakforcing.py: defines which functions to run for ensemble
- run_large_0_weakforcing.sh: slurm script to run ensemble (executes run_script.py)

Strong forcing files are very similar. The following files for the simulation of the strong forcing run are different, however: 

- constants_strongforcing.py: constants for model
- msw_model_expanded_4_strongforcing.py: functions to run simulation
- run_script_strongforcing.py: defines which functions to run for ensemble
- run_large_0_strongforcing.sh: slurm script to run ensemble (executes run_script.py)

