# Shallow water simple model

#Developing the simple model based on the SWE from Würsch and Craig 2013 (WC13)

Contains files necessary to run SWM with a large ensemble. Allows for saving across multiple netCDF files.
Initiatialise, then begin with data assimilation (with the 'DA_init_vis_dt__' function) and then free run(with the 'model_initss__' function). Jupyter notebook provides (small) example of how to run code. 'Run-script.py' provides a script which can be used to run the model on SLURM.

Details of this version: 
- base-version
- 1D 
- saves to netcdf files 
- use parameters which are already set
- 'ss__' attached to function-ends indicates that data is saved every tf timesteps.
- There will be errors if mindex is set to anything but 0. 


