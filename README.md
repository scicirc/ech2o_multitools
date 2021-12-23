## ech2o_multitools

### A multi-purpose python script for to use with the ecohydrological model 
[EcH2O-iso](https://bitbucket.org/scicirc/ech2o-iso): 
- Sensitivity analysis 
- Calibration
- Ensemble runs

--------------

### General use

The script is always launched through the main file ``ech2o_multitools.py``, also specifying a usage-specific 
*definition file* (e.g. def_file.py) where the characteristic for the are provded:

    python ech2o_multitools.py --file=def_file.py (+other options...)
    
Notably, the class Config in the definition file specifies one of 5 *modes* for the script : 
1. *calib_MCsampling*: brute-force Monte Carlo sampling, writing parameters sets in text files
2. *calib_MCruns*: runs the EcH2O-iso model reading the generated parameters, and storing goodness-of-fits for time series
3. *calib_SPOTPY*: a more advanced calibration mode using the [SPOTPY package](https://spotpy.readthedocs.io/en/latest/). 
For now, only the Differential Evolution Adaptative Metropolis (DREAM; Vrugt et al., 2011) has been adapted and is still under testing 
4. *forward_runs*: runs the model for an ensemble of parameter sets (e.g. resulting from a prior calibration) and store time series and stacked time-varying maps (netCDF files).
                 Allows for looking at observations not used in the calibration.
5. *sensi_morris*: sensitivity analysis based on the [Morris method](https://en.wikipedia.org/wiki/Morris_method)

