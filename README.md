# Flat Line Gaussian Tests
Tests flat lines and gaussian feature fits to data to check for evidence of spectral features.

This code requires [PyMultinest](http://johannesbuchner.github.io/PyMultiNest/) to run, so make sure you have this downloaded first. 

## Requirements 
This code is developed using Python 3 and requires the following packages:

- PyMultinest
- numpy
- matplotlib

## Setting up 
In the ```setup.py``` file, you need to give the code the data file, which includes the wavelength, transit depth, and error bars for the planet's spectrum. The file should be in the same [format as in BeAR retrievals](https://newstrangeworlds.github.io/BeAR/sections/observations.html#input-file-structure).
```
# load data. Assumes datafile is in the format for BeAR retrievals. #
datafile = 'example/WASP-15b_G395H.dat'
```

In ```run_fit.py```, you need to give the code the path for the fits to go into. 
```
# Set up output path #
retrieval_path = 'example/retrievals/'+model_name+'/'
```

## Running the fits
To perform a fit for a particular model, for example the ```slope_zero``` model, you can run
```
python3 run_fit.py slope_zero
```

To run all the models in succession, you can run
```
python3 run_all_models.py
```
This will run all the models and make a results table. The path for the results table just needs to be set inside ```make_results_table.py```
```
path = 'example/'
```

