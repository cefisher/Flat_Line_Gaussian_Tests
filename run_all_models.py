import os
import sys
import subprocess
import pdb

# list of all models to run #
models = ['slope_zero', 'slope_line', 'slope_zero_gaussian', 'slope_line_gaussian',
          'slope_zero_gaussianCH4', 'slope_line_gaussianCH4',
          'slope_zero_gaussianCO2', 'slope_line_gaussianCO2',
          'slope_zero_neg_gaussian']

# iterate through list and run nested-sampling fit #
for m in models:
    os.system('python3 run_fit.py '+ m)
    
os.system('python3 make_results_table.py')