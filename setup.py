import numpy as np
import pdb
import sys


def find_bin_width(wavelength):
    # find bin widths if only centre point is given #

    wavelength_bin_edges = 0.5 * (wavelength[:-1] + wavelength[1:])
    wavelength_bin_mins = np.append(2 * wavelength[0] - wavelength_bin_edges[0], wavelength_bin_edges)
    wavelength_bin_maxs = np.append(wavelength_bin_edges, 2 * wavelength[-1] - wavelength_bin_edges[-1])
    wave_bin_width = 0.5 * (wavelength_bin_maxs - wavelength_bin_mins)

    return wave_bin_width


nrs_split = True    # for G395H. False for G395M


# load data. Assumes datafile is in the format for BeAR retrievals. #
datafile = 'example/WASP-15b_G395H.dat'

with open(datafile, 'r') as file:
    df = file.readlines()
    
datatype = df[5].split()[0]
data = np.loadtxt(datafile, skiprows=11)

if datatype == 'spectroscopy':
    wavelength = data[:,0]
    spectrum = data[:,1]
    error = data[:,2]

    if nrs_split == True:
        # first split nrs1 and nrs2
        split = np.where(wavelength > 3.8)[0][0]
        wavelength_nrs1 = wavelength[:split]
        wavelength_nrs2 = wavelength[split:]
        wave_bin_width_nrs1 = find_bin_width(wavelength_nrs1)
        wave_bin_width_nrs2 = find_bin_width(wavelength_nrs2)
        wave_bin_width = np.append(wave_bin_width_nrs1, wave_bin_width_nrs2)
    else:
        wave_bin_width = find_bin_width(wavelength)
elif datatype == 'band-spectroscopy':
    wavelength = 0.5*(data[:,0]+data[:,1])
    spectrum = data[:,2]
    error = data[:,3]
    wave_bin_width = 0.5 * (data[:, 1] - data[:, 0])

    
# standard deviation of the spectrum #
spectrum_scatter = np.std(spectrum)

# find min and max transit depths for prior range of flat line #
transit_depth_min = min(spectrum - error)
transit_depth_max = max(spectrum + error)

# find maximum gradient for prior range of slope line #
max_gradient = (transit_depth_max-transit_depth_min)/(wavelength[-1]-wavelength[0])

# find minimum width of a feature as the biggest bin width #
min_feature_width = np.log10(np.max(wave_bin_width))

# mid points for CH4 and CO2 features #
wave_mid_ch4 = 3.3
wave_mid_co2 = 4.4

## Retrieval Info ##
model_name = sys.argv[1]
live = 1000 # nested-sampling live points

# all priors #
prior_dict = {'y_intercept': [transit_depth_min, transit_depth_max],
              'gradient': [-max_gradient, max_gradient],
              'nrs2_offset': [transit_depth_min-transit_depth_max, transit_depth_max-transit_depth_min],
              'log_amplitude': [-2, 1.0],
              'wave_mid': [wavelength[0]-0.5, wavelength[-1]+0.5],
              'log_feature_width': [min_feature_width, -0.3],
              'error_inflation': [0.1, 1.0]}

# default parameter values #
parameter_dict = {'y_intercept': transit_depth_min,
                  'gradient': 0,
                  'nrs2_offset': 0,
                  'log_amplitude': -2,
                  'wave_mid': wavelength[0],
                  'log_feature_width': min_feature_width,
                  'error_inflation': 1.0}

# parameters for each model #
model_parameters = {'slope_zero': ['y_intercept'],
                    'slope_zero_inflatederr': ['y_intercept', 'error_inflation'],
                    'slope_line': ['y_intercept', 'gradient'],
                    'slope_zero_offset': ['y_intercept', 'nrs2_offset'],
                    'slope_zero_gaussian': ['y_intercept', 'log_amplitude', 'wave_mid', 'log_feature_width'],
                    'slope_zero_neg_gaussian': ['y_intercept', 'log_amplitude', 'wave_mid', 'log_feature_width'],
                    'slope_line_gaussian': ['y_intercept', 'gradient', 'log_amplitude', 'wave_mid', 'log_feature_width'],
                    'slope_zero_offset_gaussian': ['y_intercept', 'nrs2_offset', 'log_amplitude', 'wave_mid', 'log_feature_width'],
                    'slope_zero_gaussianCH4': ['y_intercept', 'log_amplitude', 'log_feature_width'],
                    'slope_line_gaussianCH4': ['y_intercept', 'gradient', 'log_amplitude', 'log_feature_width'],
                    'slope_zero_offset_gaussianCH4': ['y_intercept', 'nrs2_offset', 'log_amplitude', 'log_feature_width'],
                    'slope_zero_gaussianCO2': ['y_intercept', 'log_amplitude', 'log_feature_width'],
                    'slope_line_gaussianCO2': ['y_intercept', 'gradient', 'log_amplitude', 'log_feature_width'],
                    'slope_zero_offset_gaussianCO2': ['y_intercept', 'nrs2_offset', 'log_amplitude', 'log_feature_width']}

parameters = model_parameters[model_name]

if 'CH4' in model_name:
    wmm = wave_mid_ch4
elif 'CO2' in model_name:
    wmm = wave_mid_co2
else:
    wmm = None