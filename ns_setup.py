import numpy as np
import model
from priors import Priors2
import pdb

pri=Priors2()

class Priors:
    '''
    Priors class for parameter priors and the likelihood calculation 
    '''

    def __init__(self, ndim, nparams, prior_dict, parameter_dict, parameters, x_array, ydata, yerr, model_name, spectrum_scatter, wave_mid_mol=None):
        '''
            Setting up the priors for the parameters and computing the likelihood 

            Args:
                ndim (float):
                    Number of dimensions (usually 1).
                nparams (float):
                    Number of parameters.
                prior_dict (dictionary):
                    Dictionary with the minimum and maximum of the prior for all the parameters. 
                parameter_dict (dictionary):
                    Dictionary of all parameters and their default values.
                parameters (list):
                    List of parameters used in current model.
                x_array (np.array):
                    Wavelength array for the data.
                ydata (np.array):
                    Transit depth array for the data.
                yerr (np.array):
                    Transit depth error array for the data.
                model_name (string):
                    Name of the model being run.
                spectrum_scatter (float):
                    Standard deviation of the spectrum.
                wave_mid_mol (float, optional):
                    Mid point of the feature for a specific molecule.

            '''

        self.ndim = ndim
        self.nparams = nparams
        self.prior_dict = prior_dict
        self.parameter_dict = parameter_dict
        self.parameters = parameters
        self.x_array = x_array
        self.ydata = ydata
        self.yerr = yerr
        self.model_name = model_name
        self.spectrum_scatter = spectrum_scatter
        self.wave_mid_mol = wave_mid_mol


    def prior(self, cube, ndim, n_params):
        # setting up the uniform priors #

        for i,param in enumerate(self.parameters):
                cube[i] = pri.UniformPrior(cube[i], self.prior_dict[param][0], self.prior_dict[param][1])

    def loglike(self, cube, ndim, n_params):
        ## running the model and computing the likelihood ##

        # set the parameter value #
        for i,param in enumerate(self.parameters):
            self.parameter_dict[param] = cube[i]

        # set up the model object #
        y_obj = model.Model(self.x_array, self.parameter_dict, self.spectrum_scatter, self.wave_mid_mol)

        # compute the model depending on the model_name used #
        if self.model_name == 'slope_zero':
            ymodel = y_obj.slope_zero()
        elif self.model_name == 'slope_zero_inflatederr':
            ymodel = y_obj.slope_zero()
        elif self.model_name == 'slope_line':
            ymodel = y_obj.slope_line()
        elif self.model_name == 'slope_zero_offset':
            ymodel = y_obj.slope_zero_offset()
        elif self.model_name == 'slope_zero_gaussian':
            ymodel = y_obj.slope_zero_gaussian()
        elif self.model_name == 'slope_zero_neg_gaussian':
            ymodel = y_obj.slope_zero_neg_gaussian()
        elif self.model_name == 'slope_line_gaussian':
            ymodel = y_obj.slope_line_gaussian()
        elif self.model_name == 'slope_zero_offset_gaussian':
            ymodel = y_obj.slope_zero_offset_gaussian()
        elif self.model_name in ['slope_zero_gaussianCH4', 'slope_zero_gaussianCO2']:
            ymodel = y_obj.slope_zero_gaussianMOL()
        elif self.model_name in ['slope_line_gaussianCH4', 'slope_line_gaussianCO2']:
            ymodel = y_obj.slope_line_gaussianMOL()
        elif self.model_name in ['slope_zero_offset_gaussianCH4', 'slope_zero_offset_gaussianCO2']:
            ymodel = y_obj.slope_zero_offset_gaussianMOL()
        else:
            print('Unknown model selected')

        # inflate the error if the inflatederr model is used #
        if self.model_name == 'slope_zero_inflatederr':
            yerr_inflated = self.parameter_dict['error_inflation']*self.yerr

            # compute the chi2 #
            chi2 = ((ymodel - self.ydata)**2/(yerr_inflated**2)).sum()/(len(self.ydata)-len(self.parameters))
            # evaluate likelihood #
            loglikelihood = (-0.5 * ((ymodel - self.ydata) / yerr_inflated) ** 2 - np.log(
                    abs(self.yerr) * np.sqrt(2 * np.pi))).sum()
        else:
            # evaluate loglikelihood #
            loglikelihood = (-0.5 * ((ymodel - self.ydata) / self.yerr) ** 2 - np.log(
                abs(self.yerr) * np.sqrt(2 * np.pi))).sum()

        return loglikelihood
