import numpy as np
import model
from setup import *


def eval_all(samples):
    # evaluate the model at every sample point from the posteriors #

    len_x = len(wavelength)

    # set up numpy array for the models
    sample_evaluations = np.zeros((len(samples),len_x+1))

    new_param_dict = parameter_dict

    # assign parameter values
    for i in range(len(samples)):
        for j in range(len(samples[0,:-1])):
            new_param_dict[parameters[j]] = samples[i,j]

        # set up model object
        y_obj = model.Model(wavelength, new_param_dict, spectrum_scatter, wmm)

        # compute model based on model_name
        if model_name == 'slope_zero':
            ymodel = y_obj.slope_zero()
        elif model_name == 'slope_zero_inflatederr':
            ymodel = y_obj.slope_zero()
        elif model_name == 'slope_line':
            ymodel = y_obj.slope_line()
        elif model_name == 'slope_zero_offset':
            ymodel = y_obj.slope_zero_offset()
        elif model_name == 'slope_zero_gaussian':
            ymodel = y_obj.slope_zero_gaussian()
        elif model_name == 'slope_zero_neg_gaussian':
            ymodel = y_obj.slope_zero_neg_gaussian()
        elif model_name == 'slope_zero_neg_gaussian47':
            ymodel = y_obj.slope_zero_neg_gaussian47()
        elif model_name == 'slope_line_gaussian':
            ymodel = y_obj.slope_line_gaussian()
        elif model_name == 'slope_zero_offset_gaussian':
            ymodel = y_obj.slope_zero_offset_gaussian()
        elif model_name in ['slope_zero_gaussianCH4', 'slope_zero_gaussianCO2']:
            ymodel = y_obj.slope_zero_gaussianMOL()
        elif model_name in ['slope_line_gaussianCH4', 'slope_line_gaussianCO2']:
            ymodel = y_obj.slope_line_gaussianMOL()
        elif model_name in ['slope_zero_offset_gaussianCH4', 'slope_zero_offset_gaussianCO2']:
            ymodel = y_obj.slope_zero_offset_gaussianMOL()
        else:
            print('Unknown model selected')

        # set sample to model with likelihood
        sample_evaluations[i] = np.append(ymodel, samples[i,-1])  # add likelihood


    return sample_evaluations


