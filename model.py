import numpy as np
import pdb


class Model:

    def __init__(self, x_array, parameter_dict, spectrum_scatter, wave_mid_mol=None):
        '''
        Setting up the priors for the parameters and computing the likelihood 

        Args:
            x_array (np.array):
                Wavelength array for the data.
            parameter_dict (dictionary):
                Dictionary of all parameters and their default values.
            spectrum_scatter (float):
                Standard deviation of the spectrum.        
            wave_mid_mol (float, optional):
                Mid point of the feature for a specific molecule.

        '''

        self.x_array = x_array
        self.parameter_dict = parameter_dict
        self.spectrum_scatter = spectrum_scatter
        self.wave_mid_mol = wave_mid_mol


    def slope_zero(self):
        ## flat line ##

        # compute flat line
        y = self.parameter_dict['y_intercept']*np.ones(len(self.x_array))

        return y


    def slope_line(self):
        ## sloped line ##

        # compute sloped line
        intercept = self.parameter_dict['y_intercept']
        gradient = self.parameter_dict['gradient']

        y = gradient*self.x_array + intercept

        return y


    def slope_zero_offset(self):
        ## flat line with an offset between NRS1 and NRS2

        # split between NRS1 and NRS2
        split = np.where(self.x_array > 3.8)[0][0]
        x_nrs1 = self.x_array[:split]
        x_nrs2 = self.x_array[split:]

        # flat lines covering the two wavelength regions, with an offset between them
        y_nrs1 = self.parameter_dict['y_intercept']*np.ones(len(x_nrs1))
        y_nrs2 = (self.parameter_dict['y_intercept']+self.parameter_dict['nrs2_offset'])*np.ones(len(x_nrs2))

        y = np.append(y_nrs1, y_nrs2)
        
        return y



    def slope_zero_gaussian(self):
        ## flat line with a gaussian feature ##

        # compute flat line
        c = self.parameter_dict['y_intercept']*np.ones(len(self.x_array))

        # compute gaussian feature
        log_amp = self.parameter_dict['log_amplitude']
        wave_mid = self.parameter_dict['wave_mid']
        log_sigma = self.parameter_dict['log_feature_width']
        
        amp = (10**log_amp)*self.spectrum_scatter
        sigma = 10**log_sigma

        gaussian_feature = amp*np.exp(-((self.x_array - wave_mid)**2)/(2*sigma**2))

        # add gaussian feature to flat line
        y = c + gaussian_feature

        return y


    def slope_zero_neg_gaussian(self):
        ## flat line with a negative gaussian feature ##

        # compute flat line
        c = self.parameter_dict['y_intercept'] * np.ones(len(self.x_array))

        # compute gaussian feature
        log_amp = self.parameter_dict['log_amplitude']
        wave_mid = self.parameter_dict['wave_mid']
        log_sigma = self.parameter_dict['log_feature_width']

        amp = (10 ** log_amp) * self.spectrum_scatter
        sigma = 10 ** log_sigma

        gaussian_feature = amp * np.exp(-((self.x_array - wave_mid) ** 2) / (2 * sigma ** 2))

        # subtract gaussian feature from flat line
        y = c - gaussian_feature

        return y


    def slope_zero_offset_gaussian(self):
        ## flat line with an offset between NRS1 and NRS2 and a gaussian feature ##

        # split between NRS1 and NRS2
        split = np.where(self.x_array > 3.8)[0][0]
        x_nrs1 = self.x_array[:split]
        x_nrs2 = self.x_array[split:]

        # set flat lines offset between NRS1 and NRS2
        c_nrs1 = self.parameter_dict['y_intercept'] * np.ones(len(x_nrs1))
        c_nrs2 = (self.parameter_dict['y_intercept'] + self.parameter_dict['nrs2_offset']) * np.ones(len(x_nrs2))
        c = np.append(c_nrs1, c_nrs2)

        # compute gaussian feature
        log_amp = self.parameter_dict['log_amplitude']
        wave_mid = self.parameter_dict['wave_mid']
        log_sigma = self.parameter_dict['log_feature_width']

        amp = (10 ** log_amp) * self.spectrum_scatter
        sigma = 10 ** log_sigma

        gaussian_feature = amp * np.exp(-((self.x_array - wave_mid) ** 2) / (2 * sigma ** 2))

        # add gaussian feature to flat lines
        y = c + gaussian_feature

        return y

    def slope_line_gaussian(self):
        ## sloped line with a gaussian feature ##

        # compute sloped line
        intercept = self.parameter_dict['y_intercept']
        gradient = self.parameter_dict['gradient']

        c = gradient * self.x_array + intercept

        # compute gaussian feature
        log_amp = self.parameter_dict['log_amplitude']
        wave_mid = self.parameter_dict['wave_mid']
        log_sigma = self.parameter_dict['log_feature_width']
        
        amp = (10**log_amp)*self.spectrum_scatter
        sigma = 10**log_sigma

        gaussian_feature = amp * np.exp(-((self.x_array - wave_mid) ** 2) / (2*sigma ** 2))

        # add gaussian feature to slope
        y = c + gaussian_feature

        return y


    def slope_zero_gaussianMOL(self):
        ## flat line with a gaussian feature at fixed wavelength ##

        # compute flat line
        c = self.parameter_dict['y_intercept'] * np.ones(len(self.x_array))

        # compute gaussian feature with fixed mid point
        log_amp = self.parameter_dict['log_amplitude']
        wave_mid = self.wave_mid_mol
        log_sigma = self.parameter_dict['log_feature_width']

        amp = (10 ** log_amp) * self.spectrum_scatter
        sigma = 10 ** log_sigma

        gaussian_feature = amp * np.exp(-((self.x_array - wave_mid) ** 2) / (2 * sigma ** 2))

        # add gaussian feature to flat line
        y = c + gaussian_feature

        return y

    def slope_line_gaussianMOL(self):
        ## sloped line with a gaussian feature at fixed wavelength ##

        # compute sloped line
        intercept = self.parameter_dict['y_intercept']
        gradient = self.parameter_dict['gradient']

        c = gradient * self.x_array + intercept

        # compute gaussian feature with fixed mid point
        log_amp = self.parameter_dict['log_amplitude']
        wave_mid = self.wave_mid_mol
        log_sigma = self.parameter_dict['log_feature_width']

        amp = (10 ** log_amp) * self.spectrum_scatter
        sigma = 10 ** log_sigma

        gaussian_feature = amp * np.exp(-((self.x_array - wave_mid) ** 2) / (2 * sigma ** 2))

        # add gaussian feature to slope
        y = c + gaussian_feature

        return y


    def slope_zero_offset_gaussianMOL(self):
        ## flat line with an offset between NRS1 and NRS2 and a gaussian feature at fixed wavelength ##

        # split NRS1 and NRS2
        split = np.where(self.x_array > 3.8)[0][0]
        x_nrs1 = self.x_array[:split]
        x_nrs2 = self.x_array[split:]

        # flat lines covering the two wavelength regions, with an offset between them
        c_nrs1 = self.parameter_dict['y_intercept'] * np.ones(len(x_nrs1))
        c_nrs2 = (self.parameter_dict['y_intercept'] + self.parameter_dict['nrs2_offset']) * np.ones(len(x_nrs2))
        c = np.append(c_nrs1, c_nrs2)

        # compute gaussian feature with fixed mid point
        log_amp = self.parameter_dict['log_amplitude']
        wave_mid = self.wave_mid_mol
        log_sigma = self.parameter_dict['log_feature_width']

        amp = (10 ** log_amp) * self.spectrum_scatter
        sigma = 10 ** log_sigma

        gaussian_feature = amp * np.exp(-((self.x_array - wave_mid) ** 2) / (2 * sigma ** 2))

        # add gaussian feature to flat lines
        y = c + gaussian_feature

        return y
