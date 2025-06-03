import numpy as np
import cornerplot
from setup import *
import pdb
import os
import matplotlib.pyplot as plt


def plot_posteriors(posterior_samples):
    # plot cornerplot of posteriors #


    # Load aesthetic info
    parameter_labels = parameters

    # Use cornerplot from https://joss.theoj.org/papers/10.21105/joss.00024
    fig = cornerplot.corner(posterior_samples, labels=parameter_labels,
                             max_n_ticks=3, label_kwargs={"fontsize": 34},
                             show_titles=True, title_kwargs={"fontsize": 34})

    return fig


def plot_datafit(sample_fits, sigma=1):
    # plot the median and n-sigma fits to the data (default is 1-sigma) #

    # find median and 1- or 3-sigma limits for the model at each wavelength point
    if sigma == 3:
        td_min, td_mid, td_max = np.percentile(sample_fits[:, :-1], [0.15, 50, 99.85], axis=0)
    else:
        td_min, td_mid, td_max = np.percentile(sample_fits[:, :-1], [16, 50, 84], axis=0)

    # plot median and n-sigma fits to the data
    fig, ax = plt.subplots(figsize=(15, 10))

    ax.fill_between(wavelength, td_min, td_max, color='#c1adea', label='1 sigma', zorder=2)
    ax.plot(wavelength, td_mid, color='#7447d0', label='Median', zorder=3)
    ax.errorbar(wavelength, spectrum, yerr=error, fmt='o', label='Data',
                zorder=4, ms=8, mew=1, mec='k', mfc='#6495ed', ecolor='#00008b', elinewidth=1, capsize=7)

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.legend(fontsize=20)

    ax.set_xlabel(r'Wavelength ($\mu$m)', fontsize=20)
    ax.set_ylabel(r'Transit Depth (ppm)', fontsize=20)


    return fig


def plot_bestfit(bestfit_model):
    # plot best-fit (maximum likelihood) model and data #

    fig, ax = plt.subplots(figsize=(15, 10))

    ax.plot(wavelength, bestfit_model, color='#7447d0', label='Best fit', zorder=3)
    ax.errorbar(wavelength, spectrum, yerr=error, fmt='o', label='Data',
                zorder=4, ms=8, mew=1, mec='k', mfc='#6495ed', ecolor='#00008b', elinewidth=1, capsize=7)

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.legend(fontsize=20)

    ax.set_xlabel(r'Wavelength ($\mu$m)', fontsize=20)
    ax.set_ylabel(r'Transit Depth (ppm)', fontsize=20)

    return fig
