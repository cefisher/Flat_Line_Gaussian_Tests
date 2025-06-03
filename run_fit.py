import pymultinest
import numpy as np
from setup import *
import ns_setup
import model
from matplotlib import pyplot as plt
import cornerplot
import time
import json
import pdb
import evaluate_samples
import plots
import os

start = time.time()

n_params = len(parameters)

# Set up output path #
retrieval_path = 'example/retrievals/'+model_name+'/'

output_directory = retrieval_path+'out/'
os.makedirs(os.path.dirname(output_directory), exist_ok=True)


# Set up priors #
b = ns_setup.Priors(1, n_params, prior_dict, parameter_dict, parameters, wavelength, spectrum, error, model_name,
                        spectrum_scatter, wave_mid_mol=wmm)

print('*** Starting Retrieval ***')

# Run PyMultinest #
pymultinest.run(b.loglike, b.prior, n_params,
                    outputfiles_basename=output_directory,
                    resume=False, verbose=True,
                    n_live_points=live)


json.dump(parameters, open(output_directory+'params.json', 'w'))  # save parameter names

a = pymultinest.Analyzer(outputfiles_basename=output_directory, n_params=n_params)
samples_likelihood = a.get_equal_weighted_posterior()
samples = samples_likelihood[:, :-1]
bestfit_params = a.get_best_fit()

# Compute result median and 1-sigma limits #
retrieved_results = list(map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0))))

# Print results #
print("""PyMultinest result:""")
for i in range(len(parameters)):
    print(parameters[i],' = ', retrieved_results[i][0],' +',retrieved_results[i][1],' -',retrieved_results[i][2])

# Find Bayesian evidence #
a_lnZ = a.get_stats()['global evidence']
log_evidence = a_lnZ

print ('************************')
print ('MAIN RESULT: Evidence Z ')
print ('************************')
print ('  log Z for model = %.1f' % (log_evidence))

# Save evidence to output folder #
np.save(output_directory+'logevidence', np.array([log_evidence]))

# Make some nice plots #
os.makedirs(os.path.dirname(retrieval_path+'plots/'), exist_ok=True)

# corner plot for the posteriors #
fig1 = plots.plot_posteriors(samples)
fig1.savefig(retrieval_path+'plots/cornerplot.pdf')

# evaluate model for all samples in the posterior #
sample_fits = evaluate_samples.eval_all(samples_likelihood)
np.save(output_directory+ 'sample_fits', sample_fits)

# find model with maximum likelihood #
bestfit_index = np.argmax(sample_fits[:,-1])
bestfit_model = sample_fits[bestfit_index, :-1]

# compute the chi2 of the fit #
if model_name == 'slope_zero_inflatederr':
    error_inflation = samples[:,-1]
    error_bf = error_inflation[bestfit_index]*error
    chi_sq = np.sum(((bestfit_model - spectrum) ** 2) / (error_bf ** 2))
else:
    chi_sq = np.sum(((bestfit_model - spectrum)**2)/(error**2))

# compute and save reduced chi2 #
dof = len(wavelength) - len(parameters)
chi_sq_n = chi_sq/dof
print(r'  Chi^2/N = %.2f' % chi_sq_n)
np.save(output_directory+'chisq', np.array([chi_sq_n]))

# compute and save reduced chi2 values for all samples in the posterior #
chi_sq_array = np.zeros(len(sample_fits))
for i in range(len(sample_fits)):

    if model_name == 'slope_zero_inflatederr':
        error_val = error_inflation[i]*error
    else:
        error_val = error

    chi2 = np.sum(((sample_fits[i,:-1] - spectrum)**2)/error_val**2)
    chi_sq_array[i] = chi2/dof

np.save(output_directory+'chisq_all', chi_sq_array)

# plot the retrieval fit to the data #
fig2 = plots.plot_datafit(sample_fits, sigma=1)
fig2.savefig(retrieval_path+'plots/retrievalfit.pdf')

# plot the best-fit model #
fig3 = plots.plot_bestfit(bestfit_model)
fig3.savefig(retrieval_path+'plots/bestfit.pdf')

# print the time for the retrieval #
end = time.time()
print("time = ",end-start)