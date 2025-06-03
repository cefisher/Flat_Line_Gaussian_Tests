import numpy as np
import sys
import matplotlib.pyplot as plt
import pdb
from scipy.special import erfc

path = 'example/'

models = ['slope_zero', 'slope_line', 'slope_zero_gaussian', 'slope_line_gaussian',
          'slope_zero_gaussianCH4', 'slope_line_gaussianCH4',
          'slope_zero_gaussianCO2', 'slope_line_gaussianCO2',
          'slope_zero_neg_gaussian']


cell_text = []

model_dict = {}
logz_list = []

# load evidence and chi2 for each model
for m in models:
    logz = np.load(path+'retrievals/'+m+'/out/logevidence.npy')
    chisq = np.load(path+'retrievals/'+m+'/out/chisq.npy')

    model_dict[m] = [logz[0], chisq[0]]
    logz_list.append(logz[0])


# choose maximum evidence from null hyptheses
maximum_null = max(model_dict['slope_zero'][0], model_dict['slope_line'][0])
max_logz = max(logz_list)


# compute bayes factor for each model and add info to cell text
for m in models:

    bayes_factor = max_logz - model_dict[m][0]
    model_dict[m].append(bayes_factor)

    cell_text.append(['%.2f' % model_dict[m][0], '%.2f' % model_dict[m][1], '%.2f' % model_dict[m][2]])


# table column headings
columns = ['logZ', 'Chi^2/N', 'log Bayes Factor']

# make table of data
fig, ax = plt.subplots(figsize=(10,5))

# hide axes
fig.patch.set_visible(False)
ax.axis('off')
ax.axis('tight')

ax.table(cellText=cell_text, rowLabels=models, colLabels=columns, loc='center')

ax.set_title(path, fontsize=16)

fig.tight_layout()
fig.savefig(path+'Results_Table.pdf')