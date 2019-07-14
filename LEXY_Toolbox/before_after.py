import os
import numpy as np

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from GEN_Utils import FileHandling#, scatbar_plot

from loguru import logger

logger.info("Import OK")

input_folder = ''
output_path = ''

if not os.path.exists(output_path):
    os.mkdir(output_path)

file_list = [filename for filename in os.listdir(input_folder)]
file_list = [filename for filename in file_list if '.xlsx' in filename]
# read in raw data to dictionary
raw_dict = {}
for filename in file_list:
    result_name = filename.replace('.xlsx', '')
    raw_dict[result_name] = pd.read_excel(f'{input_folder}{filename}', sheet_name=None)
# Check everything was collected correctly. NB: everey file is collected as dictionary accessed via sheetname
raw_dict.keys()
raw_dict['mean_nuclei_summary']['summary_sorted'].keys()

# Assign mean_summary data to new variable for manipulating
mean_nuclei = raw_dict['mean_nuclei_summary']['summary_sorted'].copy()


# Collect only cols of interest
cols_of_interest = ['ROI_name', 'Intensity', 'mutant', 'time', 'nuc_id', 'track_id']
# collect tracks for each mutant that have values in t0 and t7
before = mean_nuclei[mean_nuclei['time'] == 0][cols_of_interest]
after = mean_nuclei[mean_nuclei['time'] == 5][cols_of_interest]
before_after = pd.merge(before, after, how='inner', on=['mutant', 'track_id'], suffixes=('_before', '_after'), copy=True)

# add columns normalised to t0 intensity
before_after['norm_Intensity_before'] = before_after['Intensity_before'] / before_after['Intensity_before']
before_after['norm_Intensity_after'] = before_after['Intensity_after'] / before_after['Intensity_before']

# Collect useful columns, group before/after by mutant and save to excel
mutant_summary = {}
for group, df in before_after.groupby('mutant'):
    mutant_summary[group] = df[['track_id', 'ROI_name_before', 'Intensity_before', 'Intensity_after', 'norm_Intensity_before', 'norm_Intensity_after']]

FileHandling.df_to_excel(f'{output_path}first_last_summary.xlsx', sheetnames=list(mutant_summary.keys()), data_frames=list(mutant_summary.values()))

# melt to produce data for seaborn plotting normalised intensity
plottable = before_after[['mutant', 'track_id', 'norm_Intensity_before', 'norm_Intensity_after']].melt(id_vars=['mutant', 'track_id'], var_name='Activation', value_name='Intensity')

plottable['Activation'] = plottable['Activation'].str.split('_').str[-1]

#fig = scatbar_plot.scatbar_plot(x_col='mutant', y_col='Intensity', plotting_dfs=[plottable], hue_col='Activation', group_col='mutant')
#plt.savefig(f'{output_path}norm_first_last.png')

# melt to produce data for seaborn plotting raw intensity
plottable = before_after[['mutant', 'track_id', 'Intensity_before', 'Intensity_after']].melt(id_vars=['mutant', 'track_id'], var_name='Activation', value_name='Intensity')

plottable['Activation'] = plottable['Activation'].str.split('_').str[-1]

#fig = scatbar_plot.scatbar_plot(x_col='mutant', y_col='Intensity', plotting_dfs=[plottable], hue_col='Activation', group_col='mutant')
#plt.savefig(f'{output_path}first_last.png')
