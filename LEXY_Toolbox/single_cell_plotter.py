import os, re, pathlib
import pandas as pd
import numpy as np
from loguru import logger

import matplotlib.pyplot as plt
import seaborn as sns

from GEN_Utils.PlotUtils import FileHandling

logger.info("Import OK")

input_path = ''
output_path = ''

if not os.path.exists(output_path):
    os.mkdir(output_path)

raw_summary = pd.read_excel(input_path)
summary = raw_summary.copy()

normalised_data = {}

for track, grouped_data in summary.groupby('track_id'):
    norm_data = grouped_data[['ROI_name', 'Intensity', 'time']]
    norm_data['norm_Intensity'] = norm_data['Intensity'] / int(norm_data[norm_data['time'] == 0]['Intensity'])
    normalised_data[track] = norm_data



figures = []

for key, item in normalised_data.items():
    fig, ax = plt.subplots()
    sns.lineplot(x='time', y='norm_Intensity', hue=None, size=None, style=None, data=item, palette=None, hue_order=None, hue_norm=None, sizes=None, size_order=None, size_norm=None, dashes=True, markers='o', style_order=None, units=None, estimator="mean", ci=95, n_boot=1000, sort=True, err_style="band", err_kws=None, legend="brief", ax=None)
    ax.set_xlim(0, 5)
    ax.set_ylim(0, 2)
    plt.title(key, fontdict=None, loc='center', pad=None)
    ax.set_xticklabels(list(item['ROI_name']), rotation = 45)
    plt.xlabel('ROI name')

    plt.autoscale()
    plt.tight_layout()

    figures.append(fig)
    fig.savefig(f'{output_path}{key}.png')
|
