import os, re
import glob
import numpy as np

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from functools import wraps
import time

from GEN_Utils import FileHandling
from loguru import logger

logger.info('Hi new user')

def timed(func):
    """This decorator prints the execution time for the decorated function."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        logger.debug("{} ran in {}s".format(func.__name__, round(end - start, 2)))
        return result
    return wrapper

@timed
def pixel_cleaner(results_path, image_name):
    # read into pandas df
    raw_pixels = pd.read_csv(results_path)
    raw_pixels

    # Clean pixels info to generate single list with ROI, coords and intensity
    cleaning_pixels = raw_pixels.T.reset_index()
    cleaning_pixels

    cleaning_pixels['ROI_name'] = cleaning_pixels['index'].str.split('_').str[0]
    grouped_ROI = cleaning_pixels.set_index('ROI_name', drop=True).groupby('ROI_name')

    ROI_data_list = []
    for ROI, group in grouped_ROI:
        new_col_list = list(group['index'].str.split('_').str[1:3])
        ROI_data = group.T
        # Rename columns and drop previous column names
        ROI_data.columns = ['_'.join(x) for x in new_col_list]
        #ROI_data.columns = ['pos_X', 'pos_Y', 'Intensity_C1', 'Intensity_C2', 'Intensity_C3', 'Intensity_C4']
        ROI_data.drop('index', axis=0, inplace=True)
        # Remove leftover pixels from where imageJ assigns 0, 0 to the row - if all values are zero, we assume this is what happened
        ROI_data = ROI_data.replace(0, np.nan)
        ROI_data = ROI_data.dropna(how='all')
        ROI_data = ROI_data.replace(np.nan, 0)
        # add description columns for ROI name and image name
        ROI_data['ROI_name'] = ROI
        ROI_data['image_name'] = image_name
        ROI_data_list.append(ROI_data)

    ROI_pixels = pd.concat(ROI_data_list)
    ROI_pixels.reset_index(inplace=True, drop=True)
    ROI_pixels['X,Y'] = list(zip(ROI_pixels.X_pos, ROI_pixels.Y_pos))

    return ROI_pixels

@timed
def image_processor(input_folder, image_name, output_path, save_file=False):
    logger.info(f'Processing ROI pixels from {image_name}')
    results_path = input_folder+image_name+'_cells_pixels.csv'
    ROI_pixels = pixel_cleaner(results_path, image_name)
    ROI_pixels.rename(columns={'cells' : 'Intensity'}, inplace=True)
    logger.debug(ROI_pixels.columns.tolist())

    # Repeat for nuclear pixel info
    nuclear_path = input_folder+image_name+'_nuclei_pixels.csv'
    logger.info(f'Processing nuclear pixels from {image_name}')
    nuclei_pixels = pixel_cleaner(nuclear_path, image_name)
    nuclei_pixels.rename(columns={'nuclei' : 'Intensity'}, inplace=True)
    logger.debug(nuclei_pixels.columns.tolist())


    # Filter ROI_pixels to remove any pixels that are found in the nuclear pixels list
    logger.info(f'Filtering nuclear pixels from ROI pixels for {image_name}')
    cyto_pixels = ROI_pixels[~ROI_pixels['X,Y'].isin(nuclei_pixels['X,Y'])]
    logger.debug(cyto_pixels.columns.tolist())


    # Collect nuclear pixels that were excluded from whole cell pixels as nuclear pixels for each cell
    nuclear_pixels = ROI_pixels[~ROI_pixels['X,Y'].isin(cyto_pixels['X,Y'])]

    if save_file:
        # Save to excel file:
        logger.info(f'Saving individual excel files for {image_name}')
        FileHandling.df_to_excel(output_path+f'{image_name}_pixel_filtered.xlsx', ['cytoplasm', 'nucleus', 'whole_cell'], data_frames=[cyto_pixels, nuclear_pixels, ROI_pixels])

    logger.info(f'Pixels processed successfully for {image_name}! Results saved to {output_path}')

    return {'cyto':cyto_pixels, 'nuclei_pixels':nuclei_pixels, 'nucleus':nuclear_pixels, 'whole_cell':ROI_pixels}

def sorted_nicely( l ):
    """ Sorts the given iterable in the way that is expected.

    Required arguments:
    l -- The iterable to be sorted.

    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key = alphanum_key)


@timed
def roi_matcher(df_a, df_b, col='coords', x_tolerance=20, y_tolerance=10):
    roi_dict = {}
    for roi in df_a.index.tolist():
        x_pos, y_pos = df_a.loc[roi, col]
        matching_coords = [(x_new, y_new) for (x_new, y_new) in df_b[col] if abs(x_new - x_pos) < x_tolerance if abs(y_new - y_pos) < y_tolerance ]
        if len(matching_coords) == 1:
            roi_dict[roi] = df_b[df_b[col] == matching_coords[0]].index.tolist()[0]
            logger.info(f'Matching coordinates found for {roi}')
        elif len(matching_coords) > 1:
            logger.info(f'Multiple matching coordinates found. Unable to match ambiguous ROI at {roi}')
        else:
            logger.info(f'No matching coordinates found for {roi}.')
    return roi_dict


@timed
def scattbar_plotter(summary_df, group_xcol, y_col, hue_col, output_path, order=None, hue_order=None):

    fig, ax = plt.subplots()

    br = sns.barplot(x=group_xcol, y=y_col, data=summary_df, hue=hue_col, dodge=True, errwidth=1.25, alpha=0.25, ci=None, order=order, hue_order=hue_order, ax=ax)
    scat = sns.swarmplot(x=group_xcol, y=y_col, data=summary_df, hue=hue_col, dodge=True, order=order, hue_order=hue_order, ax=ax)

    leg_label = hue_col
    # To generate custom error bars
    if not order:
        order = list(summary_df[group_xcol].unique())
        logger.debug(f'Samples: {order}')

    if hue_order:
        hue_map = dict(zip(hue_order,range(len(hue_order))))
        summary_df['hue_pos'] = summary_df[hue_col].map(hue_map)
        hue_col = 'hue_pos'
    number_groups = len(list(set(summary_df[hue_col])))
    logger.debug(f'Number of groups: {number_groups}')

    bars = br.patches
    xvals = [(bar.get_x() + bar.get_width()/2) for bar in bars]
    xvals.sort()
    yvals = summary_df.groupby([group_xcol, hue_col]).mean().T[order].T[y_col]

    yerr = summary_df.groupby([group_xcol, hue_col]).std().T[order].T[y_col]

    (_, caps, _) = ax.errorbar(x=xvals,y=yvals,yerr=yerr,capsize=4,elinewidth=1.25,ecolor="black", linewidth=0)
    for cap in caps:
        cap.set_markeredgewidth(2)
    ax.set_ylabel(y_col)

    # To only label once in legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[0:number_groups], labels[0:number_groups], bbox_to_anchor=(1.26, 1.05), title=leg_label)

    # rotate tick labels

    for label in ax.get_xticklabels():
        label.set_rotation(45)

    ax.set_xlabel(group_xcol)

    plt.tight_layout()
    plt.autoscale()
    #plt.show()

    FileHandling.fig_to_pdf([fig], output_path+f'_scattbar')
    FileHandling.fig_to_svg([f'scattbar'], [fig], output_path)

    return fig


def track_plotter(df, track_col):
    fig, ax = plt.subplots(figsize=(10.24, 10.24))
    for track, track_df in df.groupby(track_col):
        sns.lineplot(track_df['X_pos'], track_df['Y_pos'], linewidth=3)
        track_num = track.split('_')[-1]
        plt.annotate(track_num, (round(track_df.mean()['X_pos'], 0)+10, round(track_df.mean()['Y_pos'], 0)+10))
    plt.xlim(0, 1024)
    plt.ylim(0, 1024)

    return fig
