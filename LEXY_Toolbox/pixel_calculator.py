import os, re
import glob
import numpy as np
import string

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from GEN_Utils import FileHandling


from analysis_scripts_HTT_edited.pixel_functions import image_processor, sorted_nicely, scattbar_plotter, pixel_cleaner, roi_matcher, track_plotter

from loguru import logger


logger.info("Import OK")

#This functionality has been moved to this 'master' file to make processing simpler
# collect list of files in the input_folder that have .csv extension
input_folder = ''
output_path = ''

trans_threshold = 700 # raw mean intensity in mCherry channel
agg_threshold = 20 # percentage of pixels in a cell needed to be marked as aggregate
GFP_threshold = 1500

# Check output folder exists - if not, create it
if not os.path.exists(output_path):
    os.mkdir(output_path)

file_list = [filename for filename in os.listdir(input_folder) if 'pixels.csv' in filename]
logger.info(f'The following files were detected from {input_folder}:\n {file_list}')

# collect list of files in the input_folder that have .csv extension
image_names = sorted_nicely(list(set([("_".join(name.split("_")[0:5])) for name in file_list if 'cells' not in name])))
logger.info(f'The following image names were detected:\n {image_names}')

mutants = list(set([name.split('_')[0] for name in image_names]))

# Clean ROI pixels for first timepoint to simple dataframe for each image
pixel_dict = {}
summary_dict = {}
for image_name in image_names:
    if 't_0' in image_name:
        image = image_name + '_mCherry'
        cleaned_pixels = image_processor(input_folder, image, output_path,  save_file=False)
        pixel_dict[image_name] = cleaned_pixels

        summary_nuclei = cleaned_pixels['nuclei_pixels'].groupby('ROI_name').mean()
        summary_nuclei['image_name'] = image_name
        summary_nuclei.rename(columns={'mCherry_nuclei' : 'Intensity'}, inplace=True)

        summary_cyto = cleaned_pixels['cyto'].groupby('ROI_name').mean()
        # Calculate % of saturated pixels as a way to mark cells with aggregates
        sat_pixels = []
        for ROI_name, df in cleaned_pixels['cyto'].groupby('ROI_name'):
            sat_pixels.append(df['mCherry_cells'].mask(df['mCherry_cells'] < 65535).count() / df.shape[0] * 100)
        summary_cyto['saturated'] = sat_pixels
        summary_cyto['agg'] = [1 if sat_percent > agg_threshold else 0 for sat_percent in summary_cyto['saturated']]
        summary_cyto['image_name'] = image_name
        # Filter cell ROIs according to transfected
        summary_cyto = summary_cyto[summary_cyto['mCherry_cells'] > trans_threshold]

        summary_dict[image_name] = {'nuclei':summary_nuclei, 'cyto':summary_cyto}
    else:
        nuclear_path = input_folder+image_name+'_mCherry_nuclei_pixels.csv'
        logger.info(f'Processing nuclear pixels from {image_name}')
        nuclei_pixels = pixel_cleaner(nuclear_path, image_name)
        pixel_dict[image_name] = {'nuclei_pixels':nuclei_pixels}

        summary_nuclei = nuclei_pixels.groupby('ROI_name').mean()
        summary_nuclei.rename(columns={'mCherry_nuclei' : 'Intensity'}, inplace=True)
        summary_nuclei['image_name'] = image_name

        summary_dict[image_name] = {'nuclei':summary_nuclei}

logger.info('All images processed successfully.')

# Additional processing for the GFP channel
GFP_summary_dict = {}
t1_images = [image for image in image_names if 't_1' in image]
for image_name in t1_images:
    image = image_name + '_GFP'
    cleaned_pixels = image_processor(input_folder, image, output_path,  save_file=False)
    summary_cyto = cleaned_pixels['cyto'].groupby('ROI_name').mean()
    # Calculate % of saturated pixels as a way to mark cells with aggregates
    sat_pixels = []
    for ROI_name, df in cleaned_pixels['cyto'].groupby('ROI_name'):
        sat_pixels.append(df['GFP_cells'].mask(df['GFP_cells'] < 65535).count() / df.shape[0] * 100)
    summary_cyto['saturated'] = sat_pixels
    summary_cyto['agg'] = [1 if sat_percent > agg_threshold else 0 for sat_percent in summary_cyto['saturated']]
    summary_cyto['image_name'] = image_name
    # Filter cell ROIs according to transfected
    summary_cyto = summary_cyto[summary_cyto['GFP_cells'] > GFP_threshold]

    GFP_summary_dict[image_name] = {'nuclei':summary_nuclei, 'cyto':summary_cyto}

# Determine cotransfected cells, generate summary cyto dataframe
cyto_summary = {}
for mutant in mutants:
    pos_list = list(set([name.split('_')[2] for name in image_names if mutant in name]))
    position_cyto = {}
    for pos in pos_list:
        mCherry_image_name = f'{mutant}_pos_{pos}_t_0'
        mCherry_summary = summary_dict[mCherry_image_name]['cyto']
        GFP_image_name = f'{mutant}_pos_{pos}_t_1'
        GFP_summary = GFP_summary_dict[GFP_image_name]['cyto']

        cotrans_summary = pd.merge(mCherry_summary, GFP_summary.drop(['X_pos', 'Y_pos'], axis=1), how='outer', left_index=True, right_index=True, suffixes=('_mCherry', '_GFP'))

        cotrans_summary['cotrans'] = [1 if roi_name in cotrans_summary.dropna().index else 0 for roi_name in cotrans_summary.index]

        cyto_summary[mCherry_image_name] = cotrans_summary

# Generate compiled summary_df for all images, after adding descriptive nuc_id to each ROI
summary_list = []
for key, summary_df in cyto_summary.items():
    # Add additional info from image names
    summary_df = summary_df.dropna()
    summary_df['mutant'] = summary_df['image_name_mCherry'].str.split('_').str[0]
    summary_df['position'] = summary_df['image_name_mCherry'].str.split('_').str[2]
    summary_df['number'] = [str(number) for number in range(0, summary_df.shape[0])]
    summary_df['cell_id'] = list(zip(summary_df['mutant'], summary_df['position'], summary_df['number']))
    summary_df['cell_id'] = ['_'.join(tuple) for tuple in summary_df['cell_id']]
    summary_dict[key]['cyto'] = summary_df
    summary_list.append(summary_df)
mean_cell_summary = pd.concat(summary_list)

FileHandling.df_to_excel(output_path+f'mean_cell_summary.xlsx', sheetnames=['summary'], data_frames=[mean_cell_summary.reset_index()])


# # Checking out the summary results
# pixel_dict.keys() # shows the image names that were processed i.e. key 1
# pixel_dict['barnaseWT+Q25_pos_001_t_0']['cyto'].keys() # Shows the available dataframes for each image i.e. key 2
# pixel_dict['Ex2_pos_1_t_0']['nuclei_pixels'] # shows the cytoplasm pixel summary for the first image
# summary_dict.keys()
# summary_dict['Ex2_pos_1_t_0'].keys()
#
# summary_dict['barnaseWT+Q25_pos_001_t_0']['cyto']


# plot individual image for t0 to test segmentation with ROI name overlayed
for image_name in [name for name in image_names if 't_0' in name]:
    pixel_cyto = pixel_dict[image_name]['cyto']
    pixel_nuc = pixel_dict[image_name]['nuclei_pixels']
    summary_df = pd.DataFrame(summary_dict[image_name]['cyto'])

    fig, ax = plt.subplots(figsize=(10.24, 10.24))
    plt.scatter(pixel_nuc['X_pos'], pixel_nuc['Y_pos'], c=pixel_nuc['mCherry_nuclei'], cmap='Blues', s=1)

    plt.scatter(pixel_cyto['X_pos'], pixel_cyto['Y_pos'], c=pixel_cyto['mCherry_cells'], cmap='Reds', s=1)

    for roi_name in summary_df.index.tolist():
        plt.annotate(roi_name, (summary_df.loc[roi_name, 'X_pos'], summary_df.loc[roi_name, 'Y_pos']))
    plt.xlim(0, 1024)
    plt.ylim(0, 1024)

    plt.savefig(f'{output_path}{image_name}_image.png')

# Generate compiled summary_df for all images, after adding descriptive nuc_id to each ROI
summary_list = []
for image_name in image_names:
    summary_df = summary_dict[image_name]['nuclei']
    # Add additional info from image names
    summary_df['mutant'] = summary_df['image_name'].str.split('_').str[0]
    summary_df['position'] = summary_df['image_name'].str.split('_').str[2]
    summary_df['time'] = summary_df['image_name'].str.split('_').str[4]
    summary_df['number'] = [str(number) for number in range(0, summary_df.shape[0])]
    summary_df['nuc_id'] = list(zip(summary_df['mutant'], summary_df['position'], summary_df['time'], summary_df['number']))
    summary_df['nuc_id'] = ['_'.join(tuple) for tuple in summary_df['nuc_id']]
    summary_dict[image_name]['nuclei'] = summary_df
    summary_list.append(summary_df)
mean_nuclei_summary = pd.concat(summary_list)


# Save relevant df to excel, hdf5
# FileHandling.df_to_excel(output_path+f'mean_summary.xlsx', ['summary_sorted'], data_frames=[mean_nuclei_summary])


# Generate nuclei tracks according to maximum moved distance, assign new name
mutant_mapped = {}
for mutant in mutants:
    pos_list = list(set([name.split('_')[2] for name in image_names if mutant in name]))
    position_mapped = {}
    for pos in pos_list:

        # Get list of nuclei for cells in t0 (only keep these nuclei and track these through time)
        #---> There are some cells without nuclei in t0, some cells that are too close to many nuclei and some that have only one. Optimising the tolerance is necessary to best capture the relevant types.
        df_a = cyto_summary[f'{mutant}_pos_{pos}_t_0'][cyto_summary[f'{mutant}_pos_{pos}_t_0']['cotrans'] == 1].copy()
        # Filter only cotransfected cells
        df_a['coords'] = list(zip(round(df_a['X_pos'], 0), round(df_a['Y_pos'], 0)))

        df_b = summary_dict[f'{mutant}_pos_{pos}_t_0']['nuclei'].copy()
        df_b['coords'] = list(zip(round(df_b['X_pos'], 0), round(df_b['Y_pos'], 0)))

        t0_cell_nuclei = roi_matcher(df_a, df_b.set_index('nuc_id'), col='coords', x_tolerance=40, y_tolerance=40)

        timepoint_names = sorted_nicely(list(set(mean_nuclei_summary['time'])))
        logger.info(f'The following timepoints were detected:\n {timepoint_names}')
        #define empty df to store tracking info, set column 1 to be cell_ids from t0
        roi_mapped = pd.DataFrame()
        roi_mapped[timepoint_names[0]] = list(t0_cell_nuclei.values())
        for x, timepoint in enumerate(timepoint_names[:-1]):
            # select ROI in first timepoint --> find matching ROI in second timepoint
            # repeat up to second last timepoint i.e. x matched with x+1 up to (len(timepoints)-1)
            df_a = summary_dict[f'{mutant}_pos_{pos}_t_{timepoint_names[x]}']['nuclei']
            df_a['coords'] = list(zip(round(df_a['X_pos'], 0), round(df_a['Y_pos'], 0)))

            df_b = summary_dict[f'{mutant}_pos_{pos}_t_{timepoint_names[x+1]}']['nuclei']
            df_b['coords'] = list(zip(round(df_b['X_pos'], 0), round(df_b['Y_pos'], 0)))
            # find matching ROIs between two timepoints - can play with the tolerance here!
            matched_rois = roi_matcher(df_a.set_index('nuc_id'), df_b.set_index('nuc_id'), col='coords', x_tolerance=30, y_tolerance=30)
            #Add new column to the dataframe for new timepoint by mapping the previous timepoint with the matched rois
            roi_mapped[timepoint_names[x+1]] = roi_mapped[timepoint_names[x]].map(matched_rois)
        # add new column with updated nuclei id for a single track
        roi_mapped['track_id'] = [f'{mutant}_pos_{pos}_nuc_{x}' for x in range(0, len(roi_mapped))]
        position_mapped[pos] = roi_mapped
    mutant_mapped[mutant] = position_mapped

# generate a dictionary mapping the cell_id to cell_track_id
track_dict = {}
for mutant in mutants:
    pos_list = list(set([name.split('_')[2] for name in image_names if mutant in name]))
    position_mapped = mutant_mapped[mutant]
    for pos in pos_list:
        mapped_rois = position_mapped[pos]
        for timepoint in timepoint_names:
            track_dict.update(dict(zip(mapped_rois[timepoint], mapped_rois['track_id'])))
# remove nan key, update background key as unnecessary
track_dict.pop(np.nan, None)

mean_nuclei_summary['track_id'] = mean_nuclei_summary['nuc_id'].map(track_dict)
mean_nuclei_summary['time'] = mean_nuclei_summary['time'].astype('int')
summary_sorted = mean_nuclei_summary.sort_values(['position','track_id', 'time']).dropna()
FileHandling.df_to_excel(output_path+f'mean_nuclei_summary.xlsx', sheetnames=['summary_sorted'], data_frames=[summary_sorted.reset_index()])


# Generate grouped lists per mutant
summary_dfs = []
summary_names = []
for group, df in summary_sorted.groupby('mutant'):
    summary_dfs.append(df.reset_index())
    summary_names.append(group)
FileHandling.df_to_excel(output_path+f'summary_per_mutant.xlsx', sheetnames=summary_names, data_frames=summary_dfs)


plottable = summary_sorted[['track_id', 'time', 'Intensity']]
plottable = plottable.pivot(index='time', columns='track_id', values='Intensity').reset_index().sort_values('time', axis=0)
# Generate grouped lists per mutant
plottable_dfs = []
plottable_names = []
for mutant in mutants:
    plottable_dfs.append(plottable[[col for col in plottable.columns.tolist() if mutant in col]].reset_index())
    plottable_names.append(group)
FileHandling.df_to_excel(output_path+f'plottable_per_mutant.xlsx', plottable_names, data_frames=plottable_dfs)


# to plot image of cyto/nucleo pixels for first timepoint, with track_id labelled
for mutant in mutants:
    for pos in pos_list:
        image_name = f'{mutant}_pos_{pos}_t_0'
        pixel_cyto = pixel_dict[image_name]['cyto']
        pixel_nuc = pixel_dict[image_name]['nuclei_pixels']
        summary_df = pd.DataFrame(summary_dict[image_name]['cyto'])
        # generate list of ROIs to annotate
        label_df = summary_sorted[summary_sorted['time'] == int(timepoint_names[0])]
        label_df = label_df[['track_id', 'X_pos', 'Y_pos']].set_index('track_id')

        fig, ax = plt.subplots(figsize=(10.24, 10.24))
        plt.scatter(pixel_nuc['X_pos'], pixel_nuc['Y_pos'], c=pixel_nuc['mCherry_nuclei'], cmap='Blues', s=1)
        plt.scatter(pixel_cyto['X_pos'], pixel_cyto['Y_pos'], c=pixel_cyto['mCherry_cells'], cmap='Reds', s=1)

        for label in label_df.index.tolist():
            text_label = label.split('_')[-1]
            plt.annotate(text_label, (round(label_df.loc[label,'X_pos'], 0), round(label_df.loc[label,'Y_pos'], 0)))
        plt.xlim(0, 1024)
        plt.ylim(0, 1024)
        plt.title(image_name)
        plt.savefig(f'{output_path}{image_name}_trackID.png')

# To plot track coords over time
for mutant in mutants:
    image_df = summary_sorted[summary_sorted['mutant'] == mutant]
    for pos in pos_list:
        image_df = image_df[image_df['position'] == pos]
        fig = track_plotter(image_df, track_col='track_id')
        plt.title(f'{mutant}_pos_{pos}')
        plt.savefig(f'{output_path}{mutant}_pos_{pos}_track_movement.png')


# To plot individual intensities for each position over time
for x, df in enumerate(plottable_dfs):
    fig, ax = plt.subplots()
    plt.plot(df)
    plt.xlabel('Timepoint')
    plt.ylabel('Mean Nuclear mCherry Intensity (A.U.)')
    plt.title(plottable_names[x])
