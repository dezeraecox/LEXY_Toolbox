import os
import re
from loguru import logger
from shutil import copyfile

logger.info('Import ok')

def jarvis(input_path, output_path):

    if not os.path.exists(output_path):
        os.mkdir(output_path)

    for filename in os.listdir(input_path):
        pathname = os.path.join(input_path, filename)
        if os.path.isdir(pathname):
            logger.info(f"Directory not processed {filename}")
        if '_pixels.csv' in filename:
            details = re.split('_', filename)
            ## Adjust channel list here
            mutant, cotrans, _, pos_number, _, time_point, image_type, channel, roi_type, _ = details
            # logger.debug(details)
            if channel == 'mCherry':
                if image_type == 'Before':
                    new_output = output_path + f'{mutant}+{cotrans}_'
                    new_name = f'pos_{pos_number}_t_{time_point}_{channel}_{roi_type}_pixels.csv'
                    # logger.debug(new_name)
                elif image_type == 'Activation':
                    time_point = str(int(time_point) + 1)
                    new_output = output_path + f'{mutant}+{cotrans}_'
                    new_name = f'pos_{pos_number}_t_{time_point}_{channel}_{roi_type}_pixels.csv'
                    # logger.debug(new_name)
                else:
                    logger.info(f'No suitable image type found. Image type listed as {image_type}')

            elif channel == 'GFP':
                if image_type == 'Activation':
                    time_point = str(int(time_point) + 1)
                    new_output = output_path + f'{mutant}+{cotrans}_'
                    new_name = f'pos_{pos_number}_t_{time_point}_{channel}_{roi_type}_pixels.csv'
                    # logger.debug(new_name)
            else:
                logger.info(f'No suitable channel information found. Channel listed as {channel}')
            copyfile(pathname, new_output+new_name)
            logger.info(f'File renamed to {new_name}')
        else:
            logger.info(f"{filename} not a pixel result, therefore not processed.")


## Adjust input path here
jarvis(input_path="", output_path='')
