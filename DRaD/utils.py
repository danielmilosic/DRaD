import sys
import os
import warnings
from erfa import ErfaWarning
from astropy.utils.exceptions import ErfaWarning
import pandas as pd
import logging

def suppress_output(func, *args, **kwargs):
    """
    Executes a function while suppressing standard output, standard error, warnings, and logging messages.
    
    Args:
        func (callable): The function to execute.
        *args: Positional arguments to pass to the function.
        **kwargs: Keyword arguments to pass to the function.
    
    Returns:
        The result of the executed function.
    """
    # Save original stdout, stderr, and logging level
    original_stdout = sys.stdout
    original_stderr = sys.stderr
    original_logging_level = logging.getLogger().level  # Get current logging level
    warnings.filterwarnings("ignore", category=ErfaWarning)
    try:
        # Redirect stdout and stderr to devnull
        sys.stdout = open(os.devnull, 'w')
        sys.stderr = open(os.devnull, 'w')

        # Suppress warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ErfaWarning)
            warnings.simplefilter("ignore", RuntimeWarning)
            warnings.simplefilter("ignore")  # Ignore all warnings

            # Suppress logging
            logging.getLogger().setLevel(logging.CRITICAL)  # Only show critical errors

            # Run the function
            result = func(*args, **kwargs)

    finally:
        # Restore stdout, stderr, and logging level
        sys.stdout = original_stdout
        sys.stderr = original_stderr
        logging.getLogger().setLevel(original_logging_level)  # Restore logging level

    return result

def spacecraft_ID(ID, ID_number=False):
    """
    Retrieve the spacecraft ID or name based on the provided DataFrame or number.

    Parameters:
        ID (DataFrame or int): DataFrame containing 'Spacecraft_ID' or an ID number.
        ID_number (bool): Whether to return the ID number instead of the name.

    Returns:
        str or int: Spacecraft name or ID number.
    """
    df = pd.DataFrame({
        'ID': ['PSP', 'SolO', 'BepiC', 'STEREO-A', 'STEREO-B', 'OMNI', 'MAVEN']
    }, index=[1, 2, 3, 4, 5, 6, 7])

    if isinstance(ID, pd.DataFrame):
        # Check if DataFrame is empty
        if ID.empty:
            raise ValueError("The input DataFrame is empty. Cannot determine spacecraft ID.")
        
        # Check for 'Spacecraft_ID' column
        if 'Spacecraft_ID' in ID.columns:
            ID.dropna(subset=['Spacecraft_ID'], inplace=True)  # Drop rows where 'Spacecraft_ID' is NaN
            if ID.empty:  # Check again if DataFrame is empty after dropping NaN
                raise ValueError("The 'Spacecraft_ID' column is empty after dropping NaN values.")
            
            # Safely access the first value
            number = ID.iloc[0]['Spacecraft_ID']
        else:
            raise ValueError("DataFrame must contain a 'Spacecraft_ID' column.")
    elif isinstance(ID, (int, float)):
        number = int(ID)
    elif isinstance(ID, str):
        match = df[df['ID'].str.upper() == ID.upper()]
        if not match.empty:
            number = match.index[0]
    else:
        raise TypeError("ID must be a DataFrame, int, or float.")

    # Return ID or name
    if ID_number:
        return number
    else:
        return df.loc[number, 'ID']

from PIL import Image


def glue_together(folder_a, folder_b, output_folder, mode="horizontal"):
    """
    Combines images from two folders by placing them side-by-side or stacked vertically,
    then saves the resulting images in the specified output folder.
    
    Args:
        folder_a (str): Path to the first folder containing images.
        folder_b (str): Path to the second folder containing images.
        output_folder (str): Path to the folder where combined images will be saved.
        mode (str): "horizontal" to place images side-by-side, "vertical" to stack them.

    Notes:
        - Assumes both folders contain PNG images.
        - If the number of images in folder_a and folder_b do not match, a warning is displayed.
        - If dimensions differ, images from folder_b are resized to match folder_a (width for vertical, height for horizontal).
    """

    # Create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    # Get sorted lists of file names from both folders
    files_a = sorted(f for f in os.listdir(folder_a) if f.endswith(".png"))
    files_b = sorted(f for f in os.listdir(folder_b) if f.endswith(".png"))

    # Ensure both folders have the same number of files
    if len(files_a) != len(files_b):
        print("Warning: The number of files in folder A and folder B do not match!")

    # Process files
    for file_a, file_b in zip(files_a, files_b):
        # Open images
        img_a = Image.open(os.path.join(folder_a, file_a))
        img_b = Image.open(os.path.join(folder_b, file_b))

        if mode == "horizontal":
            # Ensure images have the same height
            if img_a.height != img_b.height:
                img_b = img_b.resize((img_b.width, img_a.height), Image.ANTIALIAS)
            
            # Create a new image with combined width
            combined_size = (img_a.width + img_b.width, img_a.height)
            img_b_offset = (img_a.width, 0)  # Paste img_b to the right

        elif mode == "vertical":
            # Ensure images have the same width
            if img_a.width != img_b.width:
                img_b = img_b.resize((img_a.width, img_b.height), Image.ANTIALIAS)
            
            # Create a new image with combined height
            combined_size = (img_a.width, img_a.height + img_b.height)
            img_b_offset = (0, img_a.height)  # Paste img_b below
        
        else:
            raise ValueError("Invalid mode. Choose 'horizontal' or 'vertical'.")

        # Create combined image
        combined_image = Image.new("RGB", combined_size)
        combined_image.paste(img_a, (0, 0))
        combined_image.paste(img_b, img_b_offset)

        # Save combined image
        output_path = os.path.join(output_folder, file_a)
        combined_image.save(output_path)

        print(f"Saved combined image: {output_path}")

    print("All images have been processed and saved.")

import numpy as np
import pandas as pd

def pad_data_with_nans(data, before, after, cadence='H'):
    """
    Add NaNs to the start and end of the DataFrame to make it have a consistent length,
    and linearly extrapolate the `CARR_LON` column.

    Parameters:
    - data (pd.DataFrame): DataFrame to pad and extrapolate.
    - before (str): Start datetime as a string (e.g., 'YYYY-MM-DD HH:MM').
    - after (str): End datetime as a string (e.g., 'YYYY-MM-DD HH:MM').

    Returns:
    - padded_data (pd.DataFrame): Padded DataFrame with extrapolated `CARR_LON`.
    """

    # Ensure the index is a DateTimeIndex
    data = data.copy()
    if not isinstance(data.index, pd.DatetimeIndex):
        raise ValueError("The DataFrame index must be a DatetimeIndex.")

    # Create padding for before
    pad_before_index = pd.date_range(start=before, end=data.index.min(), freq=cadence)
    pad_before = pd.DataFrame(index=pad_before_index, columns=data.columns)
    pad_before[:] = np.nan

    # Create padding for after
    pad_after_index = pd.date_range(start=data.index.max(), end=after, freq=cadence)
    pad_after = pd.DataFrame(index=pad_after_index, columns=data.columns)
    pad_after[:] = np.nan

    # Concatenate padding and data
    padded_df = pd.concat([pad_before, data, pad_after])

    # Linearly extrapolate the `CARR_LON` column
    if 'CARR_LON' in data.columns:
        carr_lon = data['CARR_LON']
        slope = (carr_lon.iloc[-1] - carr_lon.iloc[0]) / (carr_lon.index[-1] - carr_lon.index[0]).total_seconds()

        # Extrapolate for padding before
        for idx in pad_before.index:
            delta = (idx - carr_lon.index[0]).total_seconds()
            pad_before.loc[idx, 'CARR_LON'] = carr_lon.iloc[0] + slope * delta

        # Extrapolate for padding after
        for idx in pad_after.index:
            delta = (idx - carr_lon.index[-1]).total_seconds()
            pad_after.loc[idx, 'CARR_LON'] = carr_lon.iloc[-1] + slope * delta

        # Update the padded DataFrame
        padded_df['CARR_LON'] = pd.concat([pad_before['CARR_LON'], carr_lon, pad_after['CARR_LON']])

    return padded_df

def plot_timeseries(timeseries_data):
    import matplotlib.pyplot as plt
    import seaborn as sns

    timeseries_data = timeseries_data.dropna(axis=1, how='all')

    fig, axes = plt.subplots(nrows=8, ncols=1, figsize=(10, 12), sharex=True)

    # Panel 0: V and V_R
    if 'V' in timeseries_data.columns:
        sns.lineplot(data=timeseries_data, x=timeseries_data.index, y='V', ax=axes[0], color='black')
    else:
        axes[0].text(0.98, 0.7, "Missing: V", transform=axes[0].transAxes,
                     ha='right', va='center', color='black', fontsize=9)
    if 'V_R' in timeseries_data.columns:
        v_data = timeseries_data[timeseries_data['V_R']>-1000]
        sns.lineplot(data=v_data, x=v_data.index, y='V_R', ax=axes[0], color='blue', alpha=0.5)
    else:
        axes[0].text(0.98, 0.5, "Missing: V_R", transform=axes[0].transAxes,
                     ha='right', va='center', color='black', fontsize=9)
    axes[0].set_ylabel('V $[km s^{-1}]$')

    # Panel 1: V_T and V_N
    if 'V_T' in timeseries_data.columns:
        sns.lineplot(data=timeseries_data, x=timeseries_data.index, y='V_T', ax=axes[1], label='V_T')
    else:
        axes[1].text(0.98, 0.7, "Missing: V_T", transform=axes[1].transAxes,
                     ha='right', va='center', color='black', fontsize=9)
    if 'V_N' in timeseries_data.columns:
        sns.lineplot(data=timeseries_data, x=timeseries_data.index, y='V_N', ax=axes[1], label='V_N')
    else:
        axes[1].text(0.98, 0.5, "Missing: V_N", transform=axes[1].transAxes,
                     ha='right', va='center', color='black', fontsize=9)
    axes[1].set_ylabel('V_TN $[km s^{-1}]$')

    # Panel 2: N and T (right y-axis)
    if 'N' in timeseries_data.columns:
        sns.lineplot(data=timeseries_data, x=timeseries_data.index, y='N', ax=axes[2], color='black')
    else:
        axes[2].text(0.98, 0.7, "Missing: N", transform=axes[2].transAxes,
                     ha='right', va='center', color='black', fontsize=9)
    axes[2].set_ylabel('N $[cm^{-3}]$')

    if 'T' in timeseries_data.columns:
        ax2 = axes[2].twinx()
        sns.lineplot(data=timeseries_data, x=timeseries_data.index, y='T', ax=ax2, color='tab:blue')
        ax2.set_ylabel('T $[K]$')
        ax2.spines['right'].set_color('tab:blue')
        ax2.yaxis.label.set_color('tab:blue')
        ax2.tick_params(axis='y', colors='tab:blue')
    else:
        axes[2].text(0.98, 0.3, "Missing: T", transform=axes[2].transAxes,
                     ha='right', va='center', color='tab:blue', fontsize=9)

    # Panel 3: P and Beta (right y-axis)
    if 'P' in timeseries_data.columns:
        sns.lineplot(data=timeseries_data, x=timeseries_data.index, y='P', ax=axes[3], color='black')
    else:
        axes[3].text(0.98, 0.7, "Missing: P", transform=axes[3].transAxes,
                     ha='right', va='center', color='black', fontsize=9)
    axes[3].set_ylabel('P $[nPa]$')

    if 'Beta' in timeseries_data.columns:
        ax3 = axes[3].twinx()
        sns.lineplot(data=timeseries_data, x=timeseries_data.index, y='Beta', ax=ax3, color='tab:blue')
        ax3.set_ylabel(r'$\beta$')
        ax3.set_yscale('log')
        ax3.spines['right'].set_color('tab:blue')
        ax3.yaxis.label.set_color('tab:blue')
        ax3.tick_params(axis='y', colors='tab:blue')
    else:
        axes[3].text(0.98, 0.3, "Missing: Beta", transform=axes[3].transAxes,
                     ha='right', va='center', color='tab:blue', fontsize=9)

    # Panel 4: B with polarity
    if 'B' in timeseries_data.columns and 'POL' in timeseries_data.columns:
        timeseries_data['polarity'] = ['+' if pol > 0 else '-' for pol in timeseries_data['POL']]
        colors = {'-': 'tab:blue', '+': 'tab:red'}
        sns.scatterplot(data=timeseries_data, x=timeseries_data.index, y='B',
                        hue='polarity', palette=colors, ax=axes[4], s=5, alpha=1)
        axes[4].get_legend().remove()
    else:
        if 'B' not in timeseries_data.columns:
            axes[4].text(0.98, 0.7, "Missing: B", transform=axes[4].transAxes,
                         ha='right', va='center', color='black', fontsize=9)
        if 'POL' not in timeseries_data.columns:
            axes[4].text(0.98, 0.5, "Missing: POL", transform=axes[4].transAxes,
                         ha='right', va='center', color='black', fontsize=9)
    axes[4].set_ylabel('B $[nT]$')

    # Panel 5: B_R, B_T, B_N
    for component, color, ypos in zip(['B_R', 'B_T', 'B_N'], ['red', 'green', 'blue'], [0.7, 0.5, 0.3]):
        if component in timeseries_data.columns:
            sns.lineplot(data=timeseries_data, x=timeseries_data.index, y=component, ax=axes[5], color=color, label=component)
        else:
            axes[5].text(0.98, ypos, f"Missing: {component}", transform=axes[5].transAxes,
                         ha='right', va='center', color=color, fontsize=9)
    axes[5].set_ylabel('B $[nT]$')

    # Panel 6: S_P and O7_O6_RATIO (right y-axis)
    if 'S_P' in timeseries_data.columns:
        sns.lineplot(data=timeseries_data, x=timeseries_data.index, y='S_P', ax=axes[6], color='black')
        axes[6].fill_between(timeseries_data.index, 2.69, 4, color='grey', alpha=0.7)
    else:
        axes[6].text(0.98, 0.7, "Missing: S_P", transform=axes[6].transAxes,
                     ha='right', va='center', color='black', fontsize=9)
    axes[6].set_ylabel('$S_p$ $[eV cm^{2}]$')

    if 'O7_O6_RATIO' in timeseries_data.columns:
        ax6 = axes[6].twinx()
        sns.scatterplot(data=timeseries_data, x=timeseries_data.index, y='O7_O6_RATIO', ax=ax6, s=5, color='tab:blue')
        sns.lineplot(data=timeseries_data, x=timeseries_data.index,
                     y=[0.145]*len(timeseries_data), ax=ax6, color='grey')
        ax6.set_ylim([0, 0.3])
        ax6.set_ylabel('$O^{7+}/O^{6+}$')
        ax6.spines['right'].set_color('tab:blue')
        ax6.yaxis.label.set_color('tab:blue')
        ax6.tick_params(axis='y', colors='tab:blue')
    else:
        axes[6].text(0.98, 0.3, "Missing: O7_O6_RATIO", transform=axes[6].transAxes,
                     ha='right', va='center', color='tab:blue', fontsize=9)

    # Panel 7: R and LAT (right y-axis)
    if 'R' in timeseries_data.columns:
        sns.lineplot(data=timeseries_data, x=timeseries_data.index, y='R', ax=axes[7], color='black')
    else:
        axes[7].text(0.98, 0.7, "Missing: R", transform=axes[7].transAxes,
                     ha='right', va='center', color='black', fontsize=9)
    axes[7].set_ylabel('r $[AU]$')

    if 'LAT' in timeseries_data.columns:
        ax7 = axes[7].twinx()
        sns.lineplot(data=timeseries_data, x=timeseries_data.index, y='LAT', ax=ax7, color='tab:blue')
        ax7.set_ylabel('LAT $[°]$')
        ax7.spines['right'].set_color('tab:blue')
        ax7.yaxis.label.set_color('tab:blue')
        ax7.tick_params(axis='y', colors='tab:blue')
    else:
        axes[7].text(0.98, 0.3, "Missing: LAT", transform=axes[7].transAxes,
                     ha='right', va='center', color='tab:blue', fontsize=9)

    plt.tight_layout(pad=1., w_pad=0.5, h_pad=.1)
