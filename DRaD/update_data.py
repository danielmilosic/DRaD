import pandas as pd
from pathlib import Path
from datetime import datetime
from DRaD.utils import spacecraft_ID
from DRaD import psp, solo, stereo_a, omni, maven

def update(timerange=None, cadence='0.1H', spacecraft='all'):
    """
    Update and download reduced data for specified spacecraft.
    """
    if spacecraft == 'all':
        spacecraft_list = ['psp', 'solo', 'stereo_a', 'omni', 'maven']
    elif isinstance(spacecraft, str):
        spacecraft_list = [normalize_spacecraft_name(spacecraft_ID(spacecraft))]
    elif isinstance(spacecraft, (list, tuple)):
        spacecraft_list = [normalize_spacecraft_name(spacecraft_ID(sc)) for sc in spacecraft]
    else:
        raise TypeError("Invalid spacecraft input.")

    for sc in spacecraft_list:
        if timerange is None:
            last_month = latest_existing_month(sc)
            now = pd.Timestamp.now()
            timerange = (last_month, now)

        start, end = normalize_timerange(timerange)
        months = pd.date_range(start=start, end=end, freq='MS')
        months = [d.strftime('%Y-%m') for d in months]

        if sc == 'psp':
            psp.download_reduce_save_space(month=months, cadence=cadence)
        elif sc == 'solo':
            solo.download_reduce_save_space(month=months, cadence=cadence)
        elif sc == 'stereo_a':
            stereo_a.download_reduce_save_space(month=months, cadence=cadence)
        elif sc == 'omni':
            omni.download_reduce_save_space(month=months, cadence=cadence)
        elif sc == 'maven':
            maven.download_reduce_save_space(month=months, cadence=cadence)
        else:
            print(f"Unknown spacecraft: {sc}")

def normalize_timerange(timerange):
    """
    Normalize timerange input to (start_date, end_date) as datetime objects.
    Accepts:
        - single string (e.g., '2020-05' or '2020-05-15')
        - single datetime or pd.Timestamp
        - tuple/list of 2 elements
    """
    if timerange is None:
        return None
    if isinstance(timerange, (str, datetime, pd.Timestamp)):
        start = pd.to_datetime(timerange)
        end = increment_month(start) - pd.Timedelta(seconds=1)
    elif isinstance(timerange, (list, tuple)) and len(timerange) == 2:
        start = pd.to_datetime(timerange[0])
        end = pd.to_datetime(timerange[1])
        if start > end:
            raise ValueError("Start of timerange must be before end.")
    else:
        raise TypeError("Invalid timerange format. Use str, datetime, or tuple of two dates.")
    return start, end

def increment_month(date):
    """
    Return the first day of the next month.
    """
    if date.month == 12:
        return date.replace(year=date.year + 1, month=1, day=1)
    else:
        return date.replace(month=date.month + 1, day=1)

def find_existing_file(spacecraft, date):
    base_dirs = [Path(r'C:/Users/14milosi/DRaD/reduced_data'), Path('reduced_data')]
    for base in base_dirs:
        path = get_file_path(spacecraft, date, base)
        if path.exists():
            return path
    return None

def get_file_path(spacecraft, date, base_dir):
    year = date.year
    month = f"{date.month:02d}"
    filename = f"{spacecraft}_data{year}-{month}.parquet"
    return base_dir / spacecraft / filename

def latest_existing_month(spacecraft):
    """
    Find the latest available month for a given spacecraft in local storage.
    """
    base_dirs = [Path(r'C:/Users/14milosi/DRaD/reduced_data'), Path('reduced_data')]
    for base in base_dirs:
        folder = base / spacecraft
        if folder.exists():
            parquet_files = sorted(folder.glob(f"{spacecraft}_data*.parquet"))
            if parquet_files:
                latest_file = parquet_files[-1]
                date_str = latest_file.stem.split('data')[-1]
                return pd.to_datetime(date_str)
    print('NO SPACECRAFT DATA FOUND; RETURNING 2014-01-01')
    return pd.Timestamp('2014-01-01')  # Fallback start date if nothing is found

def normalize_spacecraft_name(name):
    """
    Convert official spacecraft name to file system-friendly name.
    E.g., 'STEREO-A' -> 'stereo_a'
    """
    mapping = {
        'PSP': 'psp',
        'SOLO': 'solo',
        'BEPIC': 'bepic',
        'STEREO-A': 'stereo_a',
        'STEREO-B': 'stereo_b',
        'OMNI': 'omni',
        'MAVEN': 'maven'
    }
    name_upper = name.upper()
    for official, fs_name in mapping.items():
        if name_upper == official:
            return fs_name
    raise ValueError(f"Unrecognized spacecraft name: {name}")

