import os
import pandas as pd

def get_plate_paths(folder_path):
    """Gets plate data paths
    Returns:
        list: list of plate data paths
    """
    plate_files = os.listdir(path=folder_path)

    #Need to make sure we are only attempting to load .csv or .xlsx data
    plate_files = [i for i in plate_files]

    plate_files.sort()

    plate_data_paths = []

    for pf in plate_files:
        if pf != '.DS_Store':
            plate_path = folder_path + os.sep + pf
            plate_data_paths.append(plate_path)

    return plate_data_paths

def get_data_file_paths(plate_path):
    files = os.listdir(path=plate_path)

    #Need to make sure we are only attempting to load .csv or .xlsx data
    files = [i for i in files if ('.csv' in i) or ('.xlsx' in i)]

    files.sort()

    file_data_paths = []

    for pf in files:
        if pf != '.DS_Store':
            file_path = plate_path + os.sep + pf
            file_data_paths.append(plate_path)

    return file_data_paths

folder_path = '/Users/eshanking/repos/seascapes_figures/data/08312022'

plate_paths = get_plate_paths(folder_path)

data_paths0 = get_data_file_paths(plate_paths[0])