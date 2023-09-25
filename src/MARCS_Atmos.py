import numpy as np
import os
import pandas as pd

def import_NextGen_pressure(file,skiprows,skipfooter):
    """
Import a NextGen pressure data file and return it as a DataFrame.

This function reads a NextGen pressure data file in CSV format, with the option to
skip a specified number of rows at the beginning and end of the file. It uses the
'pandas' library to parse the file.

Args:
    file (str): The path to the CSV file to be imported.
    skiprows (int): The number of rows to skip at the beginning of the file.
    skipfooter (int): The number of rows to skip at the end of the file.

Returns:
    pandas.DataFrame: A DataFrame containing the data from the imported file.

Raises:
    FileNotFoundError: If the specified file is not found.

Example:
    >>> data_file = "nextgen_pressure_data.csv"
    >>> skipped_rows_start = 2
    >>> skipped_rows_end = 1
    >>> df = import_NextGen_pressure(data_file, skipped_rows_start, skipped_rows_end)

Note:
    - Make sure to have the 'pandas' library installed in your environment to use this function.
"""
    try:
        df = pd.read_csv(file,delim_whitespace=True, skiprows=skiprows,skipfooter=skipfooter,engine='python')
        return df
    except FileNotFoundError:
        print("The file might not be the right one you have specified!\n Double check and try again!")

#def import_Nextgen_abund(file, skiprows, skipfooter):

def import_MARCS_ppressure(file, skiprows, skipfooter):
    """
Import a NextGen pressure data file and return it as a DataFrame.

This function reads a NextGen pressure data file in CSV format, with the option to
skip a specified number of rows at the beginning and end of the file. It uses the
'pandas' library to parse the file.

Args:
    file (str): The path to the CSV file to be imported.
    skiprows (int): The number of rows to skip at the beginning of the file.
    skipfooter (int): The number of rows to skip at the end of the file.

Returns:
    pandas.DataFrame: A DataFrame containing the data from the imported file.

Raises:
    FileNotFoundError: If the specified file is not found.

Example:
    >>> data_file = "nextgen_pressure_data.csv"
    >>> skipped_rows_start = 2
    >>> skipped_rows_end = 1
    >>> df = import_NextGen_pressure(data_file, skipped_rows_start, skipped_rows_end)

Note:
    - Make sure to have the 'pandas' library installed in your environment to use this function.
"""
    try:
        cols = pd.read_csv(file,sep='   ', skiprows=skiprows-1, skipfooter=skiprows+55,engine='python').columns
        df = pd.read_csv(file, delim_whitespace=True,skipfooter=skipfooter, skiprows=skiprows,engine='python', header=None, names=cols)
        return df
    except FileNotFoundError:
        print("The file might not be the right one you have specified!\n Double check and try again!")

def import_model_struct(file, skiprows, skipfooter):
    """
    Import a NextGen pressure data file and return it as a DataFrame.

    This function reads a NextGen pressure data file in CSV format, with the option to
    skip a specified number of rows at the beginning and end of the file. It uses the
    'pandas' library to parse the file.

    Args:
        file (str): The path to the CSV file to be imported.
        skiprows (int): The number of rows to skip at the beginning of the file.
        skipfooter (int): The number of rows to skip at the end of the file.

    Returns:
        pandas.DataFrame: A DataFrame containing the data from the imported file.

    Raises:
        FileNotFoundError: If the specified file is not found.

    Example:
        >>> data_file = "nextgen_pressure_data.csv"
        >>> skipped_rows_start = 2
        >>> skipped_rows_end = 1
        >>> df = import_NextGen_pressure(data_file, skipped_rows_start, skipped_rows_end)

    Note:
        - Make sure to have the 'pandas' library installed in your environment to use this function.
    """
    try:
        cols = pd.read_csv(file, delim_whitespace=True, skiprows=skiprows-1, skipfooter=skipfooter+55, engine='python').columns
        df = pd.read_csv(file, delim_whitespace=True, skipfooter=skipfooter, skiprows=skiprows, engine='python', header=None, names=cols)
        df = df.drop(['lgTauR', 'lgTau5', 'Prad', 'Pturb','Pg'],axis=1)
        return df
    except FileNotFoundError:
        print("We seem to have a problem Houston!")

def import_MARCS(file):
    """
    Import MARCS model data from a file and return a combined DataFrame.

    This function reads data from a MARCS model file, which is organized into multiple
    tables. It imports three separate tables for pressure data and one table for model
    structure data from the specified file and combines them into a single DataFrame.

    Args:
        file (str): The path to the MARCS model file to be imported.

    Returns:
        pandas.DataFrame: A DataFrame containing the combined data from the MARCS model,
        including model structure and pressure data.

    Note:
        - This function relies on the 'import_MARCS_ppressure' and 'import_model_struct'
          functions for importing pressure and model structure data.
        - Ensure that these functions are properly defined and imported in your environment
          before using this function.
    """
    # First Table
    MARCS_pres_footer1 = 311-197
    MARCS_pres_skip1 = 140
    MARCS_pressure1 = import_MARCS_ppressure(file,MARCS_pres_skip1,MARCS_pres_footer1)
    MARCS_pressure1 = MARCS_pressure1.rename(columns={"k  lgPgas": "lgPgas"})

    # Second Table
    MARCS_pres_footer2 = 311-254
    MARCS_pres_skip2 = 197
    MARCS_pressure2 = import_MARCS_ppressure(file,MARCS_pres_skip2,MARCS_pres_footer2)
    MARCS_pressure2.index = MARCS_pressure2['k']
    MARCS_pressure2 = MARCS_pressure2.drop('k', axis=1)

    # Third Table
    MARCS_pres_footer3 = 0
    MARCS_pres_skip3 = 254
    MARCS_pressure3 = import_MARCS_ppressure(file,MARCS_pres_skip3,MARCS_pres_footer3)
    MARCS_pressure3.index = MARCS_pressure3['k']
    MARCS_pressure3 = MARCS_pressure3.drop('k',axis=1)

    # Model structure
    structure = import_model_struct(file, 25,229)
    structure.index = structure['k']
    structure = structure.drop('k', axis=1)

    # combine all the tables into one 
    result = pd.concat([structure, MARCS_pressure1, MARCS_pressure2, MARCS_pressure3], axis=1, join="inner")
    return result

def import_MARCS_abund(file,skiprows, skipfooter):
    """
    Import a NextGen pressure data file and return it as a DataFrame.

    This function reads a NextGen pressure data file in CSV format, with the option to
    skip a specified number of rows at the beginning and end of the file. It uses the
    'pandas' library to parse the file.

    Args:
        file (str): The path to the CSV file to be imported.
        skiprows (int): The number of rows to skip at the beginning of the file.
        skipfooter (int): The number of rows to skip at the end of the file.

    Returns:
        pandas.DataFrame: A DataFrame containing the data from the imported file.

    Raises:
        FileNotFoundError: If the specified file is not found.

    Example:
        >>> data_file = "nextgen_pressure_data.csv"
        >>> skipped_rows_start = 2
        >>> skipped_rows_end = 1
        >>> df = import_NextGen_pressure(data_file, skipped_rows_start, skipped_rows_end)

    Note:
        - Make sure to have the 'pandas' library installed in your environment to use this function.
    """
    try:
        df = pd.read_csv(file, delim_whitespace=True, skiprows=skiprows, skipfooter=skipfooter, engine='python', header=None)
        vals = 10**(np.concatenate(np.array(df))-12)
        return vals
    except FileNotFoundError:
        print("The file might not be the right one you have specified!\n Double check and try again!")



if __name__ == "__main__":
    # NextGen File
    NextGen_file = "reference_stellar_atmospheres/lte027-5.5-0.0a+0.0.BT-Settl.6"
    skiprows = 28854
    skipfooter = 30479-28984
    df_nextgen = import_NextGen_pressure(NextGen_file,skiprows,skipfooter)

    #MARCS partial pressures
    MARCS_file = "reference_stellar_atmospheres/t2700_g+5.5_z+0.00.mod"
    df = import_MARCS(MARCS_file)
    print(df)

