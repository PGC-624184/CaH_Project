import numpy as np
import os
import pandas as pd
import astropy.units as u
import astropy.constants.codata2018 as c
import matplotlib.pyplot as plt
from Prob_dist_QM import interparticle_distribution

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
Import a MARCS data file and return partial pressures as a DataFrame.

This function reads a MARCS data file in CSV format, with the option to
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
    >>> data_file = "MARCS_data.csv"
    >>> skipped_rows_start = 2
    >>> skipped_rows_end = 1
    >>> df = import_MARCS_pressure(data_file, skipped_rows_start, skipped_rows_end)

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
    Import a MARCS data file and return the model structure as a DataFrame.

    This function reads a MARCS data file in CSV format, with the option to
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
        >>> data_file = "MARCS_data.csv"
        >>> skipped_rows_start = 2
        >>> skipped_rows_end = 1
        >>> df = import_model_structure(data_file, skipped_rows_start, skipped_rows_end)

    Note:
        - Make sure to have the 'pandas' library installed in your environment to use this function.
    """
    try:
        cols = pd.read_csv(file, delim_whitespace=True, skiprows=skiprows-1, skipfooter=skipfooter+55, engine='python').columns
        df = pd.read_csv(file, delim_whitespace=True, skipfooter=skipfooter, skiprows=skiprows, engine='python', header=None, names=cols)
        df = df.drop(['lgTau5', 'Prad', 'Pturb','Pg'],axis=1)
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
    MARCS_pressure1.rename(columns=lambda x: x.strip(), inplace=True)

    # Second Table
    MARCS_pres_footer2 = 311-254
    MARCS_pres_skip2 = 197
    MARCS_pressure2 = import_MARCS_ppressure(file,MARCS_pres_skip2,MARCS_pres_footer2)
    MARCS_pressure2.index = MARCS_pressure2['k']
    MARCS_pressure2 = MARCS_pressure2.drop('k', axis=1)
    MARCS_pressure2.rename(columns=lambda x: x.strip(), inplace=True)

    # Third Table
    MARCS_pres_footer3 = 0
    MARCS_pres_skip3 = 254
    MARCS_pressure3 = import_MARCS_ppressure(file,MARCS_pres_skip3,MARCS_pres_footer3)
    MARCS_pressure3.index = MARCS_pressure3['k']
    MARCS_pressure3 = MARCS_pressure3.drop('k',axis=1)
    MARCS_pressure3.rename(columns=lambda x: x.strip(), inplace=True)

    # Model structure
    structure = import_model_struct(file, 25,229)
    structure.index = structure['k']
    structure.drop('k', axis=1,inplace=True)
    structure.rename(columns=lambda x: x.strip(), inplace=True)
    structure.rename(columns={'T':'Temp'}, inplace=True)

    # combine all the tables into one 
    result = pd.concat([structure, MARCS_pressure1, MARCS_pressure2, MARCS_pressure3], axis=1, join="inner")
    return result


def number_density(dataframe,species,abund):
    """
    Calculate number density of a specified species in a given DataFrame.

    This function calculates the number density of a specified chemical species at each
    row of a DataFrame. It uses the dataframe's 'species' column to access the abundance
    of the species. If the 'species' column is not found, it falls back to using the
    'H I' column and an abundance dictionary provided as 'abund'.

    Args:
        dataframe (pandas.DataFrame): The DataFrame containing chemical data.
        species (str): The name of the chemical species for which to calculate the number density.
        abund (dict): A dictionary containing abundance values for different species.

    Returns:
        numpy.ndarray: An array of number densities (in cm^(-3)) corresponding to each row in the DataFrame.

    Note:
        - The DataFrame is expected to have columns 'species', 'H I', and 'Temp' for calculations.
        - The function utilizes physical constants from the 'u' and 'c' objects, so ensure
          these constants are properly defined in your environment.
        - Make sure the DataFrame and abundance dictionary are correctly formatted before
          using this function.
    """
    n = np.zeros(len(dataframe.index))
    try:
        for i in range(dataframe.index[0],len(dataframe.index)):
            n[i] = (((10**dataframe[species][i])*(u.dyne/(u.cm**2)))/(c.k_B*dataframe.Temp[i]*u.K)).to(u.cm**(-3)).value
        return n * u.cm**(-3)
    except KeyError:
        for i in range(dataframe.index[0],len(dataframe.index)):
            n[i] = (((10**(dataframe['H I'][i]+abund[species].values))*(u.dyne/(u.cm**2)))/(c.k_B*dataframe.Temp[i]*u.K)).to(u.cm**(-3)).value
        return n * u.cm**(-3)

def create_abundances(file):
    """
    Create a DataFrame of elemental abundances from a file.

    This function reads elemental abundance data from a file, where each element's
    abundance is represented as a column in the file. It then creates a DataFrame
    with the elemental abundances as columns, using a predefined list of atomic symbols.

    Args:
        file (str): The path to the file containing elemental abundance data.

    Returns:
        pandas.DataFrame: A DataFrame with columns representing elemental abundances
        for a range of chemical elements.

    Note:
        - The function expects the file to have a specific format, with abundances
          in columns corresponding to atomic symbols.
        - The 'atoms' list defines the order of atomic symbols used to create DataFrame columns.
        - Ensure the file format and 'atoms' list are correctly configured for your data.
    """
    atoms = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
             'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
             'K', 'Ca', 'Sc', 'Ti', 'V','Cr', 'Mn', 'Fe','Co', 'Ni',
             'Cu', 'Zn', 'Ga', 'Ge', 'As','Se', 'Br', 'Kr','Rb','Sr',
             'Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
             'Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm',
             'Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta',
             'W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U']
    df = pd.read_csv(file,skiprows=12,header=None,skipfooter=311-23,engine='python',delim_whitespace=True)
    tmp = np.reshape(np.array(df),(100,))
    tmp = tmp[~np.isnan(tmp)]
    new_df = pd.DataFrame(data=[tmp-12],columns=atoms)
    return(new_df)

def mean_interatomic_distance(df,species,abund):
    """
    Calculate the mean interatomic distance for a given species in a DataFrame.

    This function calculates the mean interatomic distance for a specified chemical
    species using the number density, which is obtained from the 'number_density'
    function. It calculates the mean interatomic distance based on the species'
    number density and returns the result in angstroms.

    Args:
        df (pandas.DataFrame): The DataFrame containing chemical data.
        species (str): The name of the chemical species for which to calculate the mean interatomic distance.
        abund (dict): A dictionary containing abundance values for different species.

    Returns:
        astropy.Quantity: The mean interatomic distance in angstroms.

    Note:
        - The 'number_density' function is used to obtain the number density necessary
          for this calculation.
        - The function relies on physical constants defined in the 'u' object, so make
          sure these constants are properly defined in your environment.
        - Ensure the DataFrame and abundance dictionary are correctly formatted before
          using this function.
    """
    n = number_density(df,species,abund)
    r_0 = ((4*np.pi*n/3)**(-1/3))
    return (r_0).to(u.angstrom)

if __name__ == "__main__":
    # NextGen File - WIP and mostly ignore for now
    """
    NextGen_file = "reference_stellar_atmospheres/lte027-5.5-0.0a+0.0.BT-Settl.6"
    skiprows = 28854
    skipfooter = 30479-28984
    df_nextgen = import_NextGen_pressure(NextGen_file,skiprows,skipfooter)
    """

    #MARCS partial pressures
    MARCS_file = "reference_stellar_atmospheres/t2700_g+5.5_z+0.00.mod"
    df = import_MARCS(MARCS_file)
    abund = create_abundances(MARCS_file)
    H2 = mean_interatomic_distance(df,'H2', abund)
    H = mean_interatomic_distance(df,'H', abund)
    H2_number_density = number_density(df,'H2',abund)
    H_number_density = number_density(df,'H',abund)
    H2 = mean_interatomic_distance(df, 'H2',abund)
    print(H_number_density)
    #print(H[41])
    #print(H2_number_density)
    #print(H2[41])
    print(np.array(df.lgTauR))

   

