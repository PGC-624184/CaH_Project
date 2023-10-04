import numpy as np
import pandas as pd
import astropy.units as u
import matplotlib.pyplot as plt
import os

def read_data(filename):
    try:
        df = pd.read_csv(filename, sep=',',header=None)
        df.rename(columns={0:"Radius", 1:"Broadening"},inplace=True)
        df.Broadening = np.array(df.Broadening) #nm
        return df
    except FileNotFoundError:
        print("File name is wrong or does not exist!.")

def interparticle_distribution(r,N):
    """
    Calculate the interparticle distribution for a given system.

    This function computes the interparticle distribution for a system with a given
    radial distance 'r' and a specified number of particles 'N'. It applies a
    mathematical formula and returns the result in units of inverse angstroms.

    Args:
        r (astropy.Quantity): The radial distance at which to calculate the distribution.
        N (int): The number of particles in the system.

    Returns:
        astropy.Quantity: The interparticle distribution in inverse angstroms.

    Note:
        - This function is designed for a specific mathematical calculation and requires
          the 'r' parameter to be an astropy Quantity with units (e.g., u.angstrom).
        - Ensure 'r' is correctly defined with units before using this function.
    """
    return (N*4*np.pi*r**2*np.exp(-(4/3)*np.pi*r**3*N)).to(u.nm**(-1))


if __name__ == "__main__":
    N_H2 = 1.0381111448319373e+19*u.cm**(-3)
    N_H = 9.041570094762417e+17*u.cm**(-3)
    filename = "Ca_H2_output.csv"
    filename1 = "Ca_H_output.csv"
    df = read_data(filename)
    df2 = read_data(filename1)
    p_H2 =[]
    p_H = []
    radial = df.Radius*u.angstrom
    for val in radial:
        p_H2.append(interparticle_distribution(val,N_H2).value)
    radial2 = df2.Radius*u.angstrom
    for val in radial2:
        p_H.append(interparticle_distribution(val,N_H).value)
    
    average = np.sum(df.Broadening*p_H2)
    print("The average Ca H2 Broadening is:  ",average*u.nm)
    average1 = np.sum(df2.Broadening*p_H)
    print("The average Ca H Broadening is:  ",average1*u.nm)

    plt.clf()
    plt.semilogy(df.Radius, df.Broadening*p_H2,'x', label=r"Ca $H_2$")
    plt.semilogy(df2.Radius, df2.Broadening*p_H,'+', label=r"Ca $H$")
    plt.xlabel(r"Interatomic Spacing $(nm)$")
    plt.ylabel(r"Quasi-Static Broadening $(nm)$")
    plt.title("Probability Weighted Quasi-Static Calcium Broadening\n for the Photosphere of Proxima Centauri")
    plt.legend()
    plt.show()

    


