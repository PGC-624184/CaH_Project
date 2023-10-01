import numpy as np
import pandas as pd
import astropy.units as u
import matplotlib.pyplot as plt

filename = "Ca_H2_output.csv"
df = pd.read_csv(filename, sep=',',header=None)
df.rename(columns={0:"Radius", 1:"Broadening"},inplace=True)
df.Broadening = df.Broadening*10
print(df)

filename2 = "Ca_H_output.csv"
df2 = pd.read_csv(filename2, sep=',',header=None)
df2.rename(columns={0:"Radius", 1:"Broadening"},inplace=True)
df2.Broadening = df2.Broadening*10

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

N_H2 = 1.0381111448319373e+19*u.cm**(-3)
N_H = 9.041570094762417e+17*u.cm**(-3)
p_H2 =[]
p_H = []
radial = df.Radius*u.angstrom
for val in radial:
    p_H2.append(interparticle_distribution(val,N_H2).value)
radial2 = df2.Radius*u.angstrom
for val in radial2:
    p_H.append(interparticle_distribution(val,N_H).value)




plt.clf()
plt.semilogy(df.Radius, df.Broadening*p_H2,'x', label=r"Ca $H_2$")
plt.semilogy(df2.Radius, df2.Broadening*p_H,'+', label="Ca H")
plt.xlabel(r"Interatomic Spacing $(\AA)$")
plt.ylabel(r"Quasi-Static Broadening $(\AA)$")
plt.title("Probability Weighted Quasi-Static Calcium Broadening\n for the Photosphere of Proxima Centauri")
plt.legend()
plt.show()

average = np.sum(broad*p).to(u.angstrom)
print(average)


