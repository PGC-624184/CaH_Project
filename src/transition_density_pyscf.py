import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from functools import reduce
from pyscf import gto, dft, tddft
from tqdm import tqdm


def run_computation(r: float,states: int):
    """
    Perform electronic structure calculations for molecules with varying bond separations. 
    
    This function calculates electronic energies for molecules with varying bond lengths. It performs a series of electronic structure calculations using quantum chemistry methods for each specified bond length in the range [start, stop) with the given step size. For each bond length, it constructs a molecule consisting of calcium ('Ca') and hydrogen ('H') atoms separated by the specified bond length and calculates the electronic energies of the system. 
    
    It's assumed that the molecular hydrogen is bound at the most likely length, but the dependence of this on how the system behaves has not been examined. Molecular orientation of the H_2 molecule is additionally ignored for the first sets of calculations

    The electronic structure calculations include:
        - Building a molecular configuration.
        - Solving the electronic Schrödinger equation using density functional theory (DFT) with the 'wB97X_V' exchange-correlation functional.
        - Performing time-dependent density functional theory (TDDFT) calculations to compute excited states.
        - Storing the electronic energies in the 'E_mat' NumPy array.

    The function returns two numpy arrays, one containing the energies (difference between the ground and excited state in atomic energy units) and another numpy array with the oscillator strength for the particular orbital.

    Args:
        - r (float): The bond separation length in angstrom.
        - states (int): The number of excited states to compute for the molecule.

    Returns:
        - Energies (np.ndarray): A 2D NumPy array containing electronic energies.
            Rows correspond to different bond lengths.
            Columns correspond to different electronic states.
        - Oscillator (np.ndarray): Oscillator strengths for the line profile.
     """
    mol = gto.M(atom=[  ['Ca', 0,     0,     0],
                    ['H', r,0, 0.0],
                    ],
                    basis='def2-QZVPP', # This is a high order DFT basis set that reproduces the expected transitions
                    spin = 1,
                    verbose = 3)
    # Reduced Hartree Fock solution as an initial case           
    mf = dft.UKS(mol)
    # set the functional exchange
    mf.xc = 'wB97X_V'
    mf.kernel()

    mytd = tddft.TDDFT(mf) # Time Dependant Density Functional Theory
    mytd.singlet = True # Calculate the singlet states (not forbidden transitions)
    mytd.nstates = states # the number of states to calculate
    results = mytd.kernel()
    mytd.analyze()
    oscillator = mytd.oscillator_strength(gauge='length') # Other option is velocity gauge, but this seems to be the right results
    return results[0], oscillator
    


## Example https://github.com/jamesETsmith/2022_simons_collab_pyscf_workshop/blob/main/demos/05_Excited_States.ipynb

if __name__=="__main__":

    # The range of interatomic distance between the Ca and the H2 molecule
    radius = np.arange(4,20,0.25)
    """
    The number of excited states to calculate for the system. States 6,7,8 correspond to the three 4P orbitals one of the 4s electrons can jump into. This corresponds to the first excited state, with the prior excited states relating to rotational/vibrational molecule states for the system and are not the ones we are particularly interested in (for the moment).
    """
    states = 15
    
    # Run the computation
    for r in tqdm(radius):
        Energies, oscillator = run_computation(r,states)
        filename = "data/Coarse_curve_data_Ca_H2_r{r}.csv".format(r=r)
        df_E = pd.DataFrame({'Energies': Energies, 'f':oscillator})
        df_E.to_csv(filename,sep='\t')
    