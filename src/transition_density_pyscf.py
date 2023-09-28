import numpy as np
from pyscf import gto, scf, mcscf
from pyscf import gto, mp, mcscf, dft, tddft
from pyscf.mcscf import dmet_cas
from pyscf import gto, scf, cc
import matplotlib.pyplot as plt

# Bond seperation lengths in angstrom
bond = np.arange(1.0, 10.0, 0.05)
dm_init = None

radius = bond # distance to atomic units from angstrom

# Convert the below to a difference from the far field ground state
# then convert result to eV
E_h = 27.211386246012257


""" 
More working - finding the wavelength at each of the jump energy values.
This should be the normal transition level out at the far field, and then
get funkier the closer you are.

Rest of the code is WIP and so results may not be real
"""
n = 28.437353657002053
from functools import reduce
import numpy
from pyscf import gto, scf, dft, tddft

for r in reversed(bond):
    mol = gto.M(atom=[  ['Ca', 0,     0,     0],
                    ['H', n,0,0.3865],
                    ['H', n,0,-0.3865],
                    ],
                    basis='def2-QZVPP',
                    spin = 0,
                    verbose = 4)
    # Reduced Hartree Fock solution as an initial case           
    mf = dft.RKS(mol)
    mf.xc = 'wB97X_V'
    mf.kernel()

    mytd = tddft.TDDFT(mf)
    mytd.singlet = False
    mytd.nstates = 15
    mytd.kernel()
    mytd.analyze()

    # Need to collect energy values and then convert to wavelength
    
    # store results in array


## Example https://github.com/jamesETsmith/2022_simons_collab_pyscf_workshop/blob/main/demos/05_Excited_States.ipynb
from pyscf import gto, scf, dft, tddft,cc
import py3Dmol

# import data into dataframe - check if pd is row or col major
df = pd.DataFrame(data)


# plot the curves
import plotly.express as px
import nbformat
fig = px.line(df, x="Excitation Energy (eV)", y="Intensity", markers=True, color="Exchange-Correlation Functional")
fig.show()
