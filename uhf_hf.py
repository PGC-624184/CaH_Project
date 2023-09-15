'''
Unrestricted Hartree-Fock using PySCF intergrals
Author: Hung Q. Pham (UMN, phamx494@umn.edu)
'''

import numpy as np
import scipy as sp
from pyscf import gto, scf, ao2mo
from functools import reduce

mol = gto.M()
mol.atom =("""
 Ca                0.00000000    0.00000000    0.00000000
 H                 2.72539390    0.00000000    0.00000000
""")

mol.spin = 1
mol.basis = '6-31g'
mol.build()
Norb = mol.nao_nr()

#UHF convergent criteria
e_conv = 1.e-8
d_conv = 1.e-8
if (mol.nelectron %2 != 0):
	nel_a = int(0.5 * (mol.nelectron + 1))
	nel_b = int(0.5 * (mol.nelectron - 1))
else:
	nel_a = int(0.5 * mol.nelectron)
	nel_b = nel_a

damp_value = 0.20



mf = scf.UHF(mol)
mf.init_guess = '1e'
mf.max_cycle = 200
mf.diis_start_cycle=200
mf.diis = False
mf.damp = damp_value #using damping instead of DIIS
mf.kernel()


pyscf_energy = mf.e_tot
print("Energy matches PySCF %s" % np.allclose(pyscf_energy, E_total))