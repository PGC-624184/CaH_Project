import numpy
from pyscf import gto, scf, mcscf

mol = gto.M(
    atom = 'Ca 0 0 0; H 0 0 1.2',
    basis = 'ccpvdz',
    spin = 1)
myhf = scf.RHF(mol).run()

# 6 orbitals, 8 electrons
mycas = mcscf.CASSCF(myhf, 6, 3)

#
# Freeze the two innermost oxygen 1s orbitals in the orbital
# optimization
#
mycas.frozen = 16
mycas.kernel()

#
# Freeze orbitals based on the list of indices.  Two HF core orbitals and two HF
# virtual orbitals are excluded from CASSCF optimization.
#
mycas.frozen = [0,1,26,27]
mycas.kernel()

#
# Partially freeze the active space so that the frozen orbitals are always in
# the active space.  It can help CASSCF converge to reasonable solution.
#
mycas.frozen = [5,6,7,8]
mycas.kernel()