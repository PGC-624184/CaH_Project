from pyscf import gto, scf, mcscf
from pyscf.mcscf import avas
from pyscf import dmrgscf
from pyscf import mrpt
from pyscf import lib

#
# Adjust the MPI schedular and scratch directory if needed.
#
#from pyscf.dmrgscf import settings
#settings.MPIPREFIX = 'srun'
#settings.BLOCKSCRATCHDIR = '/scratch'


mol = gto.M(
    atom = '''
Ca      0.000000    0.000000    0.000000
H       8.000000    0.000000    0.386500
H       8.000000    0.000000   -0.386500
''',
    basis = 'cc-pvtz-dk',
    spin = 0,
    charge = 0,
    verbose = 4)

mf = scf.sfx2c1e(scf.ROHF(mol))
mf.chkfile = 'CaH2.chk'
mf.kernel()

#
# Generate the DMRG-CASSCF initial guess with AVAS (atomic valence active space) method.
#
norb_act, ne_act, mos = avas.avas(mf, ['Ca 4s', 'Ca 4p', 'H 1s', 'H 2s'])

#
# Pass 1
# DMRG-CASSCF calculation with small M.  In many systems, the mcscf orbitals has
# weak dependence to the quality of the active space solution.  We can use small
# M to approximate the active space wave function and increase M for DMRG-CASCI
# and DMRG-NEVPT2 in next step.
#
mc = dmrgscf.DMRGSCF(mf, norb_act, ne_act)
mc.fcisolver.maxM = 500
mc.fcisolver.extraline = ['num_thrds %d'%lib.num_threads(), 'warmup local_2site']
mc.fcisolver.memory = 100  # GB
mc.conv_tol = 1e-6
mc.state_average_([.2]*3)
mc.kernel(mos)
mc_orbs = mc.mo_coeff


#
# Pass 2
# A separated DMRG-CASCI calculation based on DMRG-CASSCF orbitals obtained and
# DMRG solver with large bond dimension to improve the multi configurational
# wave function.
#
mc = mcscf.CASCI(mf, norb_act, ne_act)
mc.fcisolver = dmrgscf.DMRGCI(mol)
mc.fcisolver.maxM = 1500
mc.fcisolver.extraline = ['num_thrds %d'%lib.num_threads(), 'warmup local_2site']
mc.fcisolver.memory = 100  # GB
mc.fcisolver.nroots = 3
mc.kernel(mc_orbs)


#
# This step computes DMRG-NEVPT2 energy state-specificly. It is recommended to
# use large bond dimension (larger than the last DMRG-CASCI calculation) to
# guarantee the accuracy.
#
print mrpt.NEVPT(mc, root=0).compress_approx(maxM=2000).kernel()
print mrpt.NEVPT(mc, root=1).compress_approx(maxM=2000).kernel()
print mrpt.NEVPT(mc, root=2).compress_approx(maxM=2000).kernel()