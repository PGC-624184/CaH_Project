import numpy as np
from pyscf import gto, scf, mcscf
from pyscf import gto, mp, mcscf, dft, tddft
from pyscf.mcscf import dmet_cas
from pyscf import gto, scf, cc
import matplotlib.pyplot as plt

# Bond seperation lengths in angstrom
bond = np.arange(1.0, 10.0, 0.05)
dm_init = None

# results list, includes ground state and excited states
e_hf = [[],[],[]]

# Main iteration loop
for j in range(0,3): # iterate over excited states
    for r in reversed(bond): # iterate over bond lengths. Start out and work in for stability
        # Define the CaH molecule with the ccpvdz basis set
        mol = gto.M(atom=[['Ca', 0.0, 0.0, 0.0],
                          #['H', r, 0, -0.3865],
                          ['H' ,r, 0.0 , 0.0]],
                      basis='ccpvdz',
                      spin = 1
                    )
        # Reduced Hartree Fock solution as an initial case           
        mf = scf.UHF(mol).run()
        # define the excited state
        state_id = j
        # For the excited states, need to change the solver
        if state_id > 0:
            ao_labels = ['Ca 4pz', 'Ca 4py', 'Ca 4s',' Ca 4px', 'Ca 5s', 'Ca 5pz', 'Ca 5py',' Ca 5px']
            ncas, nelecas, mo = dmet_cas.guess_cas(mf, mf.make_rdm1(), ao_labels)
            mycas = mcscf.CASSCF(mf, ncas, nelecas).state_specific_(state_id)
            mycas.verbose = 2
            mycas.frozen = 6 # Freeze the inner 6 orbitals to save on computational effort (1s 2s 2p[xyz] 3s 3p[xyz] 4s)
            res = mycas.kernel(mo)
            #mycas.analyze()
            e_hf[j].append(res[0]) # append results
        else:
            mycc = cc.CCSD(mf).run()
            et = mycc.ccsd_t()
            e_hf[j].append(mycc.e_tot+et)


radius = bond # distance to atomic units from angstrom
# Convert the below to a difference from the far field ground state
# then convert result to eV
E_h = 27.211386246012257

ground_state = (e_hf[0][::-1]-e_hf[0][1])*E_h
first_state = (e_hf[1][::-1]-e_hf[0][1])*E_h
second_state = (e_hf[2][::-1]-e_hf[0][1])*E_h
third_state = 2*(e_hf[3][::-1]-e_hf[0][1])*E_h

# plots of the system
from matplotlib import pyplot as plt
plt.clf()
plt.plot(radius, ground_state[0:-3],label="Ground State")
plt.plot(radius, first_state,label="1st State")
plt.plot(radius, second_state,label="2nd State")
plt.plot(radius, third_state,label="3rd State")
plt.legend()
plt.xlabel("Bond Distance (angstrom)")
plt.ylabel("Bond Energy (eV)")
plt.title("Potential Energy curves for Ca I and H_2 excited states")
plt.show()
plt.savefig("CaH_2potentials.png")


""" 
More working - finding the wavelength at each of the jump energy values.
This should be the normal transition level out at the far field, and then
get funkier the closer you are.

Rest of the code is WIP and so results may not be real
"""

# distance between lines
dist_1 = np.array(e_hf[1]) - np.array(e_hf[0])
dist_2 = np.array(e_hf[2]) - np.array(e_hf[0])
dist_3 = np.array(e_hf[3]) - np.array(e_hf[0])

# This is h*c/E_hartrees when the result is converted to nanometers
factor = 45.56335252907954

# wavelength in nm
lambda_1 = factor/(2*dist_1)
lambda_2 = factor/(2*dist_2)
#lambda_3 = factor/(2*dist_3)


# plot results
plt.clf()
plt.plot(bond[0:-4], lambda_1[0:-4][::-1],label="1st to Ground")
plt.plot(bond, lambda_2[::-1],label="2nd to Ground")
plt.plot(bond, lambda_3[::-1],label="3rd to Ground")
plt.legend()
plt.show()
plt.savefig("Lambda_spread.png")


r = 27.83897594220896
r = 10
from functools import reduce
import numpy
from pyscf import gto, scf, mcscf, fci
mol = gto.M(atom=[  ['Ca', 0,     0,     0],
                    ['H',10,     0,    -0.3865],
                    ['H',10,     0,    0.3865]],
                    basis='def2-tzvp',
                    spin = 0,
                    verbose = 4)
        # Reduced Hartree Fock solution as an initial case           
mf = scf.RHF(mol)
mf.kernel()
#mf.analyze()

#
# 1. State-average CASSCF to get optimal orbitals
#
mc = mcscf.CASSCF(mf, 6, 4)
solver_ag = fci.direct_spin0_symm.FCI(mol)
solver_b2u = fci.direct_spin0_symm.FCI(mol)
solver_b2u.wfnsym = 'B2u'
mc = mcscf.state_average_mix(mc, [solver_ag,solver_b2u], [.5,.5])
cas_list = [17,20,21,22,23,30]  # 2pz orbitals
mo = mcscf.sort_mo(mc, mf.mo_coeff, cas_list)
mc.kernel(mo)
#mc.analyze()
mc_mo = mc.mo_coeff

#
# 2. Ground state wavefunction.  This step can be passed you approximate it
# with the state-averaged CASSCF wavefunction
#
mc = mcscf.CASCI(mf, 6, 6)
mc.fcisolver.wfnsym = 'Ag'
mc.kernel(mc_mo)
ground_state = mc.ci

#
# 3. Exited states.  In this example, B2u are bright states.
#
# Here, mc.ci[0] is the first excited state.
#
mc = mcscf.CASCI(mf, 6, 6)
mc.fcisolver.wfnsym = 'B2u'
mc.fcisolver.nroots = 8
mc.kernel(mc_mo)

#
# 4. transition density matrix and transition dipole
#
# Be careful with the gauge origin of the dipole integrals
#
charges = mol.atom_charges()
coords = mol.atom_coords()
nuc_charge_center = numpy.einsum('z,zx->x', charges, coords) / charges.sum()
mol.set_common_orig_(nuc_charge_center)
dip_ints = mol.intor('cint1e_r_sph', comp=3)

def makedip(ci_id):
    # transform density matrix in MO representation
    t_dm1 = mc.fcisolver.trans_rdm1(ground_state, mc.ci[ci_id], mc.ncas, mc.nelecas)
    # transform density matrix to AO representation
    orbcas = mc_mo[:,mc.ncore:mc.ncore+mc.ncas]
    t_dm1_ao = reduce(numpy.dot, (orbcas, t_dm1, orbcas.T))
    # transition dipoles
    return numpy.einsum('xij,ji->x', dip_ints, t_dm1_ao)

# 1st and 6th excited states are B2u of D6h point group, dark states
# 3rd and 4th excited states are triplet states, dipole == 0
for i in range(8):
    print('Transition dipole between |0> and |%d>'%(i+1), makedip(i))


from pyscf.tools import cubegen
# Save Orbitals
for j in range(0,48):
    cubegen.orbital(mol, 'CaH2_mo_'+str(j)+'.cub', myhf.mo_coeff[:,j])
