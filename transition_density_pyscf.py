import numpy as np
from pyscf import gto, scf, mcscf
# Bond seperation lengths in angstrom
bond = np.arange(1.0, 10.0, 1)
dm_init = None

# results list, includes ground state and excited states
e_hf = [[],[],[],[],[]]

# Main iteration loop
for j in range(0,5): # iterate over excited states
    for r in reversed(bond): # iterate over bond lengths. Start out and work in for stability
        # Define the CaH molecule with the ccpvdz basis set
        mol = gto.M(atom=[['Ca', 0, 0, 0],
                          ['H', r, 0, 0]],
                      basis='ccpvdz',
                      spin = 1
                    )
        # Reduced Hartree Fock solution as an initial case           
        mf = scf.RHF(mol).run()
        # define the excited state
        state_id = j
        # Run the Complete Active Space SCF. 4 spaces with 3 electrons
        mc = mcscf.CASSCF(mf, 4, 3).state_specific_(state_id)
        mc.verbose = 4
        # For the excited states, need to change the solver
        if state_id > 0:
            mc.kernel() # excited state solve
            mo = mc.mo_coeff # molecular orbitals
            mc.fcisolver.nroots = 4 # solve for 4 roots
            emc = mc.casci(mo) # solver for the energy
            e_hf[j].append(emc[0]) # append results
        else:
            emc = mc.kernel()[0] # Ground state solve
            e_hf[j].append(emc) # append results



# plots of the system
from matplotlib import pyplot as plt
plt.clf()
plt.plot(bond/0.529177210903, e_ground[::-1],label="Ground State") # converts angstrom to a.u.
plt.plot(bond/0.529177210903, e_hf_1[::-1],label="1st State")
plt.plot(bond/0.529177210903, e_hf_2[::-1],label="2nd State")
plt.plot(bond/0.529177210903, e_hf_3[::-1],label="3rd State")
plt.legend()
plt.show()



""" 
More working - finding the wavelength at each of the jump energy values.
This should be the normal transition level out at the far field, and then
get funkier the closer you are.

Rest of the code is WIP and so results may not be real
"""

# distance between lines
dist_1 = np.array(e_hf_1) - np.array(e_ground)
dist_2 = np.array(e_hf_2) - np.array(e_ground)
dist_3 = np.array(e_hf_3) - np.array(e_ground)

# This is h*c/E_hartrees when the result is converted to nanometers
factor = 45.56335252907954

# wavelength in nm
lambda_1 = factor/dist_1
lambda_2 = factor/dist_2
lambda_3 = factor/dist_3


# plot results
plt.clf()
plt.plot(bond, lambda_1[::-1],label="1st to Ground")
plt.plot(bond, lambda_2[::-1],label="2nd to Ground")
plt.plot(bond, lambda_3[::-1],label="3rd to Ground")
plt.legend()
plt.show()




