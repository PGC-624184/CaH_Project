import numpy as np
from pyscf import gto, scf, mcscf
# Bond seperation lengths in angstrom
bond = np.arange(1.0, 10.0, 0.25)
dm_init = None

# results list, includes ground state and excited states
e_hf = [[],[],[],[]]

# Main iteration loop
for j in range(0,5): # iterate over excited states
    for r in reversed(bond): # iterate over bond lengths. Start out and work in for stability
        # Define the CaH molecule with the ccpvdz basis set
        mol = gto.M(atom=[['Ca', 0, 0, 0],
                          ['H', r, 0, 0],
                          ['H', r, 0, 0.737],],
                      basis='ccpvdz',
                      spin = 0
                    )
        # Reduced Hartree Fock solution as an initial case           
        mf = scf.RHF(mol).run()
        # define the excited state
        state_id = j
        # Run the Complete Active Space SCF. 6 spaces with 4 electrons
        mc = mcscf.CASSCF(mf, 6, 4).state_specific_(state_id)
        mc.verbose = 4
        # For the excited states, need to change the solver
        if state_id > 0:
            mc.kernel() # excited state solve
            mo = mc.mo_coeff # molecular orbitals
            mc.fcisolver.nroots = 6 # solve for 4 roots
            emc = mc.casci(mo) # solver for the energy
            e_hf[j].append(emc[0]) # append results
        else:
            emc = mc.kernel()[0] # Ground state solve
            e_hf[j].append(emc) # append results



radius = bond/0.529177210903 # distance to atomic units from angstrom
# Convert the below to a difference from the far field ground state
# then convert result to eV
ground_state = (e_hf[0][::-1]-e_hf[0][1])*27.211386246012257
first_state = (e_hf[1][::-1]-e_hf[0][1])*27.211386246012257
second_state = (e_hf[2][::-1]-e_hf[0][1])*27.211386246012257
third_state = (e_hf[3][::-1]-e_hf[0][1])*27.211386246012257

# plots of the system
from matplotlib import pyplot as plt
plt.clf()
plt.plot(radius, ground_state,label="Ground State")
plt.plot(radius, first_state,label="1st State")
plt.plot(radius, second_state,label="2nd State")
plt.plot(radius, third_state,label="3rd State")
plt.legend()
plt.xlabel("Bond Distance (a.u.)")
plt.ylabel("Bond Energy (eV)")
plt.title("Potential Energy curves for CaH excited states")
plt.show()



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
lambda_1 = dist_1
lambda_2 = dist_2
lambda_3 = dist_3


# plot results
plt.clf()
plt.plot(bond, lambda_1[::-1],label="1st to Ground")
plt.plot(bond, lambda_2[::-1],label="2nd to Ground")
plt.plot(bond, lambda_3[::-1],label="3rd to Ground")
plt.legend()
plt.show()




