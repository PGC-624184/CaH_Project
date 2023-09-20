import numpy as np
from pyscf import gto, scf, mcscf
from pyscf import gto, mp, mcscf, dft, tddft
from pyscf.mcscf import dmet_cas
from pyscf import gto, scf, cc

# Bond seperation lengths in angstrom
bond = np.arange(1.0, 10.0, 0.05)
dm_init = None

# results list, includes ground state and excited states
e_hf = [[],[]]

# Main iteration loop
for j in range(0,3): # iterate over excited states
    for r in reversed(bond): # iterate over bond lengths. Start out and work in for stability
        # Define the CaH molecule with the ccpvdz basis set
        mol = gto.M(atom=[['Ca', 0, 0, 0],
                          ['H', r, 0, -0.3865],
                          ['H' ,r, 0 , 0.3865]],
                      basis='ccpvdz',
                      spin = 0
                    )
        # Reduced Hartree Fock solution as an initial case           
        mf = scf.RHF(mol).run()
        # define the excited state
        state_id = j
        # For the excited states, need to change the solver
        if state_id > 0:
            ao_labels = ['Ca 4pz', 'Ca 4py', 'Ca 4s',' Ca 4px']
            ncas, nelecas, mo = dmet_cas.guess_cas(mf, mf.make_rdm1(), ao_labels)
            mycas = mcscf.CASSCF(mf, ncas, nelecas).state_specific_(state_id)
            mycas.verbose = 2
            mycas.frozen = 6 # Freeze the inner 6 orbitals to save on computational effort
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
plt.plot(radius, ground_state,label="Ground State")
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