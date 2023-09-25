import pyscf    as pyscf
import numpy    as np
import sys      as sys
from pyscf      import gto, scf, ao2mo, fci, ci, cc, mcscf
from pyscf import tdscf
from pyscf.mcscf import dmet_cas
from pyscf.lib  import logger
import matplotlib.pyplot as plt

def calculate_energy(r) -> float:
    """
    Calculate the energy of a molecular system at a specific bond radius.

    This function uses quantum chemistry methods to calculate the energy of a molecular
    system at a given bond radius. It first determines the energy using an Unrestricted
    Hartree Fock (UHF) method and then continues to refine the calculation with the CCSD(T) (Coupled Cluster with Singles and Doubles, with perturbative Triples) computation.

    Args:
        r (float): The bond radius at which to calculate the energy.

    Returns:
        float: The total energy of the molecular system, including contributions from both
               CCSD(T) and Hartree Fock calculations.

    Raises:
        SystemExit: If the optimization process for the wavefunction does not converge within
                    200 optimization cycles.

    Example:
        >>> energy = calculate_energy(1.0)
        >>> print(energy)
        -75.123456789  # Replace with the actual energy value obtained

    Note:
        This function relies on external modules like 'pyscf' and 'logger,' which should
        be properly imported and available in the environment. Make sure to set up the
        appropriate computational environment before using this function.
    """
    neutralMolecule     = pyscf.M(atom = [
                                    ['Ca', (0,0,0)], 
                                    ['H', (r,0,0)]],
                                    basis = 'sto3g', 
                                    verbose = 3, 
                                    spin = 1, 
                                    symmetry = True)
    
    neutralMoleculeHF   = scf.UHF(neutralMolecule).run()

    logObject               =   logger.new_logger(neutralMoleculeHF)                # Creates a new logger object
    orbitals,_,stability,_  =   neutralMoleculeHF.stability(return_status = True)   # Returns four variables: Internal MO, External MO, Internal Stability, External Stability
    optimizationCycles      =   0
    while (not stability and optimizationCycles < 200): 
        logObject.note("Currently at {}-th step optimizing the wavefunction".format(optimizationCycles))
    densityMatrix           =   neutralMoleculeHF.make_rdm1(orbitals, neutralMoleculeHF.mo_occ)
    neutralMoleculeHF       =   neutralMoleculeHF.run(densityMatrix) 
    orbitals,_,stability,_  =   neutralMoleculeHF.stability(return_status = True)
    optimizationCycles      =   optimizationCycles + 1
    if not stability:
        print("Unable to achieve stable wavefunction after {} steps. Exiting".format(optimizationCycles))
        sys.exit(-1)
    
    # Now start Post-HF calculations
    neutralMoleculeCC           =   cc.UCCSD(neutralMoleculeHF).set(conv_tol = 1e-7, direct = True)
    neutralMoleculeCC.run()
    neutralMoleculeCC_results   =   neutralMoleculeCC.kernel()
    neutralMoleculeCC_CCSDT     =   neutralMoleculeCC.ccsd_t()   
    return neutralMoleculeCC.e_tot + neutralMoleculeCC_CCSDT

def calculate_excited(r: float,state_id: int) -> float:
    """
    Calculate the energy of an excited state of a molecular system at a specific bond radius.

    This function computes the energy of an excited state for a molecular system with
    a given bond length 'r'. It utilizes the PySCF library to perform electronic structure
    calculations, including Hartree-Fock (HF) and Complete Active Space Configuration
    Interaction (CASCI) methods to obtain the excited state energy.

    Args:
        r (float): The bond radius (in Angstrom) at which to calculate the energy.
        state_id (int): The state-specific identifier for the excited state to compute.

    Returns:
        float: The energy of the specified excited state (in Hartrees).

    Example:
        >>> bond_length = 1.2
        >>> state_number = 2
        >>> energy = calculate_excited(bond_length, state_number)

    Note:
        - This function relies on the 'pyscf' library for quantum chemistry calculations.
        - Ensure that 'pyscf' is properly installed in your environment before using this function.
        - The 'state_id' parameter determines which excited state is calculated within the
          CASCI method, and it should be an integer corresponding to the desired state, noting 0
          is the Ground State.
    """

    excitedMolecule     = pyscf.M(atom = [
                                    ['Ca', (0,0,0)], 
                                    ['H', (r,0,0)]], 
                                    basis = 'sto3g',  
                                    verbose = 3, 
                                    spin = 1, 
                                    symmetry = True)

    excitedMoleculeHF   = scf.UHF(excitedMolecule).run()
    ao_labels = ['Ca 4pz', 'Ca 4py', 'Ca 4s',' Ca 4px', 'Ca 5s']
    ncas, nelecas, mo = dmet_cas.guess_cas(excitedMoleculeHF, excitedMoleculeHF.make_rdm1(), ao_labels)
    mc = mcscf.CASSCF(excitedMoleculeHF, ncas, nelecas).state_specific_(state_id)
    mc.kernel(mo)
    mo = mc.mo_coeff

    mc = mcscf.CASCI(excitedMoleculeHF, ncas, nelecas).state_specific_(state_id)
    emc = mc.casci(mo)[0]
    return emc


def plot_curve(bond,energy,label):
    """
    Plot the energy curve for a molecular system as a function of bond length.

    This function takes two lists, 'bond' and 'energy,' and creates a plot of the energy
    as a function of bond length for a molecular system. It is typically used to visualize
    the potential energy surface of a molecule with varying bond lengths.

    Args:
        bond (list): A list of bond lengths in Angstrom.
        energy (list): A list of corresponding energy values in Hartrees.

    Returns:
        None

    Example:
        >>> bond_lengths = [1.0, 1.1, 1.2, 1.3]
        >>> energy_values = [-75.0, -74.5, -74.2, -74.0]
        >>> plot_curve(bond_lengths, energy_values)

    Note:
        This function uses the 'matplotlib' library for plotting. Ensure that 'matplotlib'
        is properly installed in your environment before using this function.
    """
    plt.plot(bond, energy[::-1],label=label)
    plt.xlabel("Bond length (Angstrom)")
    plt.ylabel("Energy (eV)")
    plt.legend()


def run_ground_sim(start, stop, step):
    """
    Run a simulation to calculate energy as a function of bond length.

    This function performs a simulation to calculate the energy of a molecular system
    at various bond lengths within a specified range. It iterates over the bond lengths,
    calculates the energy using the 'calculate_energy' function, and stores the results.

    Args:
        start (float): The starting bond length (in Angstrom) for the simulation.
        stop (float): The ending bond length (in Angstrom) for the simulation.
        step (float): The step size (in Angstrom) between bond lengths.

    Returns:
        tuple: A tuple containing two lists:
            - bond (list): A list of bond lengths in Angstrom.
            - e_vec (list): A list of energy differences (in Hartrees) relative to the
              energy at the second bond length.

    Example:
        >>> start_length = 1.0
        >>> end_length = 2.0
        >>> step_size = 0.1
        >>> bond_lengths, energy_differences = run_sim(start_length, end_length, step_size)

    Note:
        - This function relies on the 'calculate_energy' function for energy calculations.
        - Ensure that the 'calculate_energy' function is properly defined and imported
          in your environment before using this function.
    """
    bond = np.arange(start,stop,step)
    ground_vec = []
    for r in reversed(bond):
        ground_energy = calculate_energy(r)
        ground_vec.append(ground_energy)

    return np.array(bond), np.array(ground_vec)

def run_excited_sim(start, stop, step, state):
    """
    Run a simulation to calculate energy as a function of bond length.

    This function performs a simulation to calculate the energy of a molecular system
    at various bond lengths within a specified range. It iterates over the bond lengths,
    calculates the energy using the 'calculate_energy' function, and stores the results.

    Args:
        start (float): The starting bond length (in Angstrom) for the simulation.
        stop (float): The ending bond length (in Angstrom) for the simulation.
        step (float): The step size (in Angstrom) between bond lengths.

    Returns:
        tuple: A tuple containing two lists:
            - bond (list): A list of bond lengths in Angstrom.
            - e_vec (list): A list of energy differences (in Hartrees) relative to the
              energy at the second bond length.

    Example:
        >>> start_length = 1.0
        >>> end_length = 2.0
        >>> step_size = 0.1
        >>> bond_lengths, energy_differences = run_sim(start_length, end_length, step_size)

    Note:
        - This function relies on the 'calculate_energy' function for energy calculations.
        - Ensure that the 'calculate_energy' function is properly defined and imported
          in your environment before using this function.
    """
    bond = np.arange(start,stop,step)
    excited_vec = []
    for r in reversed(bond):
        excited_energy = calculate_excited(r, state)
        excited_vec.append(excited_energy)
    return np.array(excited_vec)

if __name__=="__main__":
    # Interatomic distance range
    start = 1.75
    stop = 10.0
    step = 0.5

    # Calculate curves
    bond, ground_results = run_ground_sim(start,stop,step)
    excited_results0 = run_excited_sim(start,stop,step,0)

    
    # Hartrees Energy (a.u.) to eV
    E_h = 27.211386246012257
    # Set the far field as 0 for all curves
    zero_point = ground_results[1]
    ground_results = (ground_results - zero_point)*E_h
    excited_results0 = (excited_results0 - zero_point)*E_h
   

    # plot the results
    plt.clf()
    plot_curve(bond, ground_results,"CCSD(T) Ground State")
    plot_curve(bond, excited_results0,"MCSCF Ground State")
    plt.legend(loc='best')
    plt.show()

    difference0 = np.array(excited_results0) - np.array(ground_results)
   
    # This is the conversion factor to calculate wavelength, which is c*h/energy
    factor = 45.56335252907954    

    lambda0 = factor/difference0

    plt.clf()
    plt.plot(bond,lambda0[::-1],label = "MCSCF Ground")
    plt.xlabel("Bond Distance (Angstrom)")
    plt.ylabel("Resonant Wavelength (nm)")
    plt.legend(loc="best")
    plt.show()