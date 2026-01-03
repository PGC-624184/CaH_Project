import numpy as np
import pandas as pd
from pyscf import gto, dft, tddft
from tqdm import tqdm
import os


def rotation_matrix(theta, phi):
    """
    Create rotation matrix for H2 molecule orientation.
    
    theta: polar angle (0 to pi)
    phi: azimuthal angle (0 to 2*pi)
    
    Returns 3x3 rotation matrix
    """
    # Rotation around y-axis by theta
    Ry = np.array([
        [np.cos(theta), 0, np.sin(theta)],
        [0, 1, 0],
        [-np.sin(theta), 0, np.cos(theta)]
    ])
    
    # Rotation around z-axis by phi
    Rz = np.array([
        [np.cos(phi), -np.sin(phi), 0],
        [np.sin(phi), np.cos(phi), 0],
        [0, 0, 1]
    ])
    
    return Rz @ Ry


def generate_h2_orientations(r_ca_h2, n_theta=5, n_phi=8, r_h_h=0.74):
    """
    Generate H2 molecular orientations around Ca atom using spherical sampling.
    
    Args:
        r_ca_h2: Distance from Ca to H2 center of mass (angstrom)
        n_theta: Number of polar angle samples (5-10 recommended)
        n_phi: Number of azimuthal angle samples (8-16 recommended)
        r_h_h: H-H bond length (0.74 Å is equilibrium)
    
    Returns:
        List of molecule configurations: [[[Ca], [H1], [H2]], weight]
        Weights are for proper spherical averaging (sin(theta) factor)
    """
    configurations = []
    
    # Use Gauss-Legendre quadrature for theta (more efficient than uniform)
    theta_points, theta_weights = np.polynomial.legendre.leggauss(n_theta)
    # Map from [-1, 1] to [0, pi]
    thetas = np.arccos(theta_points)
    
    # Uniform sampling in phi
    phis = np.linspace(0, 2*np.pi, n_phi, endpoint=False)
    phi_weight = 2 * np.pi / n_phi
    
    # H2 molecule initially along x-axis, centered at origin
    h2_vector = np.array([r_h_h / 2, 0, 0])
    
    for i, (theta, theta_w) in enumerate(zip(thetas, theta_weights)):
        for phi in phis:
            # Rotate H2 molecule
            R = rotation_matrix(theta, phi)
            h1_local = R @ h2_vector
            h2_local = R @ (-h2_vector)
            
            # Translate to distance r_ca_h2 from Ca (at origin)
            # Place H2 center of mass along z-axis at distance r
            com_position = np.array([0, 0, r_ca_h2])
            h1_pos = com_position + h1_local
            h2_pos = com_position + h2_local
            
            # Create molecule configuration
            molecule = [
                ["Ca", 0, 0, 0],
                ["H", float(h1_pos[0]), float(h1_pos[1]), float(h1_pos[2])],
                ["H", float(h2_pos[0]), float(h2_pos[1]), float(h2_pos[2])]
            ]
            
            # Weight includes theta_weight (from quadrature) and phi weight
            weight = theta_w * phi_weight / (4 * np.pi)
            
            configurations.append((molecule, weight))
    
    return configurations


def run_computation(
    r=15,
    states=10,
    functional="b3lyp",
    basis="sto3g",
    atom=None,
    spin=0,
):
    """
    Perform electronic structure calculations for a molecular configuration.
    
    Args:
        r: Bond separation length (for reference, not used directly)
        states: Number of excited states to compute
        functional: DFT functional
        basis: Basis set
        atom: Atomic coordinates [[element, x, y, z], ...]
        spin: Spin multiplicity
    
    Returns:
        tuple: (energies, oscillator_strengths) or (None, None) if failed
    """
    try:
        mol = gto.M(
            atom=atom,
            basis=basis,
            spin=spin,
            verbose=0,  # Reduce verbosity for many calculations
        )
        
        if mol.spin == 1:
            mf = dft.UKS(mol)
        else:
            mf = dft.RKS(mol)
        
        mf.xc = functional
        mf.kernel()
        
        if not mf.converged:
            print(f"Warning: SCF did not converge for r={r}")
            return None, None
        
        mytd = tddft.TDDFT(mf)
        mytd.singlet = True
        mytd.nstates = states
        results = mytd.kernel()
        
        oscillator = mytd.oscillator_strength(gauge="length")
        
        return results[0], oscillator
        
    except Exception as e:
        print(f"Calculation failed for r={r}: {e}")
        return None, None


def compute_orientation_averaged(
    r_val,
    states=20,
    functional="wB97X_V",
    basis="def2-QZVPP",
    spin=1,
    n_theta=5,
    n_phi=8,
):
    """
    Compute orientation-averaged energies and oscillator strengths.
    
    Args:
        r_val: Ca to H2 center-of-mass distance
        states: Number of excited states
        functional: DFT functional
        basis: Basis set
        spin: Spin multiplicity
        n_theta: Number of polar angle samples (5-10 for good balance)
        n_phi: Number of azimuthal angle samples (8-16 recommended)
    
    Returns:
        tuple: (averaged_energies, averaged_oscillators)
    """
    # Generate orientations
    configurations = generate_h2_orientations(r_val, n_theta, n_phi)
    
    print(f"\n  Averaging over {len(configurations)} orientations...")
    
    # Storage for results
    all_energies = []
    all_oscillators = []
    all_weights = []
    
    # Compute for each orientation
    for molecule, weight in tqdm(configurations, 
                                  desc=f"  r={r_val:.2f}Å",
                                  leave=False):
        energies, osc = run_computation(
            r_val,
            states=states,
            functional=functional,
            basis=basis,
            atom=molecule,
            spin=spin,
        )
        
        if energies is not None:
            all_energies.append(energies)
            all_oscillators.append(osc)
            all_weights.append(weight)
    
    if not all_energies:
        print(f"  All calculations failed for r={r_val}")
        return None, None
    
    # Convert to arrays for weighted averaging
    all_energies = np.array(all_energies)  # Shape: (n_orientations, n_states)
    all_oscillators = np.array(all_oscillators)
    all_weights = np.array(all_weights)
    
    # Normalize weights (should already sum to 1, but just in case)
    all_weights = all_weights / all_weights.sum()
    
    # Weighted average
    avg_energies = np.average(all_energies, axis=0, weights=all_weights)
    avg_oscillators = np.average(all_oscillators, axis=0, weights=all_weights)
    
    return avg_energies, avg_oscillators


def main():
    """Main execution with orientation averaging."""
    
    # Configuration
    radius = np.arange(1.5, 4.25, 0.25)
    comp_states = 20
    base = "def2-QZVPP"
    func = "wB97X_V"
    spin_m = 1
    
    # Orientation sampling (adjust for accuracy vs speed)
    # n_theta=5, n_phi=8 gives 40 orientations (good starting point)
    # n_theta=7, n_phi=12 gives 84 orientations (better accuracy)
    n_theta = 5  
    n_phi = 8
    
    print(f"Configuration: {len(radius)} radii × {n_theta*n_phi} orientations")
    print(f"Total calculations: {len(radius) * n_theta * n_phi}")
    
    # Create output directory
    os.makedirs("data/CaH2_oriented", exist_ok=True)
    
    # Run computations
    for r_val in tqdm(radius, desc="Computing Ca-H₂ curves"):
        avg_energies, avg_osc = compute_orientation_averaged(
            r_val,
            states=comp_states,
            functional=func,
            basis=base,
            spin=spin_m,
            n_theta=n_theta,
            n_phi=n_phi,
        )
        
        if avg_energies is not None:
            filename = f"data/CaH2_oriented/Averaged_curve_data_Ca_H2_r{r_val:.2f}.csv"
            df_E = pd.DataFrame({
                "Energies": avg_energies,
                "f": avg_osc
            })
            df_E.to_csv(filename, sep="\t")
            print(f"✓ Saved: {filename}")


if __name__ == "__main__":
    main()
