import numpy as np
import pandas as pd
from pyscf import gto, dft, tddft
from tqdm import tqdm
import os
from scipy.constants import physical_constants, c, epsilon_0, pi

# Physical constants
BOHR_TO_ANGSTROM = physical_constants['Bohr radius'][0] * 1e10
HARTREE_TO_EV = physical_constants['Hartree energy in eV'][0]
HARTREE_TO_JOULE = physical_constants['Hartree energy'][0]


def rotation_matrix(theta, phi):
    """
    Create rotation matrix for H2 molecule orientation.
    
    theta: polar angle (0 to pi/2 due to symmetry)
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
    Generate H2 molecular orientations using symmetry to reduce calculations.
    
    H2 has inversion symmetry: (H1-H2) is identical to (H2-H1)
    So we only sample theta from [0, pi/2] instead of [0, pi]
    This reduces calculations by factor of 2!
    
    Args:
        r_ca_h2: Distance from Ca to H2 center of mass (angstrom)
        n_theta: Number of polar angle samples (now only 0 to pi/2)
        n_phi: Number of azimuthal angle samples
        r_h_h: H-H bond length (0.74 Å is equilibrium)
    
    Returns:
        List of molecule configurations: [[[Ca], [H1], [H2]], weight]
    """
    configurations = []
    
    # Use Gauss-Legendre quadrature for theta in [0, pi/2] only
    # Map from [-1, 1] to [0, pi/2] using symmetry
    theta_points, theta_weights = np.polynomial.legendre.leggauss(n_theta)
    thetas = (theta_points + 1) * np.pi / 4  # Maps [-1,1] to [0, pi/2]
    
    # Uniform sampling in phi
    phis = np.linspace(0, 2*np.pi, n_phi, endpoint=False)
    phi_weight = 2 * np.pi / n_phi
    
    # H2 molecule initially along x-axis, centered at origin
    h2_vector = np.array([r_h_h / 2, 0, 0])
    
    for theta, theta_w in zip(thetas, theta_weights):
        for phi in phis:
            # Rotate H2 molecule
            R = rotation_matrix(theta, phi)
            h1_local = R @ h2_vector
            h2_local = R @ (-h2_vector)
            
            # Translate to distance r_ca_h2 from Ca (at origin)
            com_position = np.array([0, 0, r_ca_h2])
            h1_pos = com_position + h1_local
            h2_pos = com_position + h2_local
            
            # Create molecule configuration
            molecule = [
                ["Ca", 0, 0, 0],
                ["H", float(h1_pos[0]), float(h1_pos[1]), float(h1_pos[2])],
                ["H", float(h2_pos[0]), float(h2_pos[1]), float(h2_pos[2])]
            ]
            
            # Weight: theta_weight * phi_weight / (2*pi) 
            # Note: integral over hemisphere (0 to pi/2) instead of full sphere
            weight = theta_w * phi_weight / (2 * np.pi)
            
            configurations.append((molecule, weight))
    
    return configurations


def calculate_cross_section(energies, oscillator_strengths, gamma_natural=1e8):
    """
    Calculate absorption cross section for each transition.
    
    Uses the quantum mechanical formula for absorption cross section:
    σ(ν) = (π e² / m_e c ε₀) × f × φ(ν)
    
    Where φ(ν) is the line profile (Lorentzian for natural broadening)
    
    Args:
        energies: Transition energies in Hartree
        oscillator_strengths: Dimensionless oscillator strengths
        gamma_natural: Natural linewidth in Hz (default ~1e8 Hz for typical transitions)
    
    Returns:
        dict with:
            - 'wavelengths': Wavelengths in nm
            - 'frequencies': Frequencies in Hz
            - 'cross_sections_peak': Peak cross sections in cm²
            - 'energies_eV': Energies in eV
    """
    # Constants
    e = physical_constants['elementary charge'][0]
    m_e = physical_constants['electron mass'][0]
    
    # Convert energies to frequencies
    frequencies = energies * HARTREE_TO_JOULE / physical_constants['Planck constant'][0]
    
    # Convert to wavelengths (nm)
    wavelengths = c / frequencies * 1e9
    
    # Classical absorption cross section at line center (Lorentzian peak)
    # σ₀ = (π e² / m_e c ε₀) × f / (π × Γ)
    # Simplifies to: σ₀ = (e² / 2 ε₀ m_e c) × (f / Γ)
    
    prefactor = (np.pi * e**2) / (m_e * c * epsilon_0)  # in SI units
    
    # Peak cross section (at resonance, assuming Lorentzian profile)
    # For Lorentzian: φ(ν₀) = 1/(π × Γ) where Γ is FWHM
    cross_sections_peak = (prefactor * oscillator_strengths / 
                          (np.pi * gamma_natural)) * 1e4  # Convert m² to cm²
    
    # Integrated cross section (more robust, doesn't depend on linewidth)
    # ∫σ(ν)dν = (π e² / m_e c ε₀) × f
    cross_sections_integrated = prefactor * oscillator_strengths * 1e4  # cm² Hz
    
    return {
        'wavelengths_nm': wavelengths,
        'frequencies_Hz': frequencies,
        'energies_eV': energies * HARTREE_TO_EV,
        'cross_sections_peak_cm2': cross_sections_peak,
        'cross_sections_integrated_cm2_Hz': cross_sections_integrated,
        'oscillator_strengths': oscillator_strengths,
    }


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
    Uses H2 symmetry to reduce calculations by factor of 2.
    
    Args:
        r_val: Ca to H2 center-of-mass distance
        states: Number of excited states
        functional: DFT functional
        basis: Basis set
        spin: Spin multiplicity
        n_theta: Number of polar angle samples (0 to π/2 only, thanks to symmetry!)
        n_phi: Number of azimuthal angle samples
    
    Returns:
        tuple: (averaged_energies, averaged_oscillators, std_energies, std_oscillators)
    """
    # Generate orientations (now only n_theta × n_phi instead of 2 × n_theta × n_phi)
    configurations = generate_h2_orientations(r_val, n_theta, n_phi)
    
    print(f"  Using H2 symmetry: {len(configurations)} orientations (2× reduction)")
    
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
        return None, None, None, None
    
    # Convert to arrays for weighted averaging
    all_energies = np.array(all_energies)
    all_oscillators = np.array(all_oscillators)
    all_weights = np.array(all_weights)
    
    # Normalize weights
    all_weights = all_weights / all_weights.sum()
    
    # Weighted average
    avg_energies = np.average(all_energies, axis=0, weights=all_weights)
    avg_oscillators = np.average(all_oscillators, axis=0, weights=all_weights)
    
    # Standard deviation (for uncertainty estimation)
    variance_energies = np.average((all_energies - avg_energies)**2, 
                                   axis=0, weights=all_weights)
    variance_oscillators = np.average((all_oscillators - avg_oscillators)**2,
                                      axis=0, weights=all_weights)
    
    std_energies = np.sqrt(variance_energies)
    std_oscillators = np.sqrt(variance_oscillators)
    
    return avg_energies, avg_oscillators, std_energies, std_oscillators


def main():
    """Main execution with orientation averaging and cross section calculation."""
    
    # Configuration
    radius = np.arange(1.5, 4.25, 0.25)
    comp_states = 20
    base = "def2-QZVPP"
    func = "wB97X_V"
    spin_m = 1
    
    # Orientation sampling (exploiting H2 symmetry for 2× speedup!)
    # n_theta=5, n_phi=8 → 40 orientations (was 80 without symmetry)
    # n_theta=7, n_phi=12 → 84 orientations (was 168 without symmetry)
    n_theta = 5  
    n_phi = 8
    
    # Natural linewidth for Ca transitions (adjust based on your specific transition)
    gamma_natural = 1e8  # Hz, typical for allowed transitions
    
    print(f"Configuration: {len(radius)} radii × {n_theta*n_phi} orientations")
    print(f"Total calculations: {len(radius) * n_theta * n_phi}")
    print(f"(Factor of 2 reduction thanks to H2 symmetry!)")
    
    # Create output directory
    os.makedirs("data/CaH2_oriented", exist_ok=True)
    
    # Run computations
    results_summary = []
    
    for r_val in tqdm(radius, desc="Computing Ca-H₂ curves"):
        avg_energies, avg_osc, std_energies, std_osc = compute_orientation_averaged(
            r_val,
            states=comp_states,
            functional=func,
            basis=base,
            spin=spin_m,
            n_theta=n_theta,
            n_phi=n_phi,
        )
        
        if avg_energies is not None:
            # Calculate cross sections
            cross_section_data = calculate_cross_section(
                avg_energies, 
                avg_osc, 
                gamma_natural
            )
            
            # Save detailed data
            filename = f"data/CaH2_oriented/Averaged_curve_data_Ca_H2_r{r_val:.2f}.csv"
            df_E = pd.DataFrame({
                "Energies_Hartree": avg_energies,
                "Energies_eV": cross_section_data['energies_eV'],
                "Wavelengths_nm": cross_section_data['wavelengths_nm'],
                "Oscillator_strength": avg_osc,
                "Oscillator_strength_std": std_osc,
                "Cross_section_peak_cm2": cross_section_data['cross_sections_peak_cm2'],
                "Cross_section_integrated_cm2_Hz": cross_section_data['cross_sections_integrated_cm2_Hz'],
            })
            df_E.to_csv(filename, sep="\t")
            
            # Store summary for strongest transitions
            # Filter for significant oscillator strengths
            significant = avg_osc > 0.01
            if np.any(significant):
                results_summary.append({
                    'r_angstrom': r_val,
                    'strongest_transition_eV': cross_section_data['energies_eV'][significant][0],
                    'strongest_wavelength_nm': cross_section_data['wavelengths_nm'][significant][0],
                    'max_oscillator_strength': np.max(avg_osc[significant]),
                    'max_cross_section_cm2': np.max(cross_section_data['cross_sections_peak_cm2'][significant]),
                })
            
            print(f"✓ Saved: {filename}")
    
    # Save summary
    if results_summary:
        summary_df = pd.DataFrame(results_summary)
        summary_file = "data/CaH2_oriented/Summary_strongest_transitions.csv"
        summary_df.to_csv(summary_file, sep="\t", index=False)
        print(f"\n✓ Summary saved: {summary_file}")
        print(f"\nStrongest transitions overview:")
        print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
