import unittest
import numpy as np
import pandas as pd
import astropy.units as u
from transition_density_pyscf import run_computation 
from Prob_dist_QM import read_data, interparticle_distribution
from MARCS_Atmos import create_abundances

class TestRunComputation(unittest.TestCase):

    def test_run_computation(self):
        # Test for a specific bond length and number of states
        r = 20.0
        states = 10
        basis = 'def2-QZVPP'
        functional = 'wB97X_V'
        atom = [['Ca',0,0,0]]
        spin=0

        # Call the function
        energies, oscillator = run_computation(r, states,
                                                basis=basis,
                                                functional=functional,
                                                atom=atom,
                                                spin=spin)

        # Assert that the returned values are NumPy arrays
        self.assertIsInstance(energies, np.ndarray)
        self.assertIsInstance(oscillator, np.ndarray)

        # Add more specific assertions based on your expectations for the output
        # For example, check the shape, values, or any other relevant criteria
        self.assertEqual(energies.shape, (states,))
        self.assertEqual(oscillator.shape, (states,))

class TestDataReading(unittest.TestCase):

    def test_file_reading(self):
        # Test that the function successfully reads a file
        filename = "../Ca_H2_output.csv"  # Replace with the actual file path
        result = read_data(filename)

        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 64)  # Replace 10 with the expected number of rows in the file

    def test_column_renaming(self):
        # Test that the function renames the columns correctly
        filename = "../Ca_H2_output.csv"  # Replace with the actual file path
        result = read_data(filename)

        self.assertIn("Radius", result.columns)
        self.assertIn("Broadening", result.columns)

    def test_broadening_conversion(self):
        # Test that the function converts the Broadening column to an array
        filename = "../Ca_H2_output.csv"  # Replace with the actual file path
        result = read_data(filename)

        self.assertTrue(isinstance(result.Broadening.values, np.ndarray))

    def test_invalid_filename(self):
        # Test that the function handles invalid filenames correctly
        filename = "nonexistent.csv"
        result = read_data(filename)

        self.assertIsNone(result)
        # You can add additional assertions to check the expected behavior of the function


class TestInterparticleDistribution(unittest.TestCase):

    def test_interparticle_distribution(self):
        # Test the interparticle distribution calculation

        # Define test inputs
        r = 5 * u.angstrom
        N = 100 * u.cm**(-3)

        # Calculate the expected result manually
        expected_result = (N*4*np.pi*r**2*np.exp(-(4/3)*np.pi*r**3*N)).to(u.nm**(-1))

        # Call the function under test
        result = interparticle_distribution(r, N)

        # Assert that the calculated result matches the expected result
        self.assertAlmostEqual(result.value, expected_result.value)
        self.assertEqual(str(result.unit), str(expected_result.unit))

class TestCreateAbundances(unittest.TestCase):

    def test_create_abundances(self):
        # Test the creation of elemental abundances DataFrame

        # Define test inputs
        file = "../reference_stellar_atmospheres/t2700_g+5.5_z+0.00.mod"

        # Call the function under test
        result = create_abundances(file)

        # Check if the result is a pandas DataFrame
        self.assertIsInstance(result, pd.DataFrame)

        # Check if the result has correct shape
        self.assertEqual(result.shape, (1, 92))

        # Check some specific values in the DataFrame
        #self.assertAlmostEqual(result.at[0, 'O'], expected_value_for_O)
        #self.assertAlmostEqual(result.at[0, 'Fe'], expected_value_for_Fe)
        # Add more assertions for other atomic symbols and expected values

if __name__ == '__main__':
    unittest.main()