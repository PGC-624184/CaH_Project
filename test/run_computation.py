import unittest
import numpy as np
from transition_density_pyscf import run_computation

class TestRunComputation(unittest.TestCase):
    def test_run_computation(self):
        # Test for a specific bond length and number of states
        r = 20.0
        states = 10
        basis = "def2-QZVPP"
        functional = "wB97X_V"
        atom = [["Ca", 0, 0, 0]]
        spin = 0

        # Call the function
        energies, oscillator = run_computation(
            r, states, basis=basis, functional=functional, atom=atom, spin=spin
        )

        # Assert that the returned values are NumPy arrays
        self.assertIsInstance(energies, np.ndarray)
        self.assertIsInstance(oscillator, np.ndarray)

        # Add more specific assertions based on your expectations for the output
        # For example, check the shape, values, or any other relevant criteria
        self.assertEqual(energies.shape, (states,))
        self.assertEqual(oscillator.shape, (states,))