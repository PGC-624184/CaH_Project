import unittest
import numpy as np
from PotentialCurves import calculate_excited


class TestSum(unittest.TestCase):

    def test_calculate_excited_basic(self):
        # Test the calculate_excited function with a simple case
        bond_length = 1.2
        state_number = 2

        energy = calculate_excited(bond_length, state_number)

        # Verify that the returned energy is of the correct type (float)
        self.assertIsInstance(energy, float)

    def test_calculate_excited_ground_state(self):
        # Test the calculate_excited function for the ground state (state_id = 0)
        bond_length = 1.2
        state_number = 0

        energy = calculate_excited(bond_length, state_number)

        # Verify that the returned energy is of the correct type (float)
        self.assertIsInstance(energy, float)

        # Verify that the energy for the ground state is less than or equal to other states
        energy_state_1 = calculate_excited(bond_length, 1)
        energy_state_2 = calculate_excited(bond_length, 2)

        self.assertLessEqual(energy, energy_state_1)
        self.assertLessEqual(energy, energy_state_2)

    def test_calculate_excited_invalid_state_id(self):
        # Test the calculate_excited function with an invalid state_id
        bond_length = 1.2
        state_number = -1  # Negative state_id, which is invalid

        # Verify that the function raises a ValueError for an invalid state_id
        with self.assertRaises(ValueError):
            calculate_excited(bond_length, state_number)

# If this script is run as the main program, run the tests
if __name__ == '__main__':
    unittest.main()