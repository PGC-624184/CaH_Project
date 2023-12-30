import pandas as pd
import numpy as np
import unittest

def read_data(filename, skiprows=None, names=range(6)):
    try:
        df = pd.read_csv(
            filename, sep="\t", header=None, skiprows=skiprows, names=names
        )
        df.rename(
            columns={
                0: "Radius",
                1: "Broadening_x",
                2: "Broadening_y",
                3: "Broadening_z",
            },
            inplace=True,
        )
        df.Broadening_x = np.array(df.Broadening_x)
        df.Broadening_y = np.array(df.Broadening_y)
        df.Broadening_z = np.array(df.Broadening_z)
        return df
    except FileNotFoundError:
        print("File name is wrong or does not exist!.")



class TestDataReading(unittest.TestCase):
    def test_file_reading(self):
        # Test that the function successfully reads a file
        filename = "../data/Ca_H2_output.csv"  # Replace with the actual file path
        result = read_data(filename)

        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(
            len(result), 64
        )  # Replace 10 with the expected number of rows in the file

    def test_column_renaming(self):
        # Test that the function renames the columns correctly
        filename = "../data/Ca_H2_output.csv"  # Replace with the actual file path
        result = read_data(filename)

        self.assertIn("Radius", result.columns)
        self.assertIn("Broadening_x", result.columns)
        self.assertIn("Broadening_y", result.columns)
        self.assertIn("Broadening_z", result.columns)

    def test_broadening_conversion(self):
        # Test that the function converts the Broadening column to an array
        filename = "../data/Ca_H2_output.csv"  # Replace with the actual file path
        result = read_data(filename)

        self.assertTrue(isinstance(result.Broadening_x.values, np.ndarray))
        self.assertTrue(isinstance(result.Broadening_y.values, np.ndarray))
        self.assertTrue(isinstance(result.Broadening_z.values, np.ndarray))

    def test_invalid_filename(self):
        # Test that the function handles invalid filenames correctly
        filename = "nonexistent.csv"
        result = read_data(filename)

        self.assertIsNone(result)
        # You can add additional assertions to check the expected behavior of the function