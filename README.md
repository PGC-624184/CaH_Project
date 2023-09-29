# Calculating Transitional Broadening of Calcium-Hydrogen Molecule System using PySCF and TDDFT

This repository provides a Python-based guide on calculating the transitional broadening of a Calcium and Hydrogen molecule system using Time-Dependent Density Functional Theory (TDDFT) with the PySCF library.

## Prerequisites

Before you begin, make sure you have the following prerequisites installed on your system:

- [Python](https://www.python.org/downloads/)
- [PySCF](http://pyscf.org/)
- [NumPy](https://numpy.org/)
- [Matplotlib](https://matplotlib.org/)

You can install these packages using Python's package manager, `pip`. For example:

```bash
pip install pyscf numpy matplotlib
```
## Usage

Most of the analysis code is located in the src folder, with the other folders providing space to store results. To get the number densities of various elements within the M dwarf you will need to run the MARCS_Atmos.py script, which is currently set up to provide the H2 number density and mean interatomic spacing at the bottom of the photosphere.

```bash
cd src && python MARCS_Atmos.py
```

To generate the potential curves data, you can run the following:

```bash
python transitional_density_pyscf.py
```
Which will write results to the data folder, for which subsequent analysis can occur.

You can use the following to extract the relevant energy differences and oscillation strengths:

```bash
cd data && awk 'NR>=7 && NR<=9' Coarse_curve_data_Ca_H2_r6.25.csv && cd ..
```

## Contributing

If you'd like to contribute to this project, please follow the standard GitHub workflow:

1. Fork the repository.
2. Create a new branch for your feature or bug fix.
3. Make your changes and commit them.
4. Push your changes to your fork.
5. Create a pull request to merge your changes into the main repository.

## License

This project is licensed under the MIT License.

## Acknowledgements

While this repo is currently private, if it becomes public then this will be updated as applicable.