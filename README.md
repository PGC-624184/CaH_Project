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
cd src &&python transitional_density_pyscf.py
```
Which will write results to the data folder, for which subsequent analysis can occur.

You can use the following to extract the relevant energy transition differences (_column 2 between 0.09 and 0.12_) and oscillation strengths (_column 3 greater than 0.4_) to briefly check the results:

```bash
echo -e "Delta_E \t\t f" && awk '$3>=0.4 && $2>=0.09 && $2<= 0.12 {print $2,$3}' Coarse_curve_data_Ca_H2_r16.25.csv
```

Or run the ./CaH_extract.sh and the ./CaH2_extract.sh to do this automatically and create csv files with the final results,

```bash
cd src && ./CaH_extract.sh
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

- PySCF: the Python-based simulations of chemistry framework, Q. Sun, T. C. Berkelbach, N. S. Blunt, G. H. Booth, S. Guo, Z. Li, J. Liu, J. McClain, S. Sharma, S. Wouters, and G. K.-L. Chan, WIREs Comput. Mol. Sci. 8, e1340 (2018)

- Spencer Collivati - for providing assistance developing this computational suite

- Mike Ireland - for providing basis of the problem and providing relevant stellar atmosphere models

- Thomas Nordlander - for providing relevant data files and providing advice on what is physical for the system

