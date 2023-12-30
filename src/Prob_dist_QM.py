import numpy as np
import pandas as pd
import astropy.units as u
import matplotlib.pyplot as plt


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


def interparticle_distribution(r, N):
    """
    Calculate the interparticle distribution for a given system.

    This function computes the interparticle distribution for a system with a given
    radial distance 'r' and a specified number of particles 'N'. It applies a
    mathematical formula and returns the result in units of inverse angstroms.

    Args:
        r (astropy.Quantity): The radial distance at which to calculate the distribution.
        N (int): The number of particles in the system.

    Returns:
        astropy.Quantity: The interparticle distribution in inverse angstroms.

    Note:
        - This function is designed for a specific mathematical calculation and requires
          the 'r' parameter to be an astropy Quantity with units (e.g., u.angstrom).
        - Ensure 'r' is correctly defined with units before using this function.
    """
    return (N * 4 * np.pi * r**2 * np.exp(-(4 / 3) * np.pi * r**3 * N)).to(
        u.nm ** (-1)
    )


if __name__ == "__main__":
    N_H2 = 1.0381111448319373e19 * u.cm ** (-3)
    n_H2 = [
        0.00000000e00,
        5.59889775e15,
        8.48542416e15,
        1.22743080e16,
        1.77580194e16,
        2.51111075e16,
        3.47064969e16,
        4.90940600e16,
        6.78687469e16,
        9.60089352e16,
        1.32690066e17,
        1.83322069e17,
        2.17948504e17,
        2.59065128e17,
        3.00857396e17,
        3.57448356e17,
        4.24569502e17,
        4.92688134e17,
        5.84850521e17,
        6.93988548e17,
        8.04415905e17,
        9.53666569e17,
        1.10435816e18,
        1.30788699e18,
        1.51269135e18,
        1.78898303e18,
        2.06602807e18,
        2.38387496e18,
        2.74635366e18,
        3.16408894e18,
        3.53436393e18,
        3.95431203e18,
        4.43268898e18,
        4.86516028e18,
        5.34750273e18,
        6.02209809e18,
        6.48309852e18,
        7.14707110e18,
        7.88470648e18,
        8.70445372e18,
        9.39500055e18,
        1.03811114e19,
        1.12136621e19,
        1.21171978e19,
        1.34026059e19,
        1.44903001e19,
        1.56702843e19,
        1.69490383e19,
        1.83355176e19,
        1.98395509e19,
        2.14694426e19,
        2.37777859e19,
        2.78599006e19,
        3.26533170e19,
        3.82783973e19,
        4.38568518e19,
    ] * u.cm ** (-3)

    n_H = [
        0.00000000e00,
        2.44400751e12,
        4.35184638e12,
        7.74456481e12,
        1.34708353e13,
        2.29016022e13,
        3.80549171e13,
        6.32454008e13,
        1.05116243e14,
        1.70730713e14,
        2.77228842e14,
        4.50002318e14,
        5.73262973e14,
        7.30144736e14,
        9.29738237e14,
        1.18362275e15,
        1.54151869e15,
        1.96142679e15,
        2.49485254e15,
        3.24602831e15,
        4.12553854e15,
        5.36286122e15,
        6.96802890e15,
        9.04836724e15,
        1.17422230e16,
        1.55813908e16,
        2.06602807e16,
        2.73705467e16,
        3.62039921e16,
        4.78904241e16,
        6.89146042e16,
        9.48573386e16,
        1.27840147e17,
        1.68693035e17,
        2.17846715e17,
        2.75262991e17,
        3.48163004e17,
        4.20850747e17,
        5.20937391e17,
        6.30581929e17,
        7.46271420e17,
        9.04157009e17,
        1.07089642e18,
        1.26882638e18,
        1.46956651e18,
        1.74211724e18,
        2.01872368e18,
        2.33961857e18,
        2.71202178e18,
        3.14435692e18,
        3.56304058e18,
        4.13210560e18,
        5.43224767e18,
        6.98115538e18,
        8.97333220e18,
        1.15355274e19,
    ] * u.cm ** (-3)

    r_Atmos = [
        -5.0,
        -4.8,
        -4.6,
        -4.4,
        -4.2,
        -4.0,
        -3.8,
        -3.6,
        -3.4,
        -3.2,
        -3.0,
        -2.9,
        -2.8,
        -2.7,
        -2.6,
        -2.5,
        -2.4,
        -2.3,
        -2.2,
        -2.1,
        -2.0,
        -1.9,
        -1.8,
        -1.7,
        -1.6,
        -1.5,
        -1.4,
        -1.3,
        -1.2,
        -1.1,
        -1.0,
        -0.9,
        -0.8,
        -0.7,
        -0.6,
        -0.5,
        -0.4,
        -0.3,
        -0.2,
        -0.1,
        -0.0,
        0.1,
        0.2,
        0.3,
        0.4,
        0.5,
        0.6,
        0.7,
        0.8,
        0.9,
        1.0,
        1.2,
        1.4,
        1.6,
        1.8,
        2.0,
    ]
    N_H = 9.041570094762417e17 * u.cm ** (-3)
    factor = 45.56335252907954
    file_name = "data/Ca_H2_test.csv"
    file_name1 = "data/Ca_H_test.csv"
    H2_df = read_data(file_name)
    H_df = read_data(file_name1, skiprows=6)

    avg_broaden_H = []
    avg_broaden_H2 = []

    for depth in n_H2:
        p_H2 = []
        radial = H2_df.Radius * u.angstrom
        for val in radial:
            p_H2.append(interparticle_distribution(val, depth).value)
        average_x = np.sum(H2_df.Broadening_x * p_H2)
        average_y = np.sum(H2_df.Broadening_y * p_H2)
        average_z = np.sum(H2_df.Broadening_z * p_H2)
        # print("The average Ca H2 x Broadening is: {:2.3}".format(average_x))
        # print("The average Ca H2 y Broadening is: {:2.3} ".format(average_y))
        # print("The average Ca H2 z Broadening is: {:2.3} ".format(average_z))
        # print("The average Ca H2 Broadening is: {:2.6} ".format(np.sqrt(average_x**2+average_y**2+average_z**2)))
        avg_broaden_H2.append(np.sqrt(average_x**2 + average_y**2 + average_z**2))

    for depth in n_H:
        p_H = []
        radial2 = H_df.Radius * u.angstrom
        for val in radial2:
            p_H.append(interparticle_distribution(val, depth).value)

        H_average_x = np.sum(H_df.Broadening_x * p_H)
        H_average_y = np.sum(H_df.Broadening_y * p_H)
        H_average_z = np.sum(H_df.Broadening_z * p_H)
        # print("The average Ca H x Broadening is: {:2.3}".format(H_average_x))
        # print("The average Ca H y Broadening is: {:2.3} ".format(H_average_y))
        # print("The average Ca H z Broadening is: {:2.3} ".format(H_average_z))
        # print("The average Ca H Broadening is: {:2.5} ".format(np.sqrt(H_average_x**2+H_average_y**2+H_average_z**2)))
        avg_broaden_H.append(
            np.sqrt(H_average_x**2 + H_average_y**2 + H_average_z**2)
        )

    plt.clf()
    plt.semilogy(r_Atmos, avg_broaden_H2, "x", label=r"Ca $H_2$")
    plt.semilogy(r_Atmos, avg_broaden_H, "+", label=r"Ca $H$")
    plt.xlabel(r"$log_{10}(\tau_R)$ Optical Depth")
    plt.ylabel(r"Average Quasi-Static Broadening $(nm)$")
    plt.title(
        "Probability Weighted Quasi-Static Calcium Broadening\n for a Proxima Centauri-like Atmosphere"
    )
    plt.legend()
    plt.show()
    print(
        "The average Ca H2 Continuum Broadening is: {:2.5} nm ".format(
            np.mean(np.array(avg_broaden_H2[-16::]))
        )
    )
    print(
        "The average Ca H Continuum Broadening is: {:2.5} nm".format(
            np.mean(np.array(avg_broaden_H[-16::]))
        )
    )
