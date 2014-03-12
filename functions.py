#!/usr/bin/env python

# Importing packages
import numpy as np
import constants as const


def kappa_func(t, rho):
    """
    This function reads the opacity.txt and returns kappa
    """

    # Calculating R
    r = float(rho) / (t / 1e6)

    # Creating empty lists
    log_t = []
    log_k = []

    # Reading opacity.txt
    infile = open("opacity.txt", "r")

    # Reading the header
    log_r = infile.readline()
    log_r = log_r.split()[1:]  # Getting rid of unwanted text
    log_r = np.array(log_r, dtype=float, copy=False)

    # Skipping empty line in file
    infile.readline()

    # Appending to lists
    for line in infile:
        log_t.append(line.strip().split()[0])
        log_k.append(line.strip().split()[1:])

    # Closing file (strictly not necessary, but good habit)
    infile.close()

    # Converting lists to NumPy arrays
    log_t = np.array(log_t, dtype=float, copy=False)
    log_k = np.array(log_k, dtype=float, copy=False)

    # Finding the nearest kappa
    t_difference = abs(10 ** log_t - t)
    r_difference = abs(10 ** log_r - r)

    i = np.argmin(t_difference)
    j = np.argmin(r_difference)

    # Calculating and returning kappa
    kappa = 10 ** (log_k[i, j])
    return kappa


def energy_production(t, rho):
    """
    This function calculates energy production from the PP chains
    """

    # Number densities
    n_p = const.X * rho / const.ATOMIC_MASS_UNIT  # H
    n_he3 = const.Y3 * rho / (3 * const.ATOMIC_MASS_UNIT)  # He3
    n_he4 = const.Y * rho / (4 * const.ATOMIC_MASS_UNIT)  # He4
    n_be7 = const.Z7BE * rho / (7 * const.ATOMIC_MASS_UNIT)  # Be7
    n_e = n_p + 2 * n_he4  # e
    n_li7 = const.Z7LI * rho / (7 * const.ATOMIC_MASS_UNIT)  # Li7

    # Reaction rates
    t /= 1e9

    # TODO Fix these equations. Too ugly.
    lambd_pp = (4.01e-15 * t ** (-2.0 / 3.0)) * np.exp(-3.38 * t ** (-1.0 / 3.0)) * (
        1 + 0.123 * t ** (1.0 / 3.0) + 1.09 * t ** (2.0 / 3.0) + 0.938 * t) * const.AVOGADRO_INV
    lambd_he3he3 = (6.04e10 * t ** (-2.0 / 3.0) * np.exp(-12.276 * t ** (-1.0 / 3.0)) * (
        1 + 0.034 * t ** (1.0 / 3.0) - 0.522 * t ** (2.0 / 3.0) - 0.124 * t + 0.353 * t ** (4.0 / 3.0) + 0.213 * t ** (
            5.0 / 3.0))) * const.AVOGADRO_INV

    lamb_he3he4 = (5.61e6 * const.T1 ** (-5.0 / 6.0) * t ** (-2.0 / 3.0) * np.exp(
        -12.826 * const.T1 ** (1. / 3) * t ** (-1.0 / 3.0))) \
                  * const.AVOGADRO_INV

    # Minimum temperature dependant equation check
    if t < 1e6:
        lamb_be7e = (1.51e-7 / n_e) * const.AVOGADRO_INV
    else:
        lamb_be7e = (1.34e-10 * t ** (-1.0 / 2.0) * (1 - 0.537 * t ** (1.0 / 3.0) + 3.86 * t ** (2.0 / 3.0)
                                                     + 0.0027 * t ** (-1.0) * np.exp(
            2.515e-3 * t ** (-1.0)))) * const.AVOGADRO_INV

    lamb_li7p = (1.096e9 * t ** (-2.0 / 3.0) * np.exp(-8.472 * t ** (-1. / 3)) - 4.830e8 * const.T2 ** (-5. / 6)
                                                                                 * t ** (-2.0 / 3.0) * np.exp(
        -8.472 * const.T2 ** (1.0 / 3.0) * t ** (-1.0 / 3.0))) * const.AVOGADRO_INV

    lamb_be7p = (3.11e5 * t ** (-2.0 / 3.0) * np.exp(-10.262 * t ** (-1.0 / 3.0))) * const.AVOGADRO_INV

    # Rates per unit mass
    r_pp = lambd_pp * (n_p * n_p) / (rho * 2)
    r_he3he3 = lambd_he3he3 * (n_he3 * n_he3) / (rho * 2)
    r_he3he4 = lamb_he3he4 * (n_he3 * n_he4) / rho
    r_be7e = lamb_be7e * (n_be7 * n_e) / rho
    r_li7p = lamb_li7p * (n_p * n_li7) / rho
    r_be7p = lamb_be7p * (n_p * n_be7) / rho

    # Total energy production from all three chains
    total_energy = ((const.QHH + const.QDP) * r_pp) + (const.QHE3HE3 * r_he3he3) + (const.QHE3HE4 * r_he3he4) + (
        (const.QBE7P + const.QB8 + const.QBE8) *
        r_be7e) + (const.QLI7P * r_li7p) + (const.QBE7P * r_be7p)
    epsilon = total_energy
    return epsilon


def rho_func(t, p):
    """
    This function calculates the density
    """

    # Average mass
    mu = 1.0 / (2 * const.X + 3.0 * const.Y / 4.0 + 9.0 * const.Z / 14.0)

    # Calculate rho
    rho = mu * const.ATOMIC_MASS_UNIT / (const.SIGMA * t) * (p - const.RAD_CONST / 3.0 * t ** 4)
    return rho


# Right hand side of our equations
def rhs_drdm(r, rho):
    """
    RHS of dr/dm
    """
    rhs = 1.0 / (4.0 * np.pi * rho * r ** 2)
    return rhs


def rhs_dpdm(r, m):
    """
    RHS of dP/dm
    """
    rhs = (const.GRAV_CONST * m) / (4.0 * np.pi * r ** 4)
    return rhs


def rhs_dldm(t, rho):
    """
    RHS of dL/dm
    """
    return energy_production(t, rho)


def rhs_dtdm(t, kappa, r, l):
    """
    RHS of dT/dm
    """
    rhs = (3.0 * kappa * l) / (256.0 * np.pi ** 2.0 * const.SIGMA * r ** 4.0 * t ** 4.0)
    return rhs