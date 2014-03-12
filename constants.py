#!/usr/bin/env python

"""
This file holds constants
"""

# Useful constants
C = 3.0e10  # [cm s^-1]
GRAV_CONST = 6.67e-8  # cm^3 g^-1 s^-2]
SIGMA = 5.67e-5  # Bolzmann [erg cm^-2 s^-1 K^4]
ATOMIC_MASS_UNIT = 1.66e-24  # [g]
AVOGADRO = 6.022e23  # [mol^-1]
AVOGADRO_INV = 1.0 / AVOGADRO
MEV_TO_ERG = 1.6e-6
RAD_CONST = 4 * SIGMA / C

# Sun constants
M_SUN = 1.988e33  # [g]
L_SUN = 3.84e33  # [erg/s]
R_SUN = 6.96e10  # [cm]

# Atomic masses
MASS_H = 1.6738e-24  # [g]
MASS_HE3 = 5.0081e-24  # [g]
MASS_HE4 = 6.6464e-24  # [g]
MASS_LI7 = 1.16503486e-23  # [g]
MASS_BE7 = 1.16518851e-23  # [g]
MASS_E = 9.10938291e-28  # [g]

# PP chain energies
# PPI
QHH = 1.77 * MEV_TO_ERG
QDP = 5.49 * MEV_TO_ERG
QHE3HE3 = 12.86 * MEV_TO_ERG

# PPII
QHE3HE4 = 1.586 * MEV_TO_ERG
QBE7E = 0.049 * MEV_TO_ERG
QLI7P = 17.35 * MEV_TO_ERG

# PPIII
QBE7P = 0.137 * MEV_TO_ERG
QB8 = 8.367 * MEV_TO_ERG
QBE8 = 2.995 * MEV_TO_ERG

# Chain efficiency temperatures
T1 = 1.0495
T2 = 1.0759

# Initial ratios
X = 0.7  # H
Y3 = 1.0e-10  # He3
Y = 0.29  # Total helium
Z = 0.01  # Metals
Z7LI = 1.0e-13  # Li7
Z7BE = 1.0e-13  # Be7

# Initial conditions
L0 = 3.8e33  # [erg s^-1]
R0 = 0.5 * R_SUN  # [cm]
M0 = 0.7 * M_SUN  # [g]
T0 = 1e5  # [K]
P0 = 1e12  # [Ba]
