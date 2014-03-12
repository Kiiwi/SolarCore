#!/usr/bin/env python

"""
This script is the main part of a program written to to and model the
solar core. This project is the first term project in AST3310 at UiO
"""

# Importing packages
import matplotlib.pyplot as plt
import integration
import plots


if __name__ == "__main__":

    # Calling the integration script
    t, m, r, l, p, rho, epsilon, T0, M0, R0, \
        L0, P0, initial_rho, initial_epsilon, i_break = integration.integrate()

    # Calling the plotting script
    plots.plot_l(m, l, M0, L0, i_break)
    plots.plot_p(m, p, M0, P0, i_break)
    plots.plot_r(m, r, M0, R0, i_break)
    plots.plot_t(m, t, M0, T0, i_break)
    plots.plot_rho(m, rho, M0, initial_rho, i_break)
    plots.plot_epsilon(m, epsilon, M0, initial_epsilon, i_break)

    # Show plots
    plt.show()

# Info
__author__ = 'Tommy Ryan'
__email__ = 'tommylry@student.matnat.uio.no'
