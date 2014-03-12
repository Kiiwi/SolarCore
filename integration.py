#!/usr/bin/env python

import numpy as np

import constants as const

import time

import sys

import functions as func


def integrate():
    """
    This function integrates the equations
    governing the structure of the Sun
    """

    # Get user input for integration steps
    #n = int(input("N: "))
    n = int()

    # Creating arrays
    t = np.zeros(n)
    m = np.zeros(n)
    r = np.zeros(n)
    l = np.zeros(n)
    p = np.zeros(n)
    rho = np.zeros(n)
    epsilon = np.zeros(n)

    # Getting initial values
    initial_rho = func.rho_func(const.T0, const.P0)
    initial_epsilon = const.L0 / const.M0

    # Adding initial values to arrays
    t[0] = const.T0
    m[0] = const.M0
    r[0] = const.R0
    l[0] = const.L0
    p[0] = const.R0
    rho[0] = initial_rho
    epsilon[0] = initial_epsilon

    # Array for increasing dm steadily by splitting N into sections
    k = 3
    i_array = [j * n / (k + 1) for j in range(1, k + 2)]
    dm_array = [1.4e23 * 10 ** j for j in range(1, len(i_array) + 1)]

    #start = time.time()
    for i in range(n - 1):
        #print(r[i]) # For debugging
        if t[i] <= 0 or m[i] <= 0 or r[i] <= 0 or l[i] <= 0 or p[i] <= 0 or rho[i] <= 0:
            i_break = i
            print("Something went to zero! Debug to find out what happened.")
            #print(r[i]) # For debugging
            break
        else:
            # Steady increase of dm
            for x in range(k+1):
                if i < i_array[x]:
                    dm = dm_array[x]

            # Forward Euler method
            t[i + 1] = t[i] - func.rhs_dtdm(t[i], func.kappa_func(t[i], rho[i]), r[i], l[i])
            m[i + 1] = m[i] - dm
            r[i + 1] = r[i] + func.rhs_drdm(r[i], rho[i]) * dm
            l[i + 1] = l[i] + func.rhs_dldm(t[i], rho[i]) * dm
            p[i + 1] = p[i] - func.rhs_dpdm(r[i], m[i]) * dm
            rho[i + 1] = func.rho_func(t[i], p[i])
            epsilon[i + 1] = func.energy_production(t[i], rho[i])

    return t, m, r, l, p, rho, epsilon, const.T0, const.M0, const.R0, \
           const.L0, const.P0, initial_rho, initial_epsilon, i_break


