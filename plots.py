#!/usr/bin/env python

import matplotlib.pyplot as plt


def plot_r(m, r, m0, r0, i_break):
    """
    Function to plot r v m
    """
    fig = plt.figure()
    fig.add_subplot(111)
    plt.title("Position")
    plt.xlabel("$m/m0$")
    plt.ylabel("$r/r0$")
    plt.plot(m[0:i_break]/m0, r[0:i_break]/r0)
    plt.grid("on")
    fig.tight_layout()


def plot_p(m, p, m0, p0, i_break):
    """
    Function to plot p v m
    """
    fig = plt.figure()
    fig.add_subplot(111)
    plt.title("Pressure")
    plt.xlabel("$m/m0$")
    plt.ylabel("$p/p0$")
    plt.plot(m[0:i_break]/m0, p[0:i_break]/p0)
    plt.grid("on")
    fig.tight_layout()


def plot_l(m, l, m0, l0, i_break):
    """
    Function to plot l v m
    """
    fig = plt.figure()
    fig.add_subplot(111)
    plt.title("Luminosity")
    plt.xlabel("$m/m0$")
    plt.ylabel("$L/l0$")
    plt.plot(m[0:i_break]/m0, l[0:i_break]/l0)
    plt.grid("on")
    fig.tight_layout()


def plot_t(m, t, m0, t0, i_break):
    """
    Function to plot t v m
    """
    fig = plt.figure()
    fig.add_subplot(111)
    plt.title("Temperature")
    plt.xlabel("$m/m0$")
    plt.ylabel("$T/t0$")
    plt.plot(m[0:i_break]/m0, t[0:i_break]/t0)
    plt.grid("on")
    fig.tight_layout()

# Unused plotting functions


def plot_rho(m, rho, m0, rho_init, i_break):
    """
    Function to plot rho v m
    """
    fig = plt.figure()
    fig.add_subplot(111)
    plt.title("Density")
    plt.xlabel("$m/M0$")
    plt.ylabel("$\\rho/\\rho_init$")
    plt.plot(m[0:i_break]/m0, rho[0:i_break]/rho_init)
    plt.grid("on")
    fig.tight_layout()


def plot_epsilon(m, epsilon, m0, epsilon_init, i_break):
    """
    Function to plot epsilon v m
    """
    fig = plt.figure()
    fig.add_subplot(111)
    plt.title("Density")
    plt.xlabel("$m/M0$")
    plt.ylabel("$\\epsilon")
    plt.plot(m[0:i_break]/m0, epsilon[0:i_break]/epsilon_init)
    plt.grid("on")
    fig.tight_layout()