'''
This python code is used to draw the graph of isentropic Equations of states, compare EPAW with NIF.
author: jz2907@columbia.edu
'''

from qha.unit_conversion import *
from scipy.optimize import curve_fit
import numpy as np
from numpy.linalg import inv
import pandas as pd
import os
import matplotlib.pyplot as plt
import palettable
import time
import random

from isothermeos import poly_fit, poly_calc

from matplotlib.axes._axes import _log as matplotlib_axes_logger
matplotlib_axes_logger.setLevel('ERROR')

# Color initialize
cm = palettable.cartocolors.qualitative.Prism_10
color_list = cm.mpl_colors
color_list.insert(0, "black")

pcut = 100 # set of pcut

def readcsv_t0_p0(t0,p0, folderName="data/epaw0"):
    '''
    read the compression curve data, default folder '/data/epaw0'
    return density, pressure, in unit of g/cm^3 and GPa respectively.
    '''

    fileName = "%.1f-%.1f.csv"%(t0,p0)
    df = pd.read_csv(os.path.join(folderName, fileName))
    return b3_to_density(df["v"]), df["p"]

def plot_nif(nif):
    '''
    Draw NIF Data on the current figure
    '''
    plt.errorbar(nif["d"], nif["p"], xerr=nif["d_dev"], yerr=nif["p_dev"], c=color_list[8], marker=".", mec=color_list[0], ecolor=color_list[6], label="NIF")

def plot_nif_0(nif):
    '''
    Draw NIF Data on 0
    '''
    plt.errorbar(nif["d"], nif["p"]* 0, xerr=nif["d_dev"], yerr=nif["p_dev"],
                 c=color_list[8], marker=".", mec=color_list[0], ecolor=color_list[6], label="NIF")

def poly_func(x,a,b,c,d,e):
    return a + b * x + c * x ** 2 + d * x ** 3 + e * x ** 4

def poly_curve_fit(d,p):
    return curve_fit(poly_func, d, p)

def poly_simu(t0,p0):
    '''
    give the initial conditions, t0, po
    return the polynomial 4th order fitting parameters and the density range
    '''
    d,p = readcsv_t0_p0(t0, p0)
    param, _ = poly_curve_fit(d, p)
    lim = (min(d),max(d))
    return param, lim


def deviation(nif, param, lim, pcut, d=[]):
    '''
    calculated dP/P on the nif's densitys, only valid of p > pcut
    return d, pratio
    '''
    if d == []:
        d = nif["d"]

    new_p = poly_func(d,*param)
    p = nif["p"]
    p_dev = nif["p_dev"]
    p_diff = (new_p - p)
    p_ratio = (new_p - p)/p_dev

    # discard outside points
    # new_p_cut = [i for i in new_p if i>= pcut]
    new_p_cut = new_p[-54:]

    return d[-len(new_p_cut):], p_diff[-len(new_p_cut):], p_ratio[-len(new_p_cut):]

def mc_uncertainty(nif,param,lim,pcut,nmc = 5000):
    '''
    '''
    d_nif = nif["d"]
    p_nif = nif["p"]
    d_dev = nif["d_dev"]
    p_dev = nif["p_dev"]

    p_diff = np.zeros((54, nmc))
    p_ratio = np.zeros((54,nmc))
    for i in range(nmc):
        d_mc = [random.normalvariate(d_nif[j] ,d_dev[j]/2 )  for j in range(len(d_nif))]
        d_mc = np.array(d_mc)
        d1, pd1, pr1 = deviation(nif, param, lim, pcut, d=d_mc)
        p_diff[:, i] = pd1
        p_ratio[:,i] = pr1
        
    return d_nif[-len(d1):], np.average(p_diff, axis=1), np.average(p_ratio, axis=1)


def plot_poly_valid(nif, t0, p0, pcut):
    '''
    Figure 0. polynomial fit!
    '''
    param, lim = poly_simu(t0, p0)
    d,p = readcsv_t0_p0(t0, p0)
    pp = poly_func(d,*param)
    
    fig=plt.figure(figsize=(5,6))
    plt.subplots_adjust(left=0.15, bottom=None, right=0.95, top=None,
                        wspace=None, hspace=0)

    plt.subplot2grid((5, 1), (0, 0), rowspan=4)

    plt.plot(d, p, label="Simulation at $T_0$ = %s, $P_0$ = %s" %
             (t0, p0), lw=2.5, alpha=0.6)
    plt.plot(d, pp, label= "Polynomial 4th order fitting curve")
    plt.legend(frameon=False)
    plt.title("Validation of Polynomial 4th order fitting\nfor EPAW simulation")
    plt.xlim(9.8, 21)
    plt.ylim(0, 1460)
    plt.xlabel("Density (g/cm$^3$)")
    plt.ylabel("Pressure (GPa)")

    left, bottom, width, height = 0.3, 0.55, 0.25, 0.21
    axn = fig.add_axes([left, bottom, width, height])
    plt.plot(d, p, label="Simulation at $T_0$ = %s, $P_0$ = %s" %
             (t0, p0), lw=2.5, alpha=0.6)
    plt.plot(d, pp, label="Polynomial 4th order fitting curve")

    axn.set_xlabel("Density (g/cm$^3$)")
    axn.set_ylabel("Pressure (GPa)")

    axn.set_xlim([18.6, 18.8])
    axn.set_ylim([1000, 1050])

    plt.subplot2grid((5, 1), (4, 0), rowspan=1)

    plt.plot(d, p - pp, lw=2, alpha = 0.6)
    plt.plot(d, pp - pp)

    plt.xlim(9.8, 21)
    plt.xlabel("Density (g/cm$^3$)")
    plt.ylabel("$\\Delta P$ (GPa)")

    plt.savefig("graph/polyverify.pdf")

def plot_nif_epaw(nif,t0,p0,pcut):
    '''
    Figure 1. NIF with uncertainty and a single EPAW simulation
    '''
    plt.figure(figsize=(5, 7))
    plt.subplots_adjust(left=0.15, bottom=None, right=0.95, top=None,
                        wspace=None, hspace=0)
    
    plt.subplot2grid((6, 1), (0, 0), rowspan=4)
    plot_nif(nif)

    param,lim = poly_simu(t0, p0)
    x = np.linspace(*lim, 100)
    plt.plot(x,poly_func(x, *param), color=color_list[3], label="Simulation at $T_0$ = %s, $P_0$ = %s" % (t0, p0),lw=2)

    plt.hlines(pcut,0,25,colors=color_list[5], ls=":",label="cut off pressure")
    plt.legend(frameon=False)
    plt.title("View of NIF data and\n one isentropic EPAW simulation")
    plt.xlim(9.8, 21)
    plt.ylim(0, 1460)
    plt.xlabel("Density (g/cm$^3$)")
    plt.ylabel("Pressure (GPa)")

    dx, p_diff, p_ratio = deviation(nif, param, lim, pcut)

    plt.subplot2grid((6, 1), (4, 0), rowspan=1)
    plot_nif_0(nif)
    plt.plot(dx, p_diff, label="Simulation at $T_0$ = %s, $P_0$ = %s" %
             (t0, p0), lw=2, color=color_list[3])
    # plt.legend(frameon=False)
    plt.xlim(9.8, 21)
    plt.xlabel("Density (g/cm$^3$)")
    plt.ylabel("$\\Delta P$ (GPa)")

    plt.subplot2grid((6, 1), (5, 0), rowspan=1)
    plt.plot(dx, p_ratio, label="Simulation at $T_0$ = %s, $P_0$ = %s" %
             (t0, p0), lw=2, color=color_list[3])
    # plt.legend(frameon=False)
    plt.xlim(9.8, 21)
    plt.xlabel("Density (g/cm$^3)$")
    plt.ylabel("$\\frac{P_{simu}-P_{NIF}}{\\sigma_{P}}$ (%)")
    plt.savefig("graph/NIF.pdf", dpi=300)


def plot_nif_epaw_mc(nif, t0, p0, pcut):
    '''
    Figure 2. NIF with uncertainty and a single EPAW simulation, 
    the bottom figure is using mc_uncertainty
    '''
    plt.figure(figsize=(5, 7))
    plt.subplots_adjust(left=0.15, bottom=None, right=0.95, top=None,
                        wspace=None, hspace=0)

    plt.subplot2grid((6, 1), (0, 0), rowspan=4)
    plot_nif(nif)

    param, lim = poly_simu(t0, p0)
    x = np.linspace(*lim, 100)
    plt.plot(x, poly_func(x, *param),
             color=color_list[3], label="Simulation at $T_0$ = %s, $P_0$ = %s" % (t0, p0), lw=2)

    plt.hlines(
        pcut, 0, 25, colors=color_list[5], ls=":", label="cut off pressure")
    plt.legend(frameon=False)
    plt.title("View of NIF data and\n one isentropic EPAW simulation")
    plt.xlim(9.8, 21)
    plt.ylim(0, 1460)
    plt.xlabel("Density (g/cm$^3$)")
    plt.ylabel("Pressure (GPa)")

    dx, p_diff, p_ratio = mc_uncertainty(nif, param, lim, pcut)

    plt.subplot2grid((6, 1), (4, 0), rowspan=1)
    plot_nif_0(nif)
    plt.plot(dx, p_diff, label="Simulation at $T_0$ = %s, $P_0$ = %s" %
             (t0, p0), lw=2, color=color_list[3])
    # plt.legend(frameon=False)
    plt.xlim(9.8, 21)
    plt.xlabel("Density (g/cm$^3$)")
    plt.ylabel("$\\Delta P$ (GPa)")

    plt.subplot2grid((6, 1), (5, 0), rowspan=1)
    plt.plot(dx, p_ratio, label="Simulation at $T_0$ = %s, $P_0$ = %s" %
             (t0, p0), lw=2, color=color_list[3])
    # plt.legend(frameon=False)
    plt.xlim(9.8, 21)
    plt.xlabel("Density (g/cm$^3)$")
    plt.ylabel("$\\frac{P_{simu}-P_{NIF}}{\\sigma_{P}}$ (%)")
    plt.savefig("graph/NIF_mc_1.pdf", dpi=300)


def plot_all_mc(nif,t0_list,p0_list,pcut):

    p_ratio_matrix = np.zeros((len(t0_list), len(p0_list)))
    for i,t0 in enumerate(t0_list):
        for j,p0 in enumerate(p0_list):
            param, lim = poly_simu(t0, p0)
            dx, p_diff, p_ratio = mc_uncertainty(nif, param, lim, pcut,nmc=1000)
            p_ratio_ttl = np.sum(np.abs(p_ratio))
            p_ratio_matrix[i, j] = p_ratio_ttl

    plt.close()
    plt.figure()
    plt.imshow(p_ratio_matrix, cmap="rainbow")
    plt.colorbar()

    init_t_list = np.linspace(500, 2000, 4)
    init_p_list = np.linspace(50, 70, 5)
    yrange = np.linspace(0, 15, 4)
    xrange = np.linspace(0, 20, 5)
    plt.xticks(xrange, init_p_list)
    plt.yticks(yrange, init_t_list)

    plt.title("Deviation map")

    best_fit = np.where(p_ratio_matrix == np.min(p_ratio_matrix))
    print(best_fit)

    plt.scatter(10, 5, marker="*", s=80, color="k")
    plt.scatter(best_fit[1], best_fit[0], marker="*", s=80, color="white")

    plt.xlabel('Initial Pressure (GPa)')
    plt.ylabel('Initial Temperature (K)')

    plt.savefig("graph/nif_all_1000.pdf")

if __name__=="__main__":

    # Load nif data
    nif = pd.read_csv(os.path.join("data", "nif.csv"))

    # Draw figure 1
    # plot_nif_epaw(nif,1000,60,100)

    # Draw figure 2: Do the monte carlo uncertainty test on a single curve
    # plot_nif_epaw_mc(nif,1000,60,100)

    # valit of polynomial 4th order fit of epaw!
    plot_poly_valid(nif, 1000, 60, 100)

    # Do the whole thing
    # t0_list = np.linspace(500, 2000, 16)
    # p0_list = np.linspace(50, 70, 21)
    # plot_all_mc(nif,t0_list,p0_list,100)

    
