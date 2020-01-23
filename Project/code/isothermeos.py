'''
This python code is used to draw the graph of isothermal Equations of states.
author: jz2907@columbia.edu
'''

from scipy.optimize import curve_fit
import numpy as np
from numpy.linalg import inv
import pandas as pd
import os
import matplotlib.pyplot as plt
import palettable

from matplotlib.axes._axes import _log as matplotlib_axes_logger
matplotlib_axes_logger.setLevel('ERROR')

# Color initialize
cm = palettable.cartocolors.qualitative.Prism_10
color_list = cm.mpl_colors
color_list.insert(0,"black")
color_exp = color_list[1]
color_ab = color_list[5]

def readcsv(fileName, folderName="data"):
    '''
    read the compression curve data, default folder '/data'
    return volume, pressure, in unit of A3 and GPa respectively.
    '''
    df = pd.read_csv(os.path.join(folderName, fileName))
    v = df['v']
    p = df['p']
    return v, p

def vinet_p(v, k0, kp0, v0):
    '''
    P(V) of Vinet equations of state.
    return pressures 
    '''
    x = (v / v0) ** (1 / 3)
    xi = 3 / 2 * (kp0 - 1)
    return 3 * k0 / x ** 2 * (1 - x) * np.exp(xi * (1 - x))

def vp_vinet(k0, kp0, v0):
    '''
    give the fitting parameters k0, kp0, and v0, 
    the default setting of volume is 6 ~ 11 A3
    return the volumes and Vinet pressures.
    '''
    v = np.linspace(6, 11, 100)
    return v, vinet_p(v, k0, kp0, v0)

def bm3_p(v, k0, kp0, v0):
    '''
    P(V) of Birch-Murnaghan 3rd order equations of state. 
    return pressures 
    '''
    eta = (v0 / v) ** (1 / 3)
    return 3 / 2 * k0 * (eta ** 7 - eta ** 5) * (1 + 3 / 4 * (kp0 - 4) * (eta ** 2 - 1))

def vp_bm3(k0, kp0, v0):
    '''
    give the fitting parameters k0, kp0, and v0, 
    the default setting of volume is 6 ~ 11 A3
    return the volumes and Birch-Murnaghan 3rd order pressures.
    '''
    v = np.linspace(6,11,100)
    return v, bm3_p(v,k0,kp0,v0)

def eos_fit(v: np.ndarray, p: np.ndarray):
    '''
    return popt, pcov of the fitting parameters;
    popt is (k0, kp0, v0)
    '''
    return curve_fit(bm3_p, v, p, p0=[200, 4, 10])


def poly_calc(x, param):
    '''
    return the poly 2nd order result: a + b x + c x^2
    '''
    ans = 0 * x
    for i,a in enumerate([*param]):
        ans += a * x ** i
    return ans

def poly_fit(v, p, n = 4):
    '''
    polynomial fitting a + b x + c x^2 + ...
    default n = 4
    return a,b,c,...
    '''

    G = np.ones((len(v),n+1))
    for i in range(n):
        G[:,i+1] = v **(i+1)
    Gt = np.transpose(G)
    invGtG = inv(np.matmul(Gt,G))
    tov = np.matmul(invGtG, Gt)
    mest = np.matmul(tov, p)

    # calculate the sig_m
    e = p - np.matmul(G, mest)
    et = np.transpose(e)
    E = np.matmul(et, e)
    sigma2d = E / (len(p) - (n+1))
    Cm = sigma2d * invGtG
    sigma2m = np.diag(Cm)
    # print(sigma2m**0.5,mest)

    return mest, sigma2d**0.5
    
def plot_raw_data(data):
    '''
    draw all the data passing into this function
    The form of data should be:
    data = [ dict1, dict2, ... ]
    and the dicts should be:
    dict = {
            "v" = [...],
            "p" = [...],
            "draw" = "plot" or "scatter",
            "color" = "...",
            "label" = "...",
            "ls" = "the_line_style"
    }
    '''

    plt.figure(figsize=(6,5))

    for item in data:
        if item["draw"] == "plot":
            plt.plot(item["v"],item["p"],c=item["color"],label=item["label"],lw=2,ls=item["ls"])
        elif item["draw"] == "scatter":
            plt.scatter(item["v"], item["p"], c=item["color"], label=item["label"], lw=1,s=60, marker=item["ls"],edgecolors="#2B2B2B",alpha=0.8)

    plt.xlabel('Volume / Atom (Å$^3$)')
    plt.ylabel('Pressure (GPa)')

    plt.xlim(6, 11.2)
    plt.ylim(-10, 400)
    
    plt.title("Isothermal equations of state")

    plt.legend(frameon=False)
    plt.savefig("graph/300K.pdf", dpi=300)

def plot_the_data(data):
    data_p = [(i["p"]) for i in data]
    data_p = np.concatenate(data_p)
    data_v = [(i["v"]) for i in data]
    data_v = np.concatenate(data_v)

    if data[0]["type"] == "exp":
        label = "Experimental data"
        color = color_exp
    else:
        label = "Ab initio data"
        color = color_ab
    plt.scatter(data_v, data_p,
                marker="o", label=label, c=color, alpha = 0.2)

def plot_poly_fitting(data):
    '''
    '''
    data_p = [(i["p"]) for i in data]
    data_p = np.concatenate(data_p)
    data_v = [(i["v"]) for i in data]
    data_v = np.concatenate(data_v)

    # calculate the bm3 fit parameters of data
    data_bm3_fit, _ = eos_fit(data_v, data_p)
    data_bm3_fit_v, data_bm3_fit_p = vp_bm3(*data_bm3_fit)

    # calculate the 4nd order poly fit of the origional data set
    data_poly_fit,data_poly_fit_sig = poly_fit(data_v, data_p)
    data_poly_fit_p = poly_calc(
        data_bm3_fit_v, data_poly_fit)
    
    # calculate the 4nd order polynomial fit using least sqare solution of the bm3 fitted curve
    data_bm3_fit_poly_fit, data_bm3_fit_poly_fit_sig = poly_fit(
        data_bm3_fit_v, data_bm3_fit_p)
    data_bm3_fit_poly_fit_p = poly_calc(
        data_bm3_fit_v, data_bm3_fit_poly_fit)

    # plot of fitting curves
    
    plt.plot(data_bm3_fit_v, data_bm3_fit_p -
             data_bm3_fit_p, label="BM3 fitting",lw=2)
    plt.plot(data_bm3_fit_v, data_bm3_fit_poly_fit_p -
             data_bm3_fit_p, label="Poly 4th fitting of BM3 curve",lw=2)
    plt.plot(data_bm3_fit_v, data_poly_fit_p-data_bm3_fit_p,
             label="Poly 4th fitting of original data",lw=2)
    plt.scatter(data_v, data_p-bm3_p(data_v, *data_bm3_fit),
                marker=".", label="Original data",c="#BABABA")

    # plt.plot(data_bm3_fit_v, (data_bm3_fit_p -
    #                           data_bm3_fit_p)/data_bm3_fit_p, label="BM3 fitting", lw=2)
    # plt.plot(data_bm3_fit_v, (data_bm3_fit_poly_fit_p -
    #                           data_bm3_fit_p)/data_bm3_fit_p, label="Poly 4th fitting of BM3 curve", lw=2)
    # plt.plot(data_bm3_fit_v, (data_poly_fit_p-data_bm3_fit_p)/data_bm3_fit_p,
    #          label="Poly 4th fitting of original data", lw=2)
    # plt.scatter(data_v, (data_p-bm3_p(data_v, *data_bm3_fit))/bm3_p(data_v, *data_bm3_fit),
    #             marker=".", label="Original data", c="#BABABA")

    plt.xlabel('Volume / Atom (Å$^3$)')
    plt.ylabel('% $\\Delta$Pressure (GPa)')
    # plt.ylabel(' $\\Delta$P/ P (%)')
    plt.xlim(6, 11)
    if data[0]["type"] =="exp":
        plt.title("In situ experimental data")
    else:
        plt.title("Ab initio calculations")
    plt.legend(frameon=False)
    
    return data_poly_fit, data_poly_fit_sig


def plot_poly_fitting_veri(data):
    '''
    '''
    data_p = [(i["p"]) for i in data]
    data_p = np.concatenate(data_p)
    data_v = [(i["v"]) for i in data]
    data_v = np.concatenate(data_v)

    # calculate the bm3 fit parameters of data
    data_bm3_fit, _ = eos_fit(data_v, data_p)
    data_bm3_fit_v, data_bm3_fit_p = vp_bm3(*data_bm3_fit)

    # calculate the 2nd order poly fit of the origional data set
    data_poly_2_fit, _ = poly_fit(data_v, data_p,n=2)
    data_poly_2_fit_p = poly_calc(
        data_bm3_fit_v, data_poly_2_fit)

    # calculate the 3nd order poly fit of the origional data set
    data_poly_3_fit, _ = poly_fit(data_v, data_p,n=3)
    data_poly_3_fit_p = poly_calc(
        data_bm3_fit_v, data_poly_3_fit)
    
    # calculate the 4nd order poly fit of the origional data set
    data_poly_4_fit, _ = poly_fit(data_v, data_p)
    data_poly_4_fit_p = poly_calc(
        data_bm3_fit_v, data_poly_4_fit)
   

    # plot of fitting curves

    plt.plot(data_bm3_fit_v, data_bm3_fit_p, label="BM3 fitting", lw=2)
    plt.plot(data_bm3_fit_v, data_poly_2_fit_p,
             label="Poly 2nd fitting of original data", lw=2, ls=":")
    plt.plot(data_bm3_fit_v, data_poly_3_fit_p,
             label="Poly 3rd fitting of original data", lw=2,ls=":")
    plt.plot(data_bm3_fit_v, data_poly_4_fit_p,
             label="Poly 4th fitting of original data", lw=2, ls="-")
             
    plt.scatter(data_v, data_p,
                marker=".", label="Original data", c="#BABABA")
    plt.xlabel('Volume / Atom (Å$^3$)')
    plt.ylabel('Pressure (GPa)')
    plt.xlim(6, 11)
    plt.legend(frameon=False)

if __name__ == "__main__":

    #--------------------------------------------------#
    #                 initialize data                  #
    #--------------------------------------------------#
    data = []

    # EPAW data
    epaw_v ,epaw_p = readcsv("epawmethod.csv")
    epaw = {"v":epaw_v[55:250],
            "p": epaw_p[55:250],
            "type": "ab",
            "draw":"plot",
            "ls":"solid",
            "label":"EPAW"}
    data.append(epaw)

    # Mao data, experiment
    mao_v, mao_p = readcsv("mao.csv")
    mao = {"v": mao_v,
            "p": mao_p,
            "type": "exp",
            "draw": "scatter",
            "ls": "o",
            "label": "Mao"}
    data.append(mao)

    # Dubrovinsky data, experiment
    dub_v, dub_p = readcsv("dub300.csv")
    dub = {"v": dub_v,
           "p": dub_p,
           "type": "exp",
           "draw": "scatter",
           "ls": "o",
           "label": "Dubrovinsky"}
    data.append(dub)

    # Dewaele data, experiment
    dewaele_v, dewaele_p = readcsv("dewaele.csv")
    dewaele = {"v": dewaele_v,
               "p": dewaele_p,
                "type": "exp",
                "draw": "scatter",
                "ls": "o",
                "label": "Dewaele"}
    data.append(dewaele)

    # Sakai data, experiment
    sakai_v, sakai_p = readcsv("sakai.csv")
    sakai = {"v": sakai_v,
             "p": sakai_p,
               "type": "exp",
               "draw": "scatter",
               "ls": "o",
               "label": "Sakai"}
    data.append(sakai)

    # Fei data, experiment
    fei_v, fei_p = readcsv("fei.csv")
    fei = {"v": fei_v,
           "p": fei_p,
             "type": "exp",
             "draw": "scatter",
             "ls": "o",
             "label": "Fei"}
    data.append(fei)

    # Yamazaki data, experiment, BM3 fitting
    Yamazaki_3BM_V0 = 22.15/2
    Yamazaki_3BM_K0 = 202
    Yamazaki_3BM_Kp = 4.5
    yamazaki_v, yamazaki_p = vp_bm3(Yamazaki_3BM_K0, Yamazaki_3BM_Kp, Yamazaki_3BM_V0)
    yamazaki = {"v": yamazaki_v,
                "p": yamazaki_p,
                "type": "exp",
                "draw": "plot",
                "ls": "--",
                "label": "Yamazaki"}
    data.append(yamazaki)

    # Alfe data, ab initio calculation, BM3 fitting
    Alfe_3BM_V0 = 10.20
    Alfe_3BM_K0 = 291
    Alfe_3BM_Kp = 4.4
    alfe_v, alfe_p = vp_bm3( Alfe_3BM_K0, Alfe_3BM_Kp, Alfe_3BM_V0)
    alfe = {"v": alfe_v,
            "p": alfe_p,
            "type": "ab",
            "draw": "plot",
            "ls": ":",
            "label": "Alfè"}
    data.append(alfe)
    
    # Sola data, ab initio calculation,
    sola_v, sola_p = readcsv("qmc.csv")
    sola = {"v": sola_v,
           "p": sola_p,
            "type": "ab",
           "draw": "plot",
           "ls": ":",
           "label": "Sola"}
    data.append(sola)

    # Ona data, ab initio calculation, Vinet fitting
    Ono_Vinet_V0 = 10.265
    Ono_Vinet_K0 = 285
    Ono_Vinet_Kp = 4.8
    ono_v, ono_p = vp_vinet(Ono_Vinet_K0, Ono_Vinet_Kp, Ono_Vinet_V0)
    ono = {"v": ono_v,
           "p": ono_p,
           "type": "ab",
            "draw": "plot",
            "ls": ":",
            "label": "Ono"}
    data.append(ono)

    # Sha data, ab initio calculation,  BM3 fitting
    Cohen_Vinet_V0 = (20.18/2)
    Cohen_Vinet_K0 = 296
    Cohen_Vinet_Kp = 4.4
    sha_v, sha_p = vp_vinet(Cohen_Vinet_K0, Cohen_Vinet_Kp, Cohen_Vinet_V0)
    sha = {"v": sha_v,
           "p": sha_p,
           "type": "ab",
           "draw": "plot",
           "ls": ":",
           "label": "Sha"}
    data.append(sha)


    # set up the color value
    for i, item in enumerate(data):
        item["color"] = color_list[i]

    # call the plot function to obtain the graph of raw data
    plot_raw_data(data)

    #--------------------------------------------------#
    #                   fit the data                   #
    #--------------------------------------------------#

    # seperate the experiments and ab initio calculations into two:
    experiments = [ i for i in data if i["type"]=="exp"]
    abinitio = [(i) for i in data if (i["type"] == "ab")] 
    
    vx = np.linspace(6, 11, 100)

    # plot the test of polynomial fit & bm3 fit
    plt.close()
    plt.figure(figsize=(10,3.8))
    plt.subplot(1,2,1)
    experiments_poly_fit, experiments_poly_fit_sige = plot_poly_fitting(
        experiments)
    plt.subplot(1, 2, 2)
    abinitio_poly_fit, abinitio_poly_fit_sige = plot_poly_fitting(
        abinitio)
    plt.suptitle("Fitting curve comparison")
    plt.savefig("graph/300K_diff.pdf", dpi=300)

    # plot the comparison between experiments and ab initio
    plt.close()
    plt.figure(figsize=(10, 5))
    plt.subplot(1, 2, 1)
    plt.plot(vx, poly_calc(vx, experiments_poly_fit),
             label="Experiment poly fit", c=color_exp, lw=2)
    plt.plot(vx, poly_calc(vx, abinitio_poly_fit),
             label="Ab initio poly fit", c=color_list[4], lw=2)
    plot_the_data(experiments)
    plot_the_data(abinitio)

    plt.xlabel('Volume / Atom (Å$^3$)')
    plt.ylabel('Pressure (GPa)')
    plt.xlim(6, 11.2)
    plt.ylim(-10, 400)
    plt.title("EoS fitting curve")
    plt.legend(frameon=False)
 
    plt.subplot(1, 2, 2)
    plt.plot(vx, poly_calc(vx, experiments_poly_fit),
             c=color_exp, label="Experiment poly fit", lw=2)
    plt.plot(vx, poly_calc(vx, experiments_poly_fit) +
             experiments_poly_fit_sige, ls="--", label="Experiment poly fit + $\\sigma$", c=color_exp)
    plt.plot(vx, poly_calc(vx, experiments_poly_fit) -
             experiments_poly_fit_sige, ls=":", label="Experiment poly fit - $\\sigma$", c=color_exp)
    
    plt.plot(vx, poly_calc(vx, abinitio_poly_fit),  c=color_ab, label="Ab initio poly fit",lw = 2)
    plt.plot(vx, poly_calc(vx, abinitio_poly_fit) +
             abinitio_poly_fit_sige, ls="--", label="Ab initio poly fit + $\\sigma$", c=color_ab)
    plt.plot(vx, poly_calc(vx, abinitio_poly_fit) -
             abinitio_poly_fit_sige, ls=":", label="Ab initio poly fit - $\\sigma$", c=color_ab)

    plt.xlabel('Volume / Atom (Å$^3$)')
    plt.ylabel('Pressure (GPa)')
    plt.xlim(6, 11.2)
    plt.ylim(-10, 400)
    plt.title("EoS fitting curve with uncertainty")
    plt.legend(frameon=False)
    plt.savefig("graph/300K_fitting.pdf" , dpi=300)


    #--------------------------------------------------#
    #                 why is poly 4th                  #
    #--------------------------------------------------#

    plt.close()
    plt.figure(figsize=(5, 4))
    plot_poly_fitting_veri(experiments)
    plt.xlabel('Volume / Atom (Å$^3$)')
    plt.ylabel('Pressure (GPa)')
    plt.xlim(6, 11.2)
    plt.ylim(-10, 400)
    plt.title("Test of polynomial orders")
    plt.legend(frameon=False)
    plt.savefig("graph/300K_poly4why.pdf", dpi=300)
