import math
from astroquery.nist import Nist
import astropy.units as u
import numpy as np
from scipy.special import voigt_profile
import matplotlib.pyplot as plt
import pandas as pd


table = Nist.query(180 * u.nm, 1140 * u.nm, linename="Fe I")[["Observed","Aki","Ei           Ek","gi   gk", "TP"]]
# print(table.info)
n = 0
spisok = []
for i in table:
    if not (table[n][1] != 0) or not (table[n][0] != 0):
        spisok.append(n)
    n += 1
table.remove_rows(spisok)

h = 6.63 * (10 ** (-34))
k = 1.38 * (10 ** (-23))
e = 1.6 * (10 ** (-19))
c = 3 * (10 ** 8)


def intens(g_i, g_k, E_i, E_k, A_ik, wave_lenght, T):
    E_i = E_i * e
    E_k = E_k * e
    wave_lenght = wave_lenght * (10 ** (-9))
    return (g_i / g_k) * A_ik * (h * c / (4 * math.pi * wave_lenght)) * np.exp((E_k - E_i)/(k * T))


intensivnosti = []
for i in table:
    gi, gk = int(i[3].split()[0]), int(i[3].split()[2])
    Ei, Ek = float(i[2].split()[0]), float(i[2].split()[2])
    T = int(''.join(filter(str.isdigit, i[4])))
    Aik = i[1]
    wave = i[0]
    intensivnosti.append(intens(gi, gk, Ei, Ek, Aik, wave, T))
#print(intensivnosti)

def voigt(lmd0, M, T, lmd_col):
    fig, ax = plt.subplots(figsize=(8, 8))
    x = np.linspace(lmd0 - 5, lmd0 + 5, 500)
    gamma = lmd_col
    lmd_dopl = 7.16 * (10 ** (-7)) * lmd0 * np.sqrt(T / M)
    sgm_dopl = lmd_dopl / (2 * np.sqrt(2 * np.log(2)))
    sgm_sum = []
    for i in np.arange(0.01, 1.5, 0.1):
        sgm_instr = i / (2 * np.sqrt(2 * np.log(2)))
        sgm_sum.append(np.sqrt((sgm_dopl ** 2) + (sgm_instr ** 2)))
    for sigma in sgm_sum:
        voigt = voigt_profile(x - lmd0, sigma, gamma)
        ax.plot(x, voigt, label=rf"$\sigma={sigma},\, \gamma={gamma}$",
                    ls="dashed")
    ax.legend()
    plt.show()


#voigt(table[0][0], 55.85, int(''.join(filter(str.isdigit, table[0][4]))), 0.02)


def shoooom(lmd, T):
    const = 10 ** (-15)
    hc = 1240
    T = T / 11600
    return const * np.exp((-1) * hc / (lmd * T))


def to_pandas(table):
    data = {'wavelength': [],
            'temperature': [],
            'intensity': [],
            'noise': []}
    for i in table:
        data['wavelength'].append(i[0])
        data['temperature'].append(int(''.join(filter(str.isdigit, i[4]))))
        data['noise'].append(shoooom(i[0], int(''.join(filter(str.isdigit, i[4])))))
    data['intensity'] = intensivnosti
    df = pd.DataFrame(data=data)
    return df


#print(to_pandas(table))