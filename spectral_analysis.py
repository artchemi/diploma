import cclib
import scipy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import Rbf
from scipy.constants import physical_constants


def _5Lorentzian(x, amp: list, cen: list, wid1, wid2, wid3, wid4, wid5):
    """_summary_

    Args:
        x (_type_): _description_
        amp (list): _description_
        cen (list): _description_
        wid1 (_type_): _description_
        wid2 (_type_): _description_
        wid3 (_type_): _description_
        wid4 (_type_): _description_
        wid5 (_type_): _description_

    Returns:
        _type_: _description_
    """
    return (amp[0]*wid1**2/((x-cen[0])**2+wid1**2)) +\
            (amp[1]*wid2**2/((x-cen[1])**2+wid2**2)) +\
            (amp[2]*wid3**2/((x-cen[2])**2+wid3**2)) +\
            (amp[3]*wid4**2/((x-cen[3])**2+wid4**2)) +\
            (amp[4]*wid5**2/((x-cen[4])**2+wid5**2))


def plot_spectral_methods(path: str, compound: str, solvent: str, energy_type: str) -> None:
    """График испускания/поглощения для одной молекулы 
    в одном растворителе, но с разными уровнями теории.

    Args:
        path (str): Путь до папки с результатами
        compound (str): Молекула для анализа (mc1, mc2, sp1, sp2)
        solvent (str): Растворитель (acetone, acetonitrile, chloroform, noscrf?)
        energy (str): Размерность энергии (eV, nm)
    """

    fig, ax = plt.subplots()
    palette = ['red', 'green', 'blue']  # добавить палитру

    x_range = np.linspace(200, 800, num=1000)

    wid_dict = {'wid1': 10, 'wid2': 10, 'wid3': 10, 'wid4': 10, 'wid5': 10}

    for method, color in zip(method_parse_lst, palette):
        filename = path + f'/{compound}_{solvent}/' + f'{compound}_s1_{method}.log'
        parser = cclib.io.ccopen(filename)
        data = parser.parse()

        if energy_type =='eV':
            x = data.etenergies / 10 * physical_constants['electron volt-inverse meter relationship'][0]
        else:
            x = np.power(data.etenergies, -1) * np.power(10, 7)

        y = data.etoscs

        for x_val, y_val in zip(x, y):
            ax.plot([x_val, x_val], [0, y_val], color=color)

        if method == 'b3lyp':
            label = 'CAM-B3LYP'
        elif method == 'pbe0':
            label = 'PBE0'
        else:
            label = 'M06-2X'

        ax.plot(x_range, _5Lorentzian(x_range, y, x, **wid_dict), color=color, label=label)

    ax.set(**{'xlim': (x_range.min(), x_range.max()), 'title': 'Absorption spectra', 
              'xlabel': energy_type, 'ylabel': 'Oscillator strenght'})

    plt.grid()
    plt.legend()
    plt.show()


def plot_spectral_comparison(path: str, method: str, compounds: list, solvent: str, save=False):
    fig, ax = plt.subplots()

    x_range = np.linspace(200, 800, num=1000)

    wid_dict = {'wid1': 10, 'wid2': 10, 'wid3': 10, 'wid4': 10, 'wid5': 10}

    color = ['red', 'blue']

    for compound, color in zip(compounds, color):
        filename = path + f'/{compound}_{solvent}/' + f'{compound}_s1_{method}.log'

        parser = cclib.io.ccopen(filename)
        data = parser.parse()

        x = np.power(data.etenergies, -1) * np.power(10, 7)
        y = data.etoscs

        for x_val, y_val in zip(x, y):
            ax.plot([x_val, x_val], [0, y_val], color=color)

        if method == 'b3lyp':
            label = 'CAM-B3LYP'
        elif method == 'pbe0':
            label = 'PBE0'
        else:
            label = 'M06-2X'

        ax.plot(x_range, _5Lorentzian(x_range, y, x, **wid_dict), color=color, label=compound)

        ax.set(**{'xlim': (x_range.min(), x_range.max()), 'title': f'Absorption spectra ({label}, {solvent})', 
              'xlabel': 'Wave lenght, nm', 'ylabel': 'Oscillator strenght'})

    plt.grid()
    plt.legend()

    if save == True:
        fig.savefig(f'comparison_{compounds[0]}_{compounds[1]}_{solvent}_{method}.png')

    plt.show()



def plot_spectral_solvent(path: str, compound: str, method: str) -> None:
    """_summary_

    Args:
        path (str): _description_
        method (str): _description_
    """

    fig, ax = plt.subplots()
    palette = ['red', 'green', 'blue']  # добавить палитру


def main():
    global method_parse_lst
    method_parse_lst = ['b3lyp', 'm062x', 'pbe0']
    solvent_parse_lst = ['acetone', 'acetonitrile', 'chloroform']
    path = '/home/daniil_artamonov/hpc4_kurchatov/diploma_gaussian'
    filename = path + '/mc1_acetone/' + 'mc1_s1_b3lyp.log'

    plot_spectral_comparison(path, 'b3lyp', ['mc1', 'sp1'], 'acetonitrile')

    # parser = cclib.io.ccopen(filename)
    # data = parser.parse()

    # print(data.etsecs)

    # global x_nm
    # x_nm = np.power(data.etenergies, -1) * np.power(10, 7)

    # global y_os
    # y_os = data.etoscs

    # fig, ax = plt.subplots()

    # for x_nm_val, y_os_val in zip(x_nm, y_os):
    #     ax.plot([x_nm_val, x_nm_val], [0, y_os_val], color='black')


    # x_range = np.linspace(x_nm.min() * 0.9, x_nm.max() * 1.1, num=1000)

    # wid_dict = {'wid1': 10, 'wid2': 10, 'wid3': 10, 'wid4': 10, 'wid5': 10}

    # ax.plot(x_range, _5Lorentzian(x_range, **wid_dict), '--', color='red', label='mc1_acetone')
    # ax.set(**{'xlim': (x_range.min(), x_range.max()), 'title': 'Absorption spectra'})

    # plt.show()

if __name__ == '__main__':
    main()