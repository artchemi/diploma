import cclib
import scipy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import Rbf
from scipy.constants import physical_constants
from energydiagram import ED
from matplotlib.gridspec import GridSpec


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


def plot_spectral_methods(path: str, methods: list, compound: str, solvent: str, energy_type: str) -> None:
    """График испускания/поглощения для одной молекулы 
    в одном растворителе, но с разными уровнями теории.

    Args:
        methods (list): Список из уровней теории
        path (str): Путь до папки с результатами
        compound (str): Молекула для анализа (mc1, mc2, sp1, sp2)
        solvent (str): Растворитель (acetone, acetonitrile, chloroform, noscrf?)
        energy (str): Размерность энергии (eV, nm)
    """

    fig, ax = plt.subplots()
    palette = ['red', 'green', 'blue']  # добавить палитру

    x_range = np.linspace(200, 800, num=1000)

    wid_dict = {'wid1': 10, 'wid2': 10, 'wid3': 10, 'wid4': 10, 'wid5': 10}

    for method, color in zip(methods, palette):
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

    palette = ['red', 'blue']

    for compound, color in zip(compounds, palette):
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


def plot_spectral_solvent(path: str, solvents: list, compound: str, method: str, save=False, ret_fig=False):
    """_summary_

    Args:
        path (str): _description_
        method (str): _description_
    """

    fig, ax = plt.subplots()
    palette = ['red', 'green', 'blue']  # добавить палитру

    x_range = np.linspace(200, 800, num=1000)

    wid_dict = {'wid1': 10, 'wid2': 10, 'wid3': 10, 'wid4': 10, 'wid5': 10}

    for solvent, color in zip(solvents, palette):
        filename = path + f'/{compound}_{solvent}/' + f'{compound}_s1_{method}.log'

        parser = cclib.io.ccopen(filename)
        data = parser.parse()

        x = np.power(data.etenergies, -1) * np.power(10, 7)
        y = data.etoscs

        for x_val, y_val in zip(x, y):
            ax.plot([x_val, x_val], [0, y_val], color=color)

        ax.plot(x_range, _5Lorentzian(x_range, y, x, **wid_dict), color=color, label=solvent)

        if method == 'b3lyp':
            label = 'CAM-B3LYP'
        elif method == 'pbe0':
            label = 'PBE0'
        else:
            label = 'M06-2X'

        ax.set(**{'xlim': (x_range.min(), x_range.max()), 'title': f'Absorption spectra ({label}, {compound})', 
              'xlabel': 'Wave lenght, nm', 'ylabel': 'Oscillator strenght'})
        
    plt.grid()
    plt.legend()

    if save == True:
        fig.savefig(f'solvents_{compound}_{method}.png')

    # Для построения спектров вместе с диаграммой МО
    if ret_fig == True:
        plt.close(fig)
        return np.array([x_range, _5Lorentzian(x_range, y, x, **wid_dict)]), data.moenergies, data.etsecs, x, y

    plt.show()


def plot_mo() -> None:

    #fig, axs = plt.subplots(1, 2, figsize=(9, 3), sharey=True)
    # gs = GridSpec(1, 2, width_ratios=None, wspace=0)
    # ax1 = fig.add_subplot(gs[0, 0])
    # ax2 = fig.add_subplot(gs[0,1])

    fig = plt.figure(figsize=(15, 10), layout='constrained')
    axs = fig.subplot_mosaic([["spectra", "mo_diagram"]])

    data, mo_data, transition_data, x, y = plot_spectral_solvent(path, solvent_parse_lst, 'mc1', 'm062x', ret_fig=True)
    

    count_state = 1
    for x_val, y_val in zip(x, y):
        axs['spectra'].plot([x_val, x_val], [0, y_val], color='black')
        axs['spectra'].text(x_val, 0, f'{count_state}', 
                            bbox=dict(facecolor='white', alpha=1, edgecolor='white'))
        count_state += 1

    axs['spectra'].plot(*data, color='black')
    axs['spectra'].set(**{'xlabel': 'Wave lenght / $nm$', 'ylabel': 'Oscillator strenght', 
                          'title': 'Absorption spectra'})
    
    trans_1 = transition_data[0][0]
    print(transition_data[1])

    diagram = ED()
    count_trans = 0
    for etrans in transition_data:
        new_col = True

        for trans in etrans:
            if new_col:
                diagram.add_level(np.round(mo_data[0][trans[0][0]], 4), 
                                  left_text=f'MO #{trans[0][0]}', top_text='')
            else:
                diagram.add_level(np.round(mo_data[0][trans[0][0]], 4), left_text=f'MO #{trans[0][0]}', position='l')

            start_id = count_trans
            count_trans += 1

            diagram.add_level(np.round(mo_data[0][trans[1][0]], 4), left_text=f'MO #{trans[1][0]}', position='l')

            end_id = count_trans
            count_trans += 1

            diagram.add_link(start_id, end_id, line_order=2)
            # diagram.add_arrow(start_id, end_id, position='center', text=None)
            # diagram.add_arrow(np.round(mo_data[0][trans[0][0]], 4), 
            #                   np.round(mo_data[0][trans[1][0]], 4), position='center')
            new_col =False
    
    # axs['spectra'].set(**{'xlim':})
    
    # diagram.offset *= 1.5
    
    diagram.plot(ylabel="Energy / $eV$", 
                 show_IDs=False, ax=axs['mo_diagram'])
    diagram.fig.set_figheight(20)
    diagram.ax.axes.get_xaxis().set_visible(True)
    
    diagram.ax.spines['bottom'].set_visible(True)
    diagram.ax.set_ylim(-10, 2)
    diagram.ax.plot([2, 2], [-5, 2])
    diagram.fig.show()
    plt.show()


def main():
    global method_parse_lst, solvent_parse_lst, path
    method_parse_lst = ['b3lyp', 'm062x', 'pbe0']
    solvent_parse_lst = ['acetone', 'acetonitrile', 'chloroform']
    path = '/home/daniil_artamonov/hpc4_kurchatov/diploma_gaussian'
    filename = path + '/mc1_acetone/' + 'mc1_s1_b3lyp.log'

    # plot_spectral_solvent(path, solvent_parse_lst, 'sp1', 'm062x')
    plot_mo()


if __name__ == '__main__':
    main()