#!/usr/bin/python

from pyscf import gto, scf, dft, tddft, dftd3
from pyscf.solvent import pcm
from pyscf.geomopt.geometric_solver import optimize
import numpy as np
from scipy.constants import physical_constants
import time


mol = gto.M('''
C        5.065347000     -1.663703000      0.250253000
C        4.234356000     -2.527399000     -0.453955000
C        2.975294000     -2.122210000     -0.884030000
C        2.600293000     -0.830371000     -0.573984000
C        3.411382000      0.045878000      0.131947000
C        4.660379000     -0.364717000      0.548574000
H        6.040291000     -2.007537000      0.574038000
H        4.568316000     -3.534848000     -0.670433000
H        5.314677000      0.301671000      1.098791000
H        2.325585000     -2.798353000     -1.424423000
C        2.689129000      1.355591000      0.282206000
C        1.365432000      1.041070000     -0.381340000
N        1.395552000     -0.170280000     -0.893661000
C        2.473558000      1.746263000      1.745827000
H        1.928201000      0.966319000      2.275694000
H        1.902006000      2.674169000      1.806848000
H        3.441651000      1.902832000      2.224487000
C        3.411789000      2.477520000     -0.476164000
H        3.539598000      2.227384000     -1.530472000
H        4.397921000      2.625213000     -0.033474000
H        2.861217000      3.416219000     -0.397323000
C        0.317180000      1.986508000     -0.516588000
C       -1.821098000      0.722619000      0.040636000
C       -1.307150000     -0.299294000      0.950949000
O       -0.128277000     -0.318591000      1.354371000
C       -2.261610000     -1.292023000      1.386447000
C       -3.561543000     -1.286230000      0.994420000
C       -4.029770000     -0.259603000      0.148231000
C       -3.173625000      0.728278000     -0.293286000
C        0.381210000     -0.781101000     -1.715593000
H       -0.245410000     -0.002876000     -2.142734000
H       -0.235750000     -1.454912000     -1.119261000
H        0.868279000     -1.341232000     -2.512241000
C       -1.029534000      1.812264000     -0.442392000
H        0.665196000      3.002302000     -0.679014000
H       -1.611788000      2.682834000     -0.735903000
H       -3.565333000      1.512464000     -0.929545000
H       -4.252151000     -2.049550000      1.327953000
H       -1.888036000     -2.063941000      2.049443000
N       -5.396991000     -0.234076000     -0.244867000
O       -5.789717000      0.670614000     -0.973260000
O       -6.138319000     -1.122486000      0.160680000''', basis='631g', verbose=3)


def timecalc(func):
    def wrapper(*args, **kwargs):
        print(f'Running {func.__name__}...')
        t_start = time.time()
        result = func(*args, **kwargs)
        t_end = time.time()
        print(f'{func.__name__} took {t_end - t_start} seconds')
        print(f'{func.__name__} stopped')
        return result

    return wrapper


@timecalc
def ground_state(mol, xc: str):
    global cm
    cm = pcm.PCM(mol)
    cm.eps = 20.493  # acetone dielectric constant
    cm.method = 'C-PCM'

    # дописать дисперсионные поправки
    mf = dftd3.dftd3(dft.RKS(mol, xc=xc).PCM(cm).run())

    return mf


@timecalc
def gs_optimization(mf):
    mol_eq = mf.Gradients().optimizer(solver='geomeTRIC').kernel()
    print(mol_eq.atom_coords())

    return mol_eq

@timecalc
def main():
    ha_2_ev = 1/physical_constants["electron volt-hartree relationship"][0]

    #xc = 'PBE0'
    xc = 'lda'
    
    mf = ground_state(mol, xc)
    mol_eq = gs_optimization(mf)

    mytd = tddft.TDDFT(mol_eq).PCM(cm)
    mytd.nstates = 5
    mytd.max_space = 100
    mytd.max_cycle = 200

    mytd.kernel()
    mytd.analyze()

    print(mytd.transition_dipole())



if __name__ == '__main__':
    main()