"""
This module contains relatively simple functions needed for calculation of mean-squared displacements (MSD) of atoms from series of time snapshots.
The "simple" means that functions do not use sophisticated algorithms for recognition of different diffusion modes,
and can be correctly applied only if the dependence of MSD from modeling time is linear.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import copy
from md_format_converter_mi import structure

# atom data from phonopy (https://github.com/phonopy/phonopy/blob/develop/phonopy/structure/atoms.py)
atom_data = [
    [0, "X", "X", None],  # 0
    [1, "H", "Hydrogen", 1.00794],  # 1
    [2, "He", "Helium", 4.002602],  # 2
    [3, "Li", "Lithium", 6.941],  # 3
    [4, "Be", "Beryllium", 9.012182],  # 4
    [5, "B", "Boron", 10.811],  # 5
    [6, "C", "Carbon", 12.0107],  # 6
    [7, "N", "Nitrogen", 14.0067],  # 7
    [8, "O", "Oxygen", 15.9994],  # 8
    [9, "F", "Fluorine", 18.9984032],  # 9
    [10, "Ne", "Neon", 20.1797],  # 10
    [11, "Na", "Sodium", 22.98976928],  # 11
    [12, "Mg", "Magnesium", 24.3050],  # 12
    [13, "Al", "Aluminium", 26.9815386],  # 13
    [14, "Si", "Silicon", 28.0855],  # 14
    [15, "P", "Phosphorus", 30.973762],  # 15
    [16, "S", "Sulfur", 32.065],  # 16
    [17, "Cl", "Chlorine", 35.453],  # 17
    [18, "Ar", "Argon", 39.948],  # 18
    [19, "K", "Potassium", 39.0983],  # 19
    [20, "Ca", "Calcium", 40.078],  # 20
    [21, "Sc", "Scandium", 44.955912],  # 21
    [22, "Ti", "Titanium", 47.867],  # 22
    [23, "V", "Vanadium", 50.9415],  # 23
    [24, "Cr", "Chromium", 51.9961],  # 24
    [25, "Mn", "Manganese", 54.938045],  # 25
    [26, "Fe", "Iron", 55.845],  # 26
    [27, "Co", "Cobalt", 58.933195],  # 27
    [28, "Ni", "Nickel", 58.6934],  # 28
    [29, "Cu", "Copper", 63.546],  # 29
    [30, "Zn", "Zinc", 65.38],  # 30
    [31, "Ga", "Gallium", 69.723],  # 31
    [32, "Ge", "Germanium", 72.64],  # 32
    [33, "As", "Arsenic", 74.92160],  # 33
    [34, "Se", "Selenium", 78.96],  # 34
    [35, "Br", "Bromine", 79.904],  # 35
    [36, "Kr", "Krypton", 83.798],  # 36
    [37, "Rb", "Rubidium", 85.4678],  # 37
    [38, "Sr", "Strontium", 87.62],  # 38
    [39, "Y", "Yttrium", 88.90585],  # 39
    [40, "Zr", "Zirconium", 91.224],  # 40
    [41, "Nb", "Niobium", 92.90638],  # 41
    [42, "Mo", "Molybdenum", 95.96],  # 42
    [43, "Tc", "Technetium", 98],  # 43 (mass is from wikipedia)
    [44, "Ru", "Ruthenium", 101.07],  # 44
    [45, "Rh", "Rhodium", 102.90550],  # 45
    [46, "Pd", "Palladium", 106.42],  # 46
    [47, "Ag", "Silver", 107.8682],  # 47
    [48, "Cd", "Cadmium", 112.411],  # 48
    [49, "In", "Indium", 114.818],  # 49
    [50, "Sn", "Tin", 118.710],  # 50
    [51, "Sb", "Antimony", 121.760],  # 51
    [52, "Te", "Tellurium", 127.60],  # 52
    [53, "I", "Iodine", 126.90447],  # 53
    [54, "Xe", "Xenon", 131.293],  # 54
    [55, "Cs", "Caesium", 132.9054519],  # 55
    [56, "Ba", "Barium", 137.327],  # 56
    [57, "La", "Lanthanum", 138.90547],  # 57
    [58, "Ce", "Cerium", 140.116],  # 58
    [59, "Pr", "Praseodymium", 140.90765],  # 59
    [60, "Nd", "Neodymium", 144.242],  # 60
    [61, "Pm", "Promethium", 145],  # 61 (mass is from wikipedia)
    [62, "Sm", "Samarium", 150.36],  # 62
    [63, "Eu", "Europium", 151.964],  # 63
    [64, "Gd", "Gadolinium", 157.25],  # 64
    [65, "Tb", "Terbium", 158.92535],  # 65
    [66, "Dy", "Dysprosium", 162.500],  # 66
    [67, "Ho", "Holmium", 164.93032],  # 67
    [68, "Er", "Erbium", 167.259],  # 68
    [69, "Tm", "Thulium", 168.93421],  # 69
    [70, "Yb", "Ytterbium", 173.054],  # 70
    [71, "Lu", "Lutetium", 174.9668],  # 71
    [72, "Hf", "Hafnium", 178.49],  # 72
    [73, "Ta", "Tantalum", 180.94788],  # 73
    [74, "W", "Tungsten", 183.84],  # 74
    [75, "Re", "Rhenium", 186.207],  # 75
    [76, "Os", "Osmium", 190.23],  # 76
    [77, "Ir", "Iridium", 192.217],  # 77
    [78, "Pt", "Platinum", 195.084],  # 78
    [79, "Au", "Gold", 196.966569],  # 79
    [80, "Hg", "Mercury", 200.59],  # 80
    [81, "Tl", "Thallium", 204.3833],  # 81
    [82, "Pb", "Lead", 207.2],  # 82
    [83, "Bi", "Bismuth", 208.98040],  # 83
    [84, "Po", "Polonium", None],  # 84
    [85, "At", "Astatine", None],  # 85
    [86, "Rn", "Radon", None],  # 86
    [87, "Fr", "Francium", None],  # 87
    [88, "Ra", "Radium", None],  # 88
    [89, "Ac", "Actinium", 227],  # 89 (mass is from wikipedia)
    [90, "Th", "Thorium", 232.03806],  # 90
    [91, "Pa", "Protactinium", 231.03588],  # 91
    [92, "U", "Uranium", 238.02891],  # 92
    [93, "Np", "Neptunium", 237],  # 93 (mass is from wikipedia)
    [94, "Pu", "Plutonium", None],  # 94
    [95, "Am", "Americium", None],  # 95
    [96, "Cm", "Curium", None],  # 96
    [97, "Bk", "Berkelium", None],  # 97
    [98, "Cf", "Californium", None],  # 98
    [99, "Es", "Einsteinium", None],  # 99
    [100, "Fm", "Fermium", None],  # 100
    [101, "Md", "Mendelevium", None],  # 101
    [102, "No", "Nobelium", None],  # 102
    [103, "Lr", "Lawrencium", None],  # 103
    [104, "Rf", "Rutherfordium", None],  # 104
    [105, "Db", "Dubnium", None],  # 105
    [106, "Sg", "Seaborgium", None],  # 106
    [107, "Bh", "Bohrium", None],  # 107
    [108, "Hs", "Hassium", None],  # 108
    [109, "Mt", "Meitnerium", None],  # 109
    [110, "Ds", "Darmstadtium", None],  # 110
    [111, "Rg", "Roentgenium", None],  # 111
    [112, "Cn", "Copernicium", None],  # 112
    [113, "Uut", "Ununtrium", None],  # 113
    [114, "Uuq", "Ununquadium", None],  # 114
    [115, "Uup", "Ununpentium", None],  # 115
    [116, "Uuh", "Ununhexium", None],  # 116
    [117, "Uus", "Ununseptium", None],  # 117
    [118, "Uuo", "Ununoctium", None],  # 118
]

def convert_structure_to_numpy(st):
    """ 
    This funciton converts structure object with properties described 
    by standard python types to analogous structure object but with properties
    described by numpy arrays.
    
    """
    st_numpy = structure()
    st_numpy.n_at = st.n_at
    st_numpy.n_mark_at = st.n_mark_at
    st_numpy.mark_at = np.array(st.mark_at)
    st_numpy.i_at = np.array(st.i_at)
    st_numpy.r_at = np.array(st.r_at)
    st_numpy.f_at = np.array(st.f_at)
    st_numpy.v_at = np.array(st.v_at)
    st_numpy.sizex = st.sizex
    st_numpy.sizey = st.sizey
    st_numpy.sizez = st.sizez
    st_numpy.a_lattice3 = np.array(st.a_lattice3)
    st_numpy.n_type_at = st.n_type_at
    st_numpy.i_type_at = np.array(st.i_type_at)
    st_numpy.i_mass_at = np.array(st.i_mass_at)
    st_numpy.type_at = np.array(st.type_at)
    st_numpy.mass_at = np.array(st.mass_at)

    return st_numpy

def get_cm_corrected_structure(st_numpy):
    """ Function for correction of atoms coordinates to the center of mass"""
    cm_x = (st_numpy.i_mass_at*st_numpy.r_at[:,0]).sum()/st_numpy.i_mass_at.sum()
    cm_y = (st_numpy.i_mass_at*st_numpy.r_at[:,1]).sum()/st_numpy.i_mass_at.sum()
    cm_z = (st_numpy.i_mass_at*st_numpy.r_at[:,2]).sum()/st_numpy.i_mass_at.sum()

    # print(f"cm_x = {cm_x}")
    # print(f"cm_y = {cm_y}")
    # print(f"cm_z = {cm_z}")

    st_numpy_cm = copy.deepcopy(st_numpy)
    st_numpy_cm.r_at[:,0] = st_numpy.r_at[:,0] - cm_x
    st_numpy_cm.r_at[:,1] = st_numpy.r_at[:,1] - cm_y
    st_numpy_cm.r_at[:,2] = st_numpy.r_at[:,2] - cm_z

    return st_numpy_cm

def get_unwrapped_structures(sts_numpy):
    """ 
    Function for unwrapping the atomic coodinates 
    for correct calculation of atomic displacements
    get: sts_numpy - list of structure objects with attributes represented as numpy arrays
    return: sts_numpy_new - list of structure objects with unwrapped coordinates
    """

    sts_numpy_unwrapped = copy.deepcopy(sts_numpy) # deep copy of list with structure objects
    num_sts = len(sts_numpy) # number of structures
    # sts_numpy_new[0] = copy.deepcopy(sts_numpy[0]) # copying the first structure in list
    for i in range(1,num_sts):
        # st_numpy_i = sts_numpy[i]
        # st_numpy_i
        shift = np.zeros((sts_numpy[i].n_at,3)) # initializing array of atom shift with zeros

        sts_numpy_unwrapped[i] = copy.deepcopy(sts_numpy[i]) # copying the input structure to the output structure

        r_at_wrapped = sts_numpy[i].r_at  # initial (wrapped within the periodic boundary conditions) coordinates of atoms

        r_at_unwrapped = copy.deepcopy(r_at_wrapped)

        dx_arr = sts_numpy[i].r_at[:,0] - sts_numpy[i-1].r_at[:,0] # array of diferrences of coordinates in x direction
        dy_arr = sts_numpy[i].r_at[:,1] - sts_numpy[i-1].r_at[:,1] # array of diferrences of coordinates in y direction
        dz_arr = sts_numpy[i].r_at[:,2] - sts_numpy[i-1].r_at[:,2] # array of diferrences of coordinates in z direction
        # for iat, dx, dy, dz in zip(range(sts_numpy[i].n_at), dx_arr, dy_arr, dz_arr):                              
        for iat in range(sts_numpy[i].n_at):                              
            if (dx_arr[iat] >  sts_numpy[i-1].sizex/2): shift[iat,0] = shift[iat,0] - 0.5*(sts_numpy[i].sizex + sts_numpy[i-1].sizex)
            if (dx_arr[iat] < -sts_numpy[i-1].sizex/2): shift[iat,0] = shift[iat,0] + 0.5*(sts_numpy[i].sizex + sts_numpy[i-1].sizex)
            if (dy_arr[iat] >  sts_numpy[i-1].sizey/2): shift[iat,1] = shift[iat,1] - 0.5*(sts_numpy[i].sizey + sts_numpy[i-1].sizey)
            if (dy_arr[iat] < -sts_numpy[i-1].sizey/2): shift[iat,1] = shift[iat,1] + 0.5*(sts_numpy[i].sizey + sts_numpy[i-1].sizey)
            if (dz_arr[iat] >  sts_numpy[i-1].sizez/2): shift[iat,2] = shift[iat,2] - 0.5*(sts_numpy[i].sizez + sts_numpy[i-1].sizez)
            if (dz_arr[iat] < -sts_numpy[i-1].sizez/2): shift[iat,2] = shift[iat,2] + 0.5*(sts_numpy[i].sizez + sts_numpy[i-1].sizez) 

        r_at_unwrapped[:,0] = sts_numpy_unwrapped[i-1].r_at[:,0] + (sts_numpy[i].r_at[:,0] + shift[:,0] - sts_numpy[i-1].r_at[:,0])     
        r_at_unwrapped[:,1] = sts_numpy_unwrapped[i-1].r_at[:,1] + (sts_numpy[i].r_at[:,1] + shift[:,1] - sts_numpy[i-1].r_at[:,1])     
        r_at_unwrapped[:,2] = sts_numpy_unwrapped[i-1].r_at[:,2] + (sts_numpy[i].r_at[:,2] + shift[:,2] - sts_numpy[i-1].r_at[:,2])     

        sts_numpy_unwrapped[i].r_at = r_at_unwrapped

    return sts_numpy_unwrapped


def calc_non_averaged_msd(sts_numpy, dt=7.5**(-11)):
    """
    This function calculates the mean-squared displacements (msd) of atoms in structures from sts_numpy list
    with respect to the first structure in sts_numpy list.
    get: sts_numpy, dt - list of structure objects with unwrapped coordinates and time difference in seconds between structures, respectively
    return: pandas dataframe, where the first column is time, 
            the second column is non-averaged msd of all atoms,
            the third column is non-averaged msd of all atoms in x direction,
            the fourth column is non-averaged msd of all atoms in y direction,
            the fifth column is non-averaged msd of all atoms in z direction,
            the sixth and other columns are non-averaged msds of atoms of i_mass i

    """

    msd = {}

    msd['time, s'] = []
    msd['msd_x_all, m^2'] = []
    msd['msd_y_all, m^2'] = []
    msd['msd_z_all, m^2'] = []
    msd['msd_all, m^2'] = []

    for i_mass in sts_numpy[0].mass_at:
        element = mass2element(i_mass)
        msd[f'msd_x_{element}, m^2'] = []
        msd[f'msd_y_{element}, m^2'] = []
        msd[f'msd_z_{element}, m^2'] = []
        msd[f'msd_{element}, m^2'] = []

    num_sts = len(sts_numpy)
    for i in range(num_sts):
        msd['time, s'].append(i*dt)

        # calculating MSD for all atoms in structure
        msd_x_i = ((sts_numpy[i].r_at[:,0] - sts_numpy[0].r_at[:,0])**2).sum()/sts_numpy[i].n_at
        msd_y_i = ((sts_numpy[i].r_at[:,1] - sts_numpy[0].r_at[:,1])**2).sum()/sts_numpy[i].n_at
        msd_z_i = ((sts_numpy[i].r_at[:,2] - sts_numpy[0].r_at[:,2])**2).sum()/sts_numpy[i].n_at
        msd_r_i = msd_x_i + msd_y_i + msd_z_i

        msd['msd_x_all, m^2'].append(msd_x_i/10**20)
        msd['msd_y_all, m^2'].append(msd_y_i/10**20)
        msd['msd_z_all, m^2'].append(msd_z_i/10**20)
        msd['msd_all, m^2'].append(msd_r_i/10**20)

        # calculating MSD for each type of atom in structure
        for i_mass in sts_numpy[0].mass_at:
            element = mass2element(i_mass)

            mask = sts_numpy[i].i_mass_at == i_mass

            msd_x_i = ((sts_numpy[i].r_at[:,0][mask] - sts_numpy[0].r_at[:,0][mask])**2).sum()/mask.sum()
            msd_y_i = ((sts_numpy[i].r_at[:,1][mask] - sts_numpy[0].r_at[:,1][mask])**2).sum()/mask.sum()
            msd_z_i = ((sts_numpy[i].r_at[:,2][mask] - sts_numpy[0].r_at[:,2][mask])**2).sum()/mask.sum()
            msd_r_i = msd_x_i + msd_y_i + msd_z_i

            msd[f'msd_x_{element}, m^2'].append(msd_x_i/10**20)
            msd[f'msd_y_{element}, m^2'].append(msd_y_i/10**20)
            msd[f'msd_z_{element}, m^2'].append(msd_z_i/10**20)           
            msd[f'msd_{element}, m^2'].append(msd_r_i/10**20)

    msd_df = pd.DataFrame(msd)
    # msd_df.set_index('time', inplace = True)
    msd_df.set_index('time, s')

    return msd_df

def mass2element(mass):
    """
    This function maps the atomic mass to the element name 
    """

    for row in atom_data:
        if row[3] != None:
            if np.round(row[3],2) == np.round(mass,2):
                element = row[1]
                break

    return element

def calc_averaged_msd(sts_numpy, dt=7.5E-11):

    """
    This function calculates the averaged mean-squared displacements (msd) of atoms in structures from sts_numpy list
    with respect to the first structure in sts_numpy list.
    get: sts_numpy, dt - list of structure objects with unwrapped coordinates and time difference in seconds between structure snapshots, respectively
    return: pandas dataframe, where the first column is time, 
            the second column is non-averaged msd of all atoms,
            the third column is non-averaged msd of all atoms in x direction,
            the fourth column is non-averaged msd of all atoms in y direction,
            the fifth column is non-averaged msd of all atoms in z direction,
            the sixth and other columns are non-averaged msds of atoms of i_mass i
            the MSD are outputed in m^2/atom units

    """

    msd_av = {}

    msd_av['time, s'] = []
    msd_av['msd_all, m^2'] = []
    msd_av['msd_x_all, m^2'] = []
    msd_av['msd_y_all, m^2'] = []
    msd_av['msd_z_all, m^2'] = []

    for i_mass in sts_numpy[0].mass_at:
        element = mass2element(i_mass)
        msd_av[f'msd_{element}, m^2'] = []
        msd_av[f'msd_x_{element}, m^2'] = []
        msd_av[f'msd_y_{element}, m^2'] = []
        msd_av[f'msd_z_{element}, m^2'] = []

    n_write = len(sts_numpy)
    n_time_diff=n_write

    for i_time_diff in range(n_time_diff):
        msd_av['time, s'].append(i_time_diff*dt)

        # temporary working variables
        msd_x_i_w = 0.0
        msd_y_i_w = 0.0
        msd_z_i_w = 0.0
        msd_r_i_w = 0.0

        # dictionaries with temporary working variables
        msd_x_i_mass_w = {}
        msd_y_i_mass_w = {}
        msd_z_i_mass_w = {}
        msd_r_i_mass_w = {}

        # assigning zero values to the dictionaries with temporary working variables
        for i_mass in sts_numpy[0].mass_at:
            element = mass2element(i_mass)
            msd_x_i_mass_w[f'msd_x_{element}'] = 0.0
            msd_y_i_mass_w[f'msd_y_{element}'] = 0.0
            msd_z_i_mass_w[f'msd_z_{element}'] = 0.0
            msd_r_i_mass_w[f'msd_{element}'] = 0.0

        for i_start in range(n_write-i_time_diff): 
            i_end = i_start + i_time_diff

            # gaining the sum of MSD between the same times differences for all atoms
            msd_x_i_ww = ((sts_numpy[i_end].r_at[:,0] - sts_numpy[i_start].r_at[:,0])**2).sum()/sts_numpy[i_end].n_at
            msd_y_i_ww = ((sts_numpy[i_end].r_at[:,1] - sts_numpy[i_start].r_at[:,1])**2).sum()/sts_numpy[i_end].n_at
            msd_z_i_ww = ((sts_numpy[i_end].r_at[:,2] - sts_numpy[i_start].r_at[:,2])**2).sum()/sts_numpy[i_end].n_at

            msd_x_i_w += msd_x_i_ww
            msd_y_i_w += msd_y_i_ww
            msd_z_i_w += msd_z_i_ww
            msd_r_i_w += msd_x_i_ww + msd_y_i_ww + msd_z_i_ww

            # gaining the sum of differences between the same times differences for each type of atom
            for i_mass in sts_numpy[0].mass_at:
                element = mass2element(i_mass)
                mask = sts_numpy[i_end].i_mass_at == i_mass

                msd_x_i_mass_ww = ((sts_numpy[i_end].r_at[:,0][mask] - sts_numpy[i_start].r_at[:,0][mask])**2).sum()/mask.sum()
                msd_y_i_mass_ww = ((sts_numpy[i_end].r_at[:,1][mask] - sts_numpy[i_start].r_at[:,1][mask])**2).sum()/mask.sum()
                msd_z_i_mass_ww = ((sts_numpy[i_end].r_at[:,2][mask] - sts_numpy[i_start].r_at[:,2][mask])**2).sum()/mask.sum()

                msd_x_i_mass_w[f'msd_x_{element}'] += msd_x_i_mass_ww
                msd_y_i_mass_w[f'msd_y_{element}'] += msd_y_i_mass_ww
                msd_z_i_mass_w[f'msd_z_{element}'] += msd_z_i_mass_ww
                msd_r_i_mass_w[f'msd_{element}'] += msd_x_i_mass_ww + msd_y_i_mass_ww + msd_z_i_mass_ww

        # calculating MSD for all atoms in structure averaged over the same time differences
        msd_x_i_av = msd_x_i_w/float(n_write-i_time_diff)
        msd_y_i_av = msd_y_i_w/float(n_write-i_time_diff)
        msd_z_i_av = msd_z_i_w/float(n_write-i_time_diff)
        msd_r_i_av = msd_r_i_w/float(n_write-i_time_diff)

        msd_av['msd_all, m^2'].append(msd_r_i_av/10**20)
        msd_av['msd_x_all, m^2'].append(msd_x_i_av/10**20)
        msd_av['msd_y_all, m^2'].append(msd_y_i_av/10**20)
        msd_av['msd_z_all, m^2'].append(msd_z_i_av/10**20)

        # calculating MSD for each type of atom in structure averaged over the same time differences 
        for i_mass in sts_numpy[0].mass_at:
            element = mass2element(i_mass)
            msd_x_i_mass_av = msd_x_i_mass_w[f'msd_x_{element}']/float(n_write-i_time_diff)
            msd_y_i_mass_av = msd_y_i_mass_w[f'msd_y_{element}']/float(n_write-i_time_diff)
            msd_z_i_mass_av = msd_z_i_mass_w[f'msd_z_{element}']/float(n_write-i_time_diff)
            msd_r_i_mass_av = msd_r_i_mass_w[f'msd_{element}']/float(n_write-i_time_diff)        

            msd_av[f'msd_{element}, m^2'].append(msd_r_i_mass_av/10**20)  
            msd_av[f'msd_x_{element}, m^2'].append(msd_x_i_mass_av/10**20)
            msd_av[f'msd_y_{element}, m^2'].append(msd_y_i_mass_av/10**20)
            msd_av[f'msd_z_{element}, m^2'].append(msd_z_i_mass_av/10**20)

    msd_av_df = pd.DataFrame(msd_av)
    # msd_av_df.set_index('time', inplace = True)
    msd_av_df.set_index('time, s')

    return msd_av_df

def fit_x_y_linear(x,y):
    """
    This function linearly fits the dependence y on x.
    get: x, y - numpy arrays with x and y values
    return: a, b - coefficients of the expression y = a*x + b
    """

    A = np.vstack([x, np.ones(len(x))]).T

    a, b = np.linalg.lstsq(A, y, rcond=None)[0]

    return a, b

def fit_msd_linear(msd_df):
    """
    This function linearly fits the time dependence of MSD of atoms. 
    get: msd_df - pandas dataframe with MSD
    return: coeffs - list of tuples with coeffs a and b of linear fit (msd = a*t + b), 
                     each tuple contain coeffs from fit of msd in each column in msd_df
    """

    coeffs = []
    t = msd_df['time, s']
    for col in msd_df:
        if col == 'Unnamed: 0' or col == 'time, s':
            pass
        else:
            msd = msd_df[col]
            a, b = fit_x_y_linear(t, msd)
            coeffs.append((a,b))

    return coeffs

def fit_and_plot_msd_linear(msd_df):
    """ 
    This function fits and plots linear dependence of msd from the modeling time. 
    """

    coeffs = fit_msd_linear(msd_df)
    t = msd_df['time, s']

    # print(coeffs)
    for coef, col in zip(coeffs,msd_df.drop(columns=['Unnamed: 0','time, s'])):
        a = coef[0]
        b = coef[1]
        figname = col.split(',')[0]+'.png'
        plt.title(f'{col}, D = {a/6:.2}')
        plt.plot(t,msd_df[col],linestyle='',marker='D',label='data')
        plt.plot(t, a*t + b, label=f'fit: msd = {a:.2}*t + {b:.2}')
        plt.xlabel('time, $s$')
        plt.ylabel('MSD, ${m^2}$')
        plt.legend()
        plt.savefig(figname, format = 'png')
        plt.clf()
        plt.cla()

