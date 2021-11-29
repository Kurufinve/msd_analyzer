"""
This module contains relatively simple functions for calculation of mean-squared displacements (MSD) of atoms.
The "simple" means that functions do not use sophisticated algorithms for recognition of different diffusion modes,
and can be correctly applied only if the dependence of MSD from modeling time is linear.
"""

import numpy as np
import pandas as pd
import copy
from md_format_converter_mi import structure

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


def calc_non_averaged_msd(sts_numpy, dt):
    """
    This function calculates the mean-squared displacements (msd) of atoms in structures from sts_numpy list
    with respect to the first structure in sts_numpy list.
    get: sts_numpy, dt - list of structure objects with unwrapped coordinates and time difference in ps between structures, respectively
    return: pandas dataframe, where the first column is time, 
            the second column is non-averaged msd of all atoms,
            the third column is non-averaged msd of all atoms in x direction,
            the fourth column is non-averaged msd of all atoms in y direction,
            the fifth column is non-averaged msd of all atoms in z direction,
            the sixth and other columns are non-averaged msds of atoms of i_type i

    """

    msd = {}

    msd['time'] = []
    msd['r_all'] = []
    msd['x_all'] = []
    msd['y_all'] = []
    msd['z_all'] = []

    for i_type in sts_numpy[0].type_at:
        msd[f'r_{i_type}'] = []
        msd[f'x_{i_type}'] = []
        msd[f'y_{i_type}'] = []
        msd[f'z_{i_type}'] = []

    num_sts = len(sts_numpy)
    for i in range(num_sts):
        msd_x_i = ((sts_numpy[i].r_at[:,0] - sts_numpy[0].r_at[:,0])**2).sum()/sts_numpy[i].n_at
        msd_y_i = ((sts_numpy[i].r_at[:,1] - sts_numpy[0].r_at[:,1])**2).sum()/sts_numpy[i].n_at
        msd_z_i = ((sts_numpy[i].r_at[:,2] - sts_numpy[0].r_at[:,2])**2).sum()/sts_numpy[i].n_at
        msd_r_i = msd_x_i + msd_y_i + msd_z_i

        # print(msd_r_i)
        msd['r_all'].append(msd_x_i)
        msd['x_all'].append(msd_y_i)
        msd['y_all'].append(msd_z_i)
        msd['z_all'].append(msd_r_i)

        for i_type in sts_numpy[0].type_at:
            mask = sts_numpy[i].i_type_at == i_type

            msd_x_i = ((sts_numpy[i].r_at[:,0][mask] - sts_numpy[0].r_at[:,0][mask])**2).sum()/mask.sum()
            msd_y_i = ((sts_numpy[i].r_at[:,1][mask] - sts_numpy[0].r_at[:,1][mask])**2).sum()/mask.sum()
            msd_z_i = ((sts_numpy[i].r_at[:,2][mask] - sts_numpy[0].r_at[:,2][mask])**2).sum()/mask.sum()
            msd_r_i = msd_x_i + msd_y_i + msd_z_i

            msd[f'r_{i_type}'].append(msd_x_i)
            msd[f'x_{i_type}'].append(msd_y_i)
            msd[f'y_{i_type}'].append(msd_z_i)
            msd[f'z_{i_type}'].append(msd_r_i)           

    msd_df = pd.DataFrame(msd)

    return msd_df


def calc_averaged_msd(sts_numpy, dt):

    """
    This function calculates the averaged mean-squared displacements (msd) of atoms in structures from sts_numpy list
    with respect to the first structure in sts_numpy list.
    get: sts_numpy, dt - list of structure objects with unwrapped coordinates and time difference in seconds between structure snapshots, respectively
    return: pandas dataframe, where the first column is time, 
            the second column is non-averaged msd of all atoms,
            the third column is non-averaged msd of all atoms in x direction,
            the fourth column is non-averaged msd of all atoms in y direction,
            the fifth column is non-averaged msd of all atoms in z direction,
            the sixth and other columns are non-averaged msds of atoms of i_type i

    """

    msd_av = {}

    msd_av['time'] = []
    msd_av['r_all'] = []
    msd_av['x_all'] = []
    msd_av['y_all'] = []
    msd_av['z_all'] = []

    for i_type in sts_numpy[0].type_at:
        msd_av[f'r_{i_type}'] = []
        msd_av[f'x_{i_type}'] = []
        msd_av[f'y_{i_type}'] = []
        msd_av[f'z_{i_type}'] = []

    n_write = len(sts_numpy)
    n_time_diff=n_write-1

    for i_time_diff in range(n_time_diff):

        # temporary working variables
        msd_x_i_w = 0.0
        msd_y_i_w = 0.0
        msd_z_i_w = 0.0
        msd_r_i_w = 0.0

        # dictionaries with temporary working variables
        msd_x_i_type_w = {}
        msd_y_i_type_w = {}
        msd_z_i_type_w = {}
        msd_r_i_type_w = {}

        # assigning zero values to the dictionaries with temporary working variables
        for i_type in sts_numpy[0].type_at:
            msd_x_i_type_w[f'r_{i_type}'] = 0.0
            msd_y_i_type_w[f'x_{i_type}'] = 0.0
            msd_z_i_type_w[f'y_{i_type}'] = 0.0
            msd_r_i_type_w[f'z_{i_type}'] = 0.0

        for i_start in range(n_write-i_time_diff): 
            i_end = i_start + i_time_diff

            msd_x_i_w += ((sts_numpy[i_end].r_at[:,0] - sts_numpy[i_start].r_at[:,0])**2).sum()/sts_numpy[i_end].n_at
            msd_y_i_w += ((sts_numpy[i_end].r_at[:,1] - sts_numpy[i_start].r_at[:,1])**2).sum()/sts_numpy[i_end].n_at
            msd_z_i_w += ((sts_numpy[i_end].r_at[:,2] - sts_numpy[i_start].r_at[:,2])**2).sum()/sts_numpy[i_end].n_at
            msd_r_i_w += msd_x_i_w + msd_y_i_w + msd_z_i_w

            # print(msd_r_i)


            for i_type in sts_numpy[0].type_at:
                mask = sts_numpy[i_end].i_type_at == i_type

                msd_x_i_type_w[f'r_{i_type}'] = ((sts_numpy[i_end].r_at[:,0][mask] - sts_numpy[i_start].r_at[:,0][mask])**2).sum()/mask.sum()
                msd_y_i_type_w[f'x_{i_type}'] = ((sts_numpy[i_end].r_at[:,1][mask] - sts_numpy[i_start].r_at[:,1][mask])**2).sum()/mask.sum()
                msd_z_i_type_w[f'y_{i_type}'] = ((sts_numpy[i_end].r_at[:,2][mask] - sts_numpy[i_start].r_at[:,2][mask])**2).sum()/mask.sum()
                msd_r_i_type_w[f'z_{i_type}'] = msd_x_i_type_w[f'r_{i_type}'] + msd_y_i_type_w[f'x_{i_type}'] + msd_z_i_type_w[f'y_{i_type}']


        msd_x_i_av = msd_x_i_w/float(n_write-i_time_diff)
        msd_y_i_av = msd_y_i_w/float(n_write-i_time_diff)
        msd_z_i_av = msd_z_i_w/float(n_write-i_time_diff)
        msd_r_i_av = msd_r_i_w/float(n_write-i_time_diff)

        msd_av['r_all'].append(msd_x_i_av)
        msd_av['x_all'].append(msd_y_i_av)
        msd_av['y_all'].append(msd_z_i_av)
        msd_av['z_all'].append(msd_r_i_av)

        for i_type in sts_numpy[0].type_at:
            msd_x_i_type_av = msd_x_i_type_w[f'r_{i_type}']/float(n_write-i_time_diff)
            msd_y_i_type_av = msd_y_i_type_w[f'x_{i_type}']/float(n_write-i_time_diff)
            msd_z_i_type_av = msd_z_i_type_w[f'y_{i_type}']/float(n_write-i_time_diff)
            msd_r_i_type_av = msd_r_i_type_w[f'z_{i_type}']/float(n_write-i_time_diff)        

            msd_av[f'r_{i_type}'].append(msd_x_i_type_av)
            msd_av[f'x_{i_type}'].append(msd_y_i_type_av)
            msd_av[f'y_{i_type}'].append(msd_z_i_type_av)
            msd_av[f'z_{i_type}'].append(msd_r_i_type_av)  

    msd_av_df = pd.DataFrame(msd_av)
    msd_av_df.set_index

    return msd_av_df