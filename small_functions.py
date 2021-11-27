import numpy as np
import pandas as pd
import copy
from md_format_converter_mi import structure

# class numpy_structure(structure):
#     def __init__(self):
#         super().__init__()

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
    return: msd_arr - num_sts*(5+4*n_type_at) array, where the first column is time, 
            the second column is non-averaged msd of all atoms,
            the third column is non-averaged msd of all atoms in x direction,
            the fourth column is non-averaged msd of all atoms in y direction,
            the fifth column is non-averaged msd of all atoms in z direction,
            the sixth and other columns are non-averaged msds of atoms of i_type i

    return pandas dataframe        
    """

    msd = {}

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


    # n_time_diff=n_write-1

	# do i_time_diff=1,n_time_diff

	# ww=0.0D0
	# ww_sort=0.0D0
	# do i_start=1,n_write-i_time_diff
	# i_end=i_start+i_time_diff

	#     w=0.0D0
	#     w_sort=0.0D0
	#     n_at_sort=0
	#     do i_at=1,n_at
	#     d =     (r_at_set(1,i_at,i_end)-r_at_set(1,i_at,i_start))**2 + 
    #  :		(r_at_set(2,i_at,i_end)-r_at_set(2,i_at,i_start))**2 + 
    #  :		(r_at_set(3,i_at,i_end)-r_at_set(3,i_at,i_start))**2
	#     w=w+d
	#     i_sort=i_sort_at(i_at)
	#     w_sort(i_sort)=w_sort(i_sort)+d
	#     n_at_sort(i_sort)=n_at_sort(i_sort)+1
    # 	    enddo
    # 	w=w/dfloat(n_at)
	# ww=ww+w
	#     do i_sort=1,n_sort
	#     if(n_at_sort(i_sort).gt.0) then
    # 	    w_sort(i_sort)=w_sort(i_sort)/dfloat(n_at_sort(i_sort))
	#     ww_sort(i_sort)=ww_sort(i_sort)+w_sort(i_sort)
	#     endif
	#     enddo
	# enddo
	# ww=ww/dfloat(n_write-i_time_diff)
    # 	dr_at(i_time_diff)=ww
	#     do i_sort=1,n_sort
	#     if(n_at_sort(i_sort).gt.0) then
	#     dr_at_sort(i_time_diff,i_sort)=ww_sort(i_sort)/dfloat(n_write-i_time_diff)
	#     endif
	#     enddo
	# enddo