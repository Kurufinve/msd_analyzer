import numpy as np
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

    print(f"cm_x = {cm_x}")
    print(f"cm_y = {cm_y}")
    print(f"cm_z = {cm_z}")

    st_numpy_new = copy.deepcopy(st_numpy)
    st_numpy_new.r_at[:,0] = st_numpy.r_at[:,0] - cm_x
    st_numpy_new.r_at[:,1] = st_numpy.r_at[:,1] - cm_y
    st_numpy_new.r_at[:,2] = st_numpy.r_at[:,2] - cm_z

    return st_numpy_new
