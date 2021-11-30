import os
import numpy as np
from md_format_converter_mi import *
from simple_msd_functions import *

# print('Hello!')

""" Start block with reading data from r_at files and writing it to csv and excel"""

# # reading the structure in rv_at fromat from /dataw directory in order to get info about masses
# path2rvat = os.getcwd()+'/dataw/'
# rvats = os.listdir(path2rvat)
# rvat = rvats[0]
# st_parent = read_rv_at(path2rvat+rvat)


# # reading the dumps of structures in r_at format (only coordinates of atoms)
# path2dumps = os.getcwd()+'/r_at/'
# dumps = os.listdir(path2dumps)
# structures = []
# for dmp in dumps:
#     st = read_r_at(path2dumps+dmp)
#     st.i_mass_at = st_parent.i_mass_at  # asigning atomic masses from parent rv_at structure
#     st.i_at = st_parent.i_at            # asigning atomic indexes from parent rv_at structure
#     st.i_type_at = st_parent.i_type_at  # asigning atomic types from parent rv_at structure
#     st.n_type_at = st_parent.n_type_at  # asigning atomic types from parent rv_at structure
#     st.type_at = st_parent.type_at  # asigning atomic types from parent rv_at structure
#     st_numpy = convert_structure_to_numpy(st) # converting attributes of structure class into the numpy array
#     structures.append(st_numpy)

# # getting the list of structures with unwrapped positions
# unwrapped_structures = get_unwrapped_structures(structures) # unwrapping the atomic coordinates in structures

# # getting the list of structures with unwrapped positions and corrected center of mass 
# cm_corrected_unwrapped_structures = []
# for uwst in unwrapped_structures:
#     cm_corrected_st = get_cm_corrected_structure(uwst)  # correcting the position of center of mass in structure 
#     cm_corrected_unwrapped_structures.append(cm_corrected_st)

# # getting the list of structures with corrected center of mass and unwrapped positions
# # cm_corrected_unwrapped_structures = get_unwrapped_structures(cm_corrected_structures) # unwrapping the atomic coordinates in structures

# # calculating the non-averaged and averaged msd 
# dt = 7.5E-11
# msd_cm_corrected = calc_non_averaged_msd(cm_corrected_unwrapped_structures, dt)
# msd_av_cm_corrected = calc_averaged_msd(cm_corrected_unwrapped_structures, dt)
# msd = calc_non_averaged_msd(unwrapped_structures, dt)
# msd_av = calc_averaged_msd(unwrapped_structures, dt)

# # writing the data into csv and excel files
# msd.to_csv('msd.csv')
# msd_av.to_csv('msd_averaged.csv')
# msd_cm_corrected.to_csv('msd_cm_corrected.csv')
# msd_av_cm_corrected.to_csv('msd_averaged_cm_corrected.csv')

# msd.to_excel('msd.xlsx')
# msd_av.to_excel('msd_averaged.xlsx')
# msd_cm_corrected.to_excel('msd_cm_corrected.xlsx')
# msd_av_cm_corrected.to_excel('msd_averaged_cm_corrected.xlsx')

""" End block with reading data from r_at files and writing it to csv and excel"""


# # reading the data from csv files
# # import pandas as pd
# # msd_cm_corrected = pd.read_csv('msd_cm_corrected.csv')
# # msd_av_cm_corrected = pd.read_csv('msd_averaged_cm_corrected.csv')
# # msd_cm_corrected.to_excel('msd_cm_corrected.csv')
# # msd_av_cm_corrected.to_excel('msd_averaged_cm_corrected.csv')


# reading the data from excel files
import pandas as pd
# msd_cm_corrected = pd.read_excel('msd_cm_corrected.xlsx')
msd_av_cm_corrected = pd.read_excel('msd_averaged_cm_corrected.xlsx')

fit_and_plot_msd_linear(msd_av_cm_corrected)

# print(msd_cm_corrected.describe())
# print(msd_av_cm_corrected.describe())
# print(msd_non_corrected)