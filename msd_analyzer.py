import os
import numpy as np
from md_format_converter_mi import *
from simple_msd_functions import *

print('Hello!')

# reading the structure in rv_at fromat from /dataw directory in order to get info about masses
path2rvat = os.getcwd()+'/dataw/'
rvats = os.listdir(path2rvat)
rvat = rvats[0]
st_parent = read_rv_at(path2rvat+rvat)


# reading the dumps of structures in r_at format (only coordinates of atoms)
path2dumps = os.getcwd()+'/r_at/'
dumps = os.listdir(path2dumps)
structures = []
for dmp in dumps:
    st = read_r_at(path2dumps+dmp)
    st.i_mass_at = st_parent.i_mass_at  # asigning atomic masses from parent rv_at structure
    st.i_at = st_parent.i_at            # asigning atomic indexes from parent rv_at structure
    st.i_type_at = st_parent.i_type_at  # asigning atomic types from parent rv_at structure
    st.n_type_at = st_parent.n_type_at  # asigning atomic types from parent rv_at structure
    st.type_at = st_parent.type_at  # asigning atomic types from parent rv_at structure
    st_numpy = convert_structure_to_numpy(st) # converting attributes of structure class into the numpy array
    structures.append(st_numpy)

# getting the list of structures with unwrapped positions
unwrapped_structures = get_unwrapped_structures(structures) # unwrapping the atomic coordinates in structures

# getting the list of structures with unwrapped positions and corrected center of mass 
cm_corrected_unwrapped_structures = []
for uwst in unwrapped_structures:
    cm_corrected_st = get_cm_corrected_structure(uwst)  # correcting the position of center of mass in structure 
    cm_corrected_unwrapped_structures.append(cm_corrected_st)

# getting the list of structures with corrected center of mass and unwrapped positions
# cm_corrected_unwrapped_structures = get_unwrapped_structures(cm_corrected_structures) # unwrapping the atomic coordinates in structures

dt = 75
msd_cm_corrected = calc_non_averaged_msd(cm_corrected_unwrapped_structures, dt)
msd_non_corrected = calc_non_averaged_msd(unwrapped_structures, dt)

print(msd_cm_corrected)
print(msd_non_corrected)