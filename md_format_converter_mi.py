#!/usr/bin/env python

""" 
This module contains classes and functions for reading and writing files
in different formats for MD simulations including so-named "_at_" format 
for intrinsic MD codes from BSU laboratory of theoretical investigations 
and computer modeling.

The reading is made to the structure class
The writing is made from the structure class 

The module can be imported or executed as main program

"""

class structure():
    def __init__(self):
        self.n_at = None # number of atoms (integer)
        self.n_mark_at = 10 # number of atom markers (only for intrinsic codes)
        self.mark_at = [] # n_mark_at array of atom markers (only for intrinsic codes)
        self.i_at = [] # n_at array of atomic indexes (floats)
        self.r_at = [] # 3*n_at array of atomic coordinates (floats)
        self.f_at = [] # 3*n_at array of atomic forces (floats)
        self.v_at = [] # 3*n_at array of atomic velocities (floats)
        self.sizex = None # size of simulation cell in x direction (float)
        self.sizey = None # size of simulation cell in y direction (float)
        self.sizez = None # size of simulation cell in z direction (float)
        self.a_lattice3 = [1.0,1.0,1.0] # lattice constants (usually do not matter at all as sizex, sizey, sizez is sufficient)
        self.n_type_at = None # number of atomic types (integer)
        self.i_type_at = [] # n_at array of atomic types (integers)
        self.i_mass_at = [] # n_at array of atomic masses (floats)
        self.type_at = [] # n_type_at array of different atom types (integers)
        self.mass_at = [] # n_type_at array of masses of different atom types (floats)

    def __str__(self):
        nl = '\n'
        return  f"Structure with {self.n_at} atoms and {self.n_type_at} atom types\n" + \
                f"Sizes: {self.sizex} {self.sizey} {self.sizez}\n" + \
                f"Coordinates of atoms in angstroms:\n" + \
                f"{nl.join([str(r[0])+' '+str(r[1])+' '+str(r[2]) for r in self.r_at[0:5]])}\n" + \
                f" ... \n" + \
                f"{nl.join([str(r[0])+' '+str(r[1])+' '+str(r[2]) for r in self.r_at[-5:]])}\n"

    # def mass_at

def read_lmp_data(filename=None):

    st = structure()
    print('read_lmp_data filename:',filename)
    f = open(filename) 
    lmp_data = f.readlines()
    f.close()

    mass_reading_mode = False
    atom_reading_mode = False
    for line in lmp_data[1:]:
        if len(line.split()) == 0:
            continue # skip line if it is empty
        elif 'atoms' in line:
            # write the content of the line with the number of atoms
            st.n_at = int(line.split()[0]) 
        elif 'atom types' in line:
            # write the content of the line with the number of atom types
            st.n_type_at = int(line.split()[0]) 
        elif 'xlo' and 'xhi' in line:
            # write the content of the line with xlo xhi limits
            st.sizex = float(line.split()[1]) - float(line.split()[0]) # xhi - xlo
        elif 'ylo' and 'yhi' in line:
            # write the content of the line with ylo yhi limits
            st.sizey = float(line.split()[1]) - float(line.split()[0]) # yhi - ylo
        elif 'zlo' and 'zhi' in line:
            # write the content of the line with zlo zhi limits
            st.sizez = float(line.split()[1]) - float(line.split()[0]) # zhi - zlo            
        elif 'Masses' in line:
            mass_reading_mode = True
            atom_reading_mode = False
            continue
        elif 'Atoms' in line:
            mass_reading_mode = False
            atom_reading_mode = True    
            continue    

        if mass_reading_mode:
            st.type_at.append(int(line.split()[0]))
            st.mass_at.append(float(line.split()[1]))
            # st.type_mass_at   

        if atom_reading_mode:
            i_at = int(line.split()[0]) # read atomic number
            t_at = int(line.split()[1]) # read atomic type
            x_at = float(line.split()[2]) # read atomic x coordinate
            y_at = float(line.split()[3]) # read atomic y coordinate 
            z_at = float(line.split()[4]) # read atomic z coordinate 

            st.i_at.append(i_at)
            st.i_type_at.append(t_at)
            st.r_at.append([x_at,y_at,z_at])
            st.f_at.append([0.0,0.0,0.0]) 
            st.v_at.append([0.0,0.0,0.0]) 
            st.i_mass_at.append(st.mass_at[st.type_at.index(t_at)])

    return st


def read_rv_at(filename=None):

    st = structure()
    print('read_rv_at filename:',filename)
    # f = open(filename) 
    # rv_at = f.readlines()
    # f.close()

    f = open(filename)
    line = f.readline().split()
    st.n_at = int(line[0])
    st.n_mark_at = int(line[1])
    line = f.readline().split()
    size = [float(i) for i in line]
    st.sizex = size[0]
    st.sizey = size[1]
    st.sizez = size[2]
    line = f.readline().split()
    st.a_lattice3 = [float(i) for i in line]
    fline = f.read().split()
    f.close()

    iw=0
    # r_at=[]
    for i in range(st.n_at):
        x=float(fline[iw]); iw+=1
        y=float(fline[iw]); iw+=1
        z=float(fline[iw]); iw+=1
        st.r_at.append([x,y,z])

    # v_at=[]
    for i in range(st.n_at):
        vx=float(fline[iw]); iw+=1
        vy=float(fline[iw]); iw+=1
        vz=float(fline[iw]); iw+=1
        st.v_at.append([vx,vy,vz])

    # mass_at=[]
    for i in range(st.n_at):
        mass=float(fline[iw]); iw+=1
        st.i_mass_at.append(mass)

    # i_sort_at=[]
    for i in range(st.n_at):
        i_type=int(fline[iw]); iw+=1
        st.i_type_at.append(i_type)

    # num_at_r=[]
    for i in range(st.n_at):
        num_at=int(fline[iw]); iw+=1
        st.i_at.append(num_at)

    st.type_at = list(set(st.i_type_at))
    st.n_type_at = len(st.type_at)
    st.mass_at = list(set(st.i_mass_at))

    return st        


def read_r_at(filename=None):

    """ Function for reading files in simple r_at format for MSD calculations"""
    st = structure()
    print('read_r_at filename:',filename)

    f = open(filename)
    line = f.readline().split()
    st.n_at = int(line[0])
    st.n_mark_at = 0
    line = f.readline().split()
    size = [float(i) for i in line]
    st.sizex = size[0]
    st.sizey = size[1]
    st.sizez = size[2]
    st.a_lattice3 = [st.sizex, st.sizey, st.sizez]
    fline = f.read().split()
    f.close()

    iw=0
    # r_at=[]
    for i in range(st.n_at):
        x=float(fline[iw]); iw+=1
        y=float(fline[iw]); iw+=1
        z=float(fline[iw]); iw+=1
        st.r_at.append([x,y,z])

    # v_at=[]
    for i in range(st.n_at):
        vx=0.0
        vy=0.0
        vz=0.0
        st.v_at.append([vx,vy,vz])

    # mass_at=[]
    for i in range(st.n_at):
        mass=0.0
        st.i_mass_at.append(mass)

    # i_sort_at=[]
    for i in range(st.n_at):
        i_type=0.0
        st.i_type_at.append(i_type)

    # num_at_r=[]
    for i in range(st.n_at):
        num_at=i
        st.i_at.append(num_at)

    st.type_at = list(set(st.i_type_at))
    st.n_type_at = len(st.type_at)
    st.mass_at = list(set(st.i_mass_at))

    return st    


def write_rf_at(st,filename=None):
    """
    This function takes Structure() object as input.
    Forces will be zeros as they are unknown!

    """

    if filename == None:
        filename = 'at_rf'

    with open(filename,'w') as f:
        f.write('{0:5} 10\n'.format(st.n_at))
        f.write('{0:9.6f} {1:9.6f} {2:9.6f}  \n'.format(st.sizex, st.sizey,st.sizez))
        f.write('{0:9.6f} {1:9.6f} {2:9.6f}  \n'.format(1, 1, 1))

        for r_at_i in st.r_at:
            f.write('{0:9.6f}  {1:9.6f}  {2:9.6f}  '.format(r_at_i[0],r_at_i[1],r_at_i[2]))
        f.write('\n')

        for f_at_i in st.f_at:
            f.write('{0:9.6f}  {1:9.6f}  {2:9.6f}  '.format(f_at_i[0],f_at_i[1],f_at_i[2]))
        f.write('\n')      

        for i in st.i_mass_at:
            f.write(str(i)+' ')

        f.write('\n')
        for i in st.i_type_at:
            f.write(str(i)+' ')
        f.write('\n')   

        for i in range(st.n_at):
            f.write(str(i+1)+' ')
        f.write('\n')   

        for i in range(st.n_at):
            f.write('T ')
        f.write('\n')      

        for i in range(10):
            for j in range(st.n_at):
                f.write('F ')
            f.write('\n')     

        f.close()
    return 

def write_rv_at(st,filename=None):
    """
    This function takes Structure() object as input.
    Forces will be zeros as they are unknown!

    """

    if filename == None:
        filename = 'at_rv'

    with open(filename,'w') as f:
        f.write('{0:5} 10\n'.format(st.n_at))
        f.write('{0:9.6f} {1:9.6f} {2:9.6f}  \n'.format(st.sizex, st.sizey,st.sizez))
        f.write('{0:9.6f} {1:9.6f} {2:9.6f}  \n'.format(1, 1, 1))

        for r_at_i in st.r_at:
            f.write('{0:9.6f}  {1:9.6f}  {2:9.6f}  '.format(r_at_i[0],r_at_i[1],r_at_i[2]))
        f.write('\n')

        for v_at_i in st.v_at:
            f.write('{0:9.6f}  {1:9.6f}  {2:9.6f}  '.format(v_at_i[0],v_at_i[1],v_at_i[2]))
        f.write('\n')      

        for i in st.i_mass_at:
            f.write(str(i)+' ')

        f.write('\n')
        for i in st.i_type_at:
            f.write(str(i)+' ')
        f.write('\n')   

        for i in range(st.n_at):
            f.write(str(i+1)+' ')
        f.write('\n')   

        for i in range(st.n_at):
            f.write('T ')
        f.write('\n')      

        for i in range(10):
            for j in range(st.n_at):
                f.write('F ')
            f.write('\n')     

        f.close()


    return 

def write_lmp_data(st,filename=None):

    # Checking if n_type_at == type_at[-1]
    if st.n_type_at < st.type_at[-1]:
        type_at_old = st.type_at
        type_at_new = [i+1 for i in range(st.n_type_at)] 
        delta_type_at = [i-j for i,j in zip(type_at_old,type_at_new)]
        st.type_at = type_at_new
        i_type_at_old = st.i_type_at
        i_type_at_new = []
        # i_type_at_new = [i for i in st.i_type_at]
        for i_type in i_type_at_old:
            for type_at, delta in zip(type_at_old,delta_type_at):
                if i_type == type_at: 
                    i_type_at_new.append(i_type-delta)
        st.i_type_at = i_type_at_new


    # Cheking if cell is centered around 0:
    x_centered = False
    y_centered = False
    z_centered = False
    for r in st.r_at:
        x = r[0]; y = r[1]; z = r[2]
        if x < -0.5:
            x_centered = True
        if y < -0.5:
            y_centered = True
        if z < -0.5:
            z_centered = True


    with open(filename,'w') as f:  
        f.write('Start file for LAMMPS generated by md_format_converter_mi.py\n')
        f.write('{} atoms\n'.format(st.n_at))
        f.write('{} atom types\n'.format(st.n_type_at))
        f.write('\n')
        # writing sizes of cell
        if x_centered:
            f.write('{0:9.6f} {1:9.6f} xlo xhi\n'.format(-st.sizex/2,st.sizex/2))
        elif not x_centered:
            f.write('0.0 {0:9.6f} xlo xhi\n'.format(st.sizex))

        if y_centered:
            f.write('{0:9.6f} {1:9.6f} ylo yhi\n'.format(-st.sizey/2,st.sizey/2))
        elif not y_centered:
            f.write('0.0 {0:9.6f} ylo yhi\n'.format(st.sizey))

        if z_centered:
            f.write('{0:9.6f} {1:9.6f} zlo zhi\n'.format(-st.sizez/2,st.sizez/2))
        elif not z_centered:
            f.write('0.0 {0:9.6f} zlo zhi\n'.format(st.sizez))

        f.write('\n')
        f.write('Masses \n')
        f.write('\n')
        # for type_at, mass_at in zip(range(1,len(st.type_at)+1),st.mass_at):
        for type_at, mass_at in zip(st.type_at,st.mass_at):
            f.write('{0} {1:9.6f}\n'.format(type_at,mass_at))
        f.write('\n')
        f.write('Atoms \n')
        f.write('\n')
        for i_num, i_type, i_r_at in zip(st.i_at,st.i_type_at,st.r_at):
            f.write('{0} {1} {2:9.6f} {3:9.6f} {4:9.6f}\n'.format(i_num,i_type,i_r_at[0],i_r_at[1],i_r_at[2]))
        f.write('\n')
        f.close()

    return

def print_instructions():
    print('This program converts files in different formats for MD simulation codes\
 including the so-named "_at_" format for intrinsic MD codes\
 from BSU laboratory of theoretical investigations\
 and computer modeling authored by A.G. Lipnitskii.')
    print('\n')
    print('Usage: md_format_converter_mi.py mode infile outfile')
    print('\n')
    print('modes: ')
    print(' --lmpdata2rvat or -ldarv - convert LAMMPS data file to rv_at file (r means coordinates, v means velocities)')
    print(' --lmpdata2rfat or -ldarf - convert LAMMPS data file to rf_at file (r means coordinates, f means forces)')
    print(' --rvat2lmpdata or -arvld - convert rv_at file to LAMMPS data file')
    print(' --rfat2lmpdata or -arfld - convert rf_at file to LAMMPS data file')

if __name__=='__main__':
    import sys
    try:
        mode    = sys.argv[1]
        infile  = sys.argv[2]
        outfile = sys.argv[3]
    except:
        print_instructions()
        print('\n')
        print('Error: Not all arguments are pointed out!')
        print('\n')
        sys.exit(1)

    print('mode:', mode)
    print('infile:', infile)
    print('outfile:', outfile)

    if mode == '--lmpdata2rvat' or mode == '-ldarv':
        print('Converting from LAMMPS data file to rv_at file...')
        st = read_lmp_data(infile)
        write_rv_at(st,outfile)

    elif mode == '--lmpdata2rfat' or mode == '-ldarf':
        print('Converting from LAMMPS data file to rf_at file...')
        st = read_lmp_data(infile)
        write_rf_at(st,outfile)

    elif mode == '--rvat2lmpdata' or mode == '-arvld':
        print('Converting from rv_at file to LAMMPS data file...')
        st = read_rv_at(infile)
        write_lmp_data(st,outfile)        

    elif mode == '--rfat2lmpdata' or mode == '-arfld':
        print('Converting from rf_at file to LAMMPS data file...')
        st = read_rv_at(infile)
        write_lmp_data(st,outfile) 

    else:
        print('Error in arguments!')
        print_instructions()
        sys.exit(1)