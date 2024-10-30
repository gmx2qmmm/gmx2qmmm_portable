from collections.abc import Iterable

import numpy as np


def _flatten(x):

    '''
    ------------------------------ \\
    EFFECT: \\
    --------------- \\
    XX \\
    ------------------------------ \\
    INPUT: \\
    --------------- \\
    NONE \\
    ------------------------------ \\
    RETURN: \\
    --------------- \\
    NONE \\
    ------------------------------ \\
    '''

    """Replace deprecated ``compiler.ast.flatten``"""
    for e in x:
        if not isinstance(e, Iterable) or isinstance(e, str):
            yield e
        else:
            yield from _flatten(e)

def create_dict(label, info):

    '''
    ------------------------------ \\
    EFFECT: \\
    --------------- \\
    XX \\
    ------------------------------ \\
    INPUT: \\
    --------------- \\
    NONE \\
    ------------------------------ \\
    RETURN: \\
    --------------- \\
    NONE \\
    ------------------------------ \\
    '''

    out_dict = {}
    for i in range(len(label)):
        out_dict[label[i]] = info[i]
    return out_dict

def stepper(filename, step):

    '''
    ------------------------------ \\
    EFFECT: \\
    --------------- \\
    XX \\
    ------------------------------ \\
    INPUT: \\
    --------------- \\
    NONE \\
    ------------------------------ \\
    RETURN: \\
    --------------- \\
    NONE \\
    ------------------------------ \\
    '''

    """Move to more appropriate module"""
    if step == 0:
        return filename
    elif filename[-4:] == ".g96":
        return str(filename[:-4] + "." + str(step) + ".g96")
    elif filename[-4:] == ".gro":
        return str(filename[:-4] + "." + str(step) + ".gro")
    elif filename[-13:] == ".pointcharges":
        return str(filename[:-13] + "." + str(step) + ".pointcharges")
    else:
        if len(filename) > 7:
            if filename[-7:] == ".fort.7":
                new_filename = str(filename[:-7] + "." + str(step) + ".fort.7")
                return new_filename
        for i in range(0, len(filename)):
            buffer = 0
            if filename[i] == ".":
                buffer = i
        new_filename = str(filename[:buffer] + "." + str(step) + filename[buffer:])
        return new_filename


def make_xyzq(geo, chargevec):

    '''
    ------------------------------ \\
    EFFECT: \\
    --------------- \\
    XX \\
    ------------------------------ \\
    INPUT: \\
    --------------- \\
    NONE \\
    ------------------------------ \\
    RETURN: \\
    --------------- \\
    NONE \\
    ------------------------------ \\
    '''

    xyzq=np.zeros((len(chargevec),4))
    try:
        xyzq[:,0]=geo[0::3]
        xyzq[:,1]=geo[1::3]
        xyzq[:,2]=geo[2::3]
        xyzq[:,3]=chargevec
    except:
        print("Error: Can't make XYZ-charge matrix. Maybe the number of atoms is different in the structure and the topology file ?!")
    return xyzq

def get_coordinates_linkatoms_angstrom(array_xyzq, list_atoms_qm, list_atoms_m1, list_connectivity_topology, prev_scalefacs) -> list:
    '''
    ------------------------------ \\
    EFFECT: \\
    ---------------
    Calculate The Coordinates Of Linkatoms \\
    ------------------------------ \\
    INPUT: \\
    --------------- \\
    # XX AJ what is prev_scalefacs? \\
    XX AJ add input variables when we decided where to keep this function
    ------------------------------ \\
    RETURN: \\
    --------------- \\
    list_coordinates_linkatoms: list -> List Of Coordinates Of Linkatoms And Scaling Factor \\
    ------------------------------ \\
    '''

    list_coordinates_linkatoms = []
    list_atoms_linkpairs = []

    #   Read All QM-M1 Pairs
    for index_m1 in list_atoms_m1:
        list_linkpair = []
        for i in range(0, len(list_connectivity_topology)):
            if int(index_m1) in np.array(list_connectivity_topology[i]).astype(int):
                if int(index_m1) == int(list_connectivity_topology[i][0]):
                    for j in range(0, len(list_connectivity_topology[i])):
                        if int(list_connectivity_topology[i][j]) in np.array(list_atoms_qm).astype(int):
                            list_linkpair = [int(list_connectivity_topology[i][j]), int(index_m1)]
                            list_atoms_linkpairs.append(list_linkpair)
                            break
                    break
                else:
                    if int(list_connectivity_topology[i][0]) in np.array(list_atoms_qm).astype(int):
                        list_linkpair = [int(list_connectivity_topology[i][0]), int(index_m1)]
                        list_atoms_linkpairs.append(list_linkpair)
                        break   # XX AJ I don't understand why we break here. What if m1 also is in another connectivity list connected to other QM atoms?

    #   Calculate The Linkatom Coordinates
    count = 0
    for index_m1 in list_atoms_linkpairs:
        list_coordinates_linkatom = []
        array_coordinates_atom_qm = array_xyzq[index_m1[0]-1][:3]
        array_coordinates_atom_m1 = array_xyzq[index_m1[1]-1][:3]
        array_vector_m1_qm = array_coordinates_atom_m1 - array_coordinates_atom_qm

        #   Scaling The M1-QM Vector
        scalefac = 0.71290813568205  # ratio between B3LYP/6-31G* optimum of C-C in butane vs C-Link relaxed butane
        if len(prev_scalefacs) != 0:
            scalefac = prev_scalefacs[count][3]
        array_vector_m1_qm_scaled = array_vector_m1_qm * scalefac
        list_coordinates_linkatom = list(array_vector_m1_qm_scaled + array_coordinates_atom_qm)
        list_coordinates_linkatom.append(scalefac)
        list_coordinates_linkatoms.append(list_coordinates_linkatom)
        count += 1

    return list_coordinates_linkatoms

def filter_xyzq(xyzq, atom_list, coordinates=True, charges=True):

    '''
    ------------------------------ \\
    EFFECT: \\
    --------------- \\
    XX \\
    ------------------------------ \\
    INPUT: \\
    --------------- \\
    NONE \\
    ------------------------------ \\
    RETURN: \\
    --------------- \\
    NONE \\
    ------------------------------ \\
    '''

    # XX AJ I'm not happy with how this function works (combining lists and arrays) but that's because of different inputs
    python_atom_list = [index - 1 for index in atom_list]
    if isinstance(xyzq, np.ndarray):
        xyzq = xyzq.tolist()
    filtered_xyzq = [xyzq[index] for index in python_atom_list]
    begin = 0 if coordinates else 3
    end = 4 if charges else 3
    filtered_information = np.array(filtered_xyzq)[:, begin:end]
    return filtered_information

def normalized_vector(vec):

    '''
    ------------------------------ \\
    EFFECT: \\
    normalizing a vector \\
    --------------- \\
    INPUT: \\
    vec: np.array
    ------------------------------ \\
    RETURN: \\
    np.array -> normalized vector \\
    --------------- \\
    '''
    return np.array(vec) / np.linalg.norm(np.array(vec))

def mask_atoms(list_2d, list_atoms_2keep):

    '''
    ------------------------------ \\
    EFFECT: \\
    --------------- \\
    XX \\
    ------------------------------ \\
    INPUT: \\
    --------------- \\
    NONE \\
    ------------------------------ \\
    RETURN: \\
    --------------- \\
    NONE \\
    ------------------------------ \\
    '''

    list_atoms_indices = list(np.array(list_atoms_2keep) - 1)
    masked_list = np.zeros_like(list_2d)
    masked_list[list_atoms_indices] = list_2d[list_atoms_indices]

    return masked_list
