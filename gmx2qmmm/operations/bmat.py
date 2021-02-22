import numpy as np
from scipy import linalg
from gmx2qmmm._helper import logger

#distance or angles
def length(coord_A, coord_B):
    import numpy as np
    return np.linalg.norm(coord_A - coord_B)

def angle(coords_A, coords_B, coords_C):
    vba = coord_A - coord_B
    vbc = coord_C - coord_B
    return np.arccos(vba.dot(vbc) / length(vba,vbc))

def dihedral(coord_A, coord_B, coord_C, coord_D):
    v0 = -1*(coord_B - coord_A)
    v1 = coord_C - coord_B
    v2 = coord_D - coord_C
    v0xv1 = np.cross(v0, v1)
    v1xv2 = np.cross(v2, v1)
    v0xv1_x_v1xv2 = np.cross(v0xv1, v1xv2)
    y = np.dot(v0xv1_x_v1xv2, v1) * (1.0 / np.linalg.norm(v1) )
    x = np.dot(v0xv1, v1xv2)
    return np.arctan2(y, x)

# unit vector
def e_vector(atom_1, atom_2):
    return (atom_2 - atom_1) / length(atom_1, atom_2)


def length_S_vector(two_atom_list):
    e12 = e_vector(two_atom_list[0], two_atom_list[1])
    s_vector = np.append([-e12], [e12])
    return s_vector

def angle_S_vector(three_atom_list):
    A, B, C = three_atom_list[:3]
    angle = angle(A, B, C)
    length_1_2 = length(A, B)
    length_2_3 = length(B, C)
    e_vectors = []
    for id_a, a in enumerate(three_atom_list):
        for id_b, b in enumerate(three_atom_list):
            if id_a == id_b:
                continue
            e_vec = e_vector(a, b)
            e_vectors.append(e_vec)
    s1 = ((np.cos(angle) * e_vectors[2]) - e_vectors[3]) / (length_1_2 * np.sin(angle))
    s3 = ((np.cos(angle) * e_vectors[3]) - e_vectors[2]) / (length_2_3 * np.sin(angle))
    s2 = -1 * (s1 + s3)
    s_vector = np.concatenate((s1, s2, s3), axis=None)
    return s_vector

def dihedral_S_vector(four_atom_list):
    A, B, C, D = four_atom_list[:4]
    angle123 = angle(A, B, C)
    angle234 = angle(B, C, D)   
    length_1_2 = length(A, B)
    length_2_3 = length(B, C)
    length_3_4 = length(C, D)
    e_vectors = []
    for id_a, a in enumerate(four_atom_list):
        for id_b, b in enumerate(four_atom_list):
            if id_a == id_b:
                continue
            e_vec = e_vector(a, b)
            e_vectors.append(e_vec)
    s1 = (-1 * (np.cross(e_vectors[0], e_vectors[4])) / (length_1_2*(np.sin(angle123)**2)))
    s2 = ((length_2_3 - length_1_2 * np.cos(angle123) / (length_2_3 * length_1_2 * np.sin(angle123))) * ((np.cross(e_vectors[0], e_vectors[4])) / np.sin(angle123)) + (np.cos(angle234) / (length_2_3 * np.sin(angle234))) * ((np.cross(e_vectors[11], e_vectors[7])) / np.sin(angle234)))
    s3 = ((length_2_3 - length_3_4 * np.cos(angle234) / (length_2_3 * length_3_4 * np.sin(angle234))) * ((np.cross(e_vectors[11], e_vectors[7])) / np.sin(angle234)) + (np.cos(angle123) / (length_2_3 * np.sin(angle123))) * ((np.cross(e_vectors[0], e_vectors[4])) / np.sin(angle123)))
    s4 = (-1 * (np.cross(e_vectors[11], e_vectors[7])) / (length_3_4*(np.sin(angle234)**2)))
    s_vector = np.concatenate((s1, s2, s3, s4), axis=None)
    return s_vector

def full_B_matrix(atom_coord_list, logfile):
    
    if len(atom_coord_list) < 2:
        logger(logfile, "No reflections possible for one atom. At least two are required.\n")
        return 0
    
    tot_s_vectors = 3 * len(atom_coord_list) - 6
    logger(logfile, "Number of atoms found: %d\n"%len(atom_coord_list))
    logger(logfile, "Minimum DOF for %s atoms are: %d"%(len(atom_coord_list),tot_s_vectors))
    num_s_vectors = range(1, tot_s_vectors + 1)
    if len(atom_coord_list) == 2:
        num_s_vectors = [1]
    s_vectors = []
    for id, svec in enumerate(num_s_vectors):
        if svec == 1:
            array = np.zeros((len(atom_coord_list) * 3))
            twoatomlist = atom_coord_list[:2]
            s_vec = length_S_vector(twoatomlist)
            array[0:3] = s_vec[0:3] 
            array[3:6] = s_vec[3:6]
            s_vectors.append(array)
            continue
        if svec == 2:
            array = np.zeros((len(atom_coord_list) * 3))
            threeatomlist = atom_coord_list[:3]
            s_vec = angle_S_vector(threeatomlist)
            array[0:3] = s_vec[0:3]
            array[3:6] = s_vec[3:6]
            array[6:9] = s_vec[6:]
            s_vectors.append(array)
            continue
        if svec == 3:
            array = np.zeros((len(atom_coord_list) * 3))
            twoatomlist = atom_coord_list[0:2]
            s_vec = length_S_vector(twoatomlist)
            array[3:6] = s_vec[0:3]
            array[6:9] = s_vec[3:]
            s_vectors.append(array)
            continue
        if (svec - 1) % 3 == 0:
            #angles
            array = np.zeros((len(atom_coord_list) * 3))
            atom_index = int(((svec + 2) / 3) + 1 )
            initial_index = int(svec + 8)
            threeatomlist = atom_coord_list[atom_index - 2: atom_index + 1]
            s_vec = angle_S_vector(threeatomlist)
            if svec - num_s_vectors[-3] == 0:
                array[initial_index - 9:] = s_vec[:]
            else:
                array[initial_index - 9: initial_index] = s_vec[:]
            s_vectors.append(array)
        if (svec - 2) % 3 == 0:
            #bonds ERROR
            array = np.zeros((len(atom_coord_list) * 3))
            atom_index = int(((svec + 1) / 3) + 2)
            initial_index = int(svec + 7)
            twoatomlist = atom_coord_list[atom_index - 2: atom_index]
            s_vec = length_S_vector(twoatomlist)
            if svec - num_s_vectors[-2] == 0:
                array[initial_index - 6:] = s_vec[:]
            else:
                array[initial_index - 6: initial_index] = s_vec[:]
            s_vectors.append(array)
        if svec % 3 == 0:
            #dihedrals
            array = np.zeros((len(atom_coord_list) * 3))
            atom_index = int((svec / 3) + 1)
            initial_index = int(svec + 6)
            fouratomlist = atom_coord_list[atom_index - 3: atom_index + 1]
            s_vec = dihedral_S_vector(fouratomlist)
            if svec - num_s_vectors[-1] == 0:
                array[initial_index - 12:] = s_vec[:]
            else:
                array[initial_index - 12: initial_index] = s_vec[:]
            s_vectors.append(array)
    bmat = np.stack(s_vectors[:])
    return bmat