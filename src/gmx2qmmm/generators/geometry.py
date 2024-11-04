#   // INITIAL DESCRIPTION //
"""Short Module Description; Reference To Readme"""

#   // MEATDATA //
__author__ = 'Alina Jansen'
__date__ = '2024-03-25'

#   // IMPORTS //

#   Imports Of Existing Libraries
import re
import numpy as np

#   Imports From Existing Libraries

#   Imports Of Custom Libraries

#   Imports From Custom Libraries
from gmx2qmmm.logging import Logger
from gmx2qmmm.generators._helper import filter_xyzq, _flatten


#   // TODOS & NOTES //
#   TODO:
#   - change 'make_xyzq' and 'make_xyzq_io' to be one function
#   NOTE: For now I'm just randomly keeping functions regarding geometry stuff here. Sorting into a class later (AJ)

def readg96(inp):
    '''
    ------------------------------
    EFFECT: \\
    ---------------
    reads coordinates of all atoms from a g96 file
    ------------------------------
    INPUT: \\
    ---------------
    inp: string, name of g96 file
    ------------------------------
    RETURN: \\
    ---------------
    coords: list of coordinates
    ------------------------------
    '''
    coords = []
    with open(inp) as ifile:
        count = 0
        for line in ifile:
            count += 1
            if count == 4:
                break
        count = 1
        for line in ifile:
            match = re.search(
                r"^(.{5})\s(.{5})\s(.{5})\s(.{6})\s*([-]*\d+\.*\d*)\s*([-]*\d+\.*\d*)\s*([-]*\d+\.*\d*)",
                line,
                flags=re.MULTILINE,
            )
            if match:
                coords.append(float(match.group(5)) * 10.0)
                coords.append(float(match.group(6)) * 10.0)
                coords.append(float(match.group(7)) * 10.0)
            else:
                break
    return coords

def readgeo(inp):
    coords = []
    n_a = 0
    with open(inp) as ifile:
        for line in ifile:
            break
        for line in ifile:
            match = re.search(r"^\s*(\d+)", line, flags=re.MULTILINE)
            if match:
                n_a = int(match.group(1))
                break
            else:
                print(".gro is corrupt (no number of atoms found, second line). Exiting.")
                exit(1)
        count = 1
        for line in ifile:
            match = re.search(
                r"^(.{5})(.{5})(.{5})(.{5})\s*([-]*\d+\.*\d*)\s*([-]*\d+\.*\d*)\s*([-]*\d+\.*\d*)",
                line,
                flags=re.MULTILINE,
            )
            if match:
                coords.append(float(match.group(5)) * 10.0)
                coords.append(float(match.group(6)) * 10.0)
                coords.append(float(match.group(7)) * 10.0)
            else:
                print(".gro is corrupt. Exiting.")
                print("Last line:")
                print(line)
                exit(1)
            count += 1
            if count > n_a:
                break
    return coords

def make_xyzq(geo, chargevec):
    xyzq=np.zeros((len(chargevec),4))
    try:
        xyzq[:,0]=geo[0::3]
        xyzq[:,1]=geo[1::3]
        xyzq[:,2]=geo[2::3]
        xyzq[:,3]=chargevec
    except:
        print("Error: Can't make XYZ-charge matrix. Maybe the number of atoms is different in the structure and the topology file ?!")
    return xyzq

def make_xyzq_io(geo, chargevec, outerlist):
    xyzq=np.zeros((len(chargevec),4))
    try:
        xyzq[:,0]=geo[0::3]
        xyzq[:,1]=geo[1::3]
        xyzq[:,2]=geo[2::3]
        outer=np.ones(len(chargevec))
        outer[np.array(outerlist)-1]=0          #the charge of the atoms of the outerlist is set to 0
        xyzq[:,3]=np.array(chargevec)*outer
    except:
        print("Error: Can't make XYZ-charge matrix. Maybe the number of atoms is different in the structure and the topology file ?!")
    return xyzq

def propagate_dispvec(propagator, xyzq, all_forces, float_force_max, stepsize, curr_step, scan_flag=False, scan_atoms=[]):
    # XX AJ xyzq is only used in scan and BFGS, I will only do that later
    dispvec = []
    total_force=all_forces[-1]
    # last_forces=all_forces[-2]
    # if scan_flag :
    #     scan_type = len(scan_atoms)
    #     if scan_type == 2 :
    #         clean_force = scan_adj_force(clean_force, xyzq, scan_atoms, 'R')
    #         if (propagator == "BFGS") and (len(old_clean_force) > 1):
    #             old_clean_force = scan_adj_force(old_clean_force, xyzq, scan_atoms, 'R')

    #     elif scan_type == 3 :
    #         clean_force = scan_adj_force(clean_force, xyzq, scan_atoms, 'A')

    #     elif scan_type == 4 :
    #         clean_force = scan_adj_force(clean_force, xyzq, scan_atoms, 'B')

    max_abs_value_index  = np.argmax(np.abs((np.array(total_force)).flatten()))
    maxatom = max_abs_value_index // 3
    maxcoord = max_abs_value_index % 3

    # logger(
    #     logfile,
    #     str(
    #         "Maximum force is "
    #         + str(float(maxforce))
    #         + " a.u. at coord "
    #         + str(int(maxatom) + 1)
    #         + "/"
    #         + str(int(maxcoord) + 1)
    #         + ".\n"
    #     ),
    # )

    #   Always Use Steepest Descent For The First Propagation
    if propagator == "STEEP" or curr_step < 2:
        # logger(logfile, "Propagate with steepest descent...\n") or log propagate with steepest descent for first cycle if propagator not steep
        normalized_forces = np.array(total_force)*float(stepsize)/abs(float(float_force_max))
        dispvec = normalized_forces

    elif propagator == "CONJGRAD" :
        # logger(logfile, "Propagate with conjugate gradient...\n")
        # Fletcher-Reeves
        _flattened = list(_flatten(total_force))
        corr_fac = np.array(_flattened).dot(np.array(_flattened))
        _flattened = list(_flatten(last_forces))
        corr_fac /= np.array(_flattened).dot(np.array(_flattened))

        corr_length = np.array(total_force)
        if curr_step > 1:   #0 step for the inital SP
            corr_length = np.array(total_force) + corr_fac * corr_length
        displacement = corr_length*float(stepsize)
        dispvec = displacement.tolist()

        # logger(
        #     logfile,
        #     str(
        #         "Effective step at maximum force coord is "
        #         + str(float(dispvec[maxatom][maxcoord]))
        #         + " a.u.\n"
        #     ),
        # )

    elif propagator == "BFGS":
        pass    # XX AJ will do BFGS in the end
        # logger(logfile, "Propagate with BFGS method...\n")
        # coords = np.array(new_xyzq)[:, 0:3]
        # old_coords = np.array(xyzq)[:, 0:3]
        # # old_hessian = np.loadtxt("bfgs_hessian.txt")
        # old_hessian = sp.load_npz('bfgs_hessian.npz')
        # gradient = np.array(total_force)
        # old_gradient = np.array(last_forces)
        # hessian, hesseig, warning_flag = get_approx_hessian(
        #     coords, old_coords, gradient, old_gradient, old_hessian, logfile
        # )
        # # np.savetxt("bfgs_hessian.txt", hessian)
        # sp.save_npz('bfgs_hessian.npz', hessian)

        # coords = coords.reshape(
        #     3 * len(coords), 1
        # )  # reshape coords to use in dot products
        # gradient = gradient.reshape(
        #     3 * len(gradient), 1
        # )  # reshape grad to use in dot products

        # # direc = -sp.linalg.inv(hessian).dot(gradient) #XX CSC
        # direc = -sp.linalg.spsolve(hessian.tocsc(), sp.csc_matrix(gradient)) #XX CSC
        # direc = direc.reshape(int(len(coords) / 3), 3)
        # displacement = np.array(direc)*float(stepsize)
        # dispvec = displacement.tolist()
        dispvec *= 0.052917721
    return dispvec

# def add_displacement_vector():

def read_gmx_structure_header(file):
    with open(file, 'r') as structure_file:
        header = ''.join([next(structure_file) for _ in range(4)])

    return header

def read_gmx_structure_atoms(file):
    list_atom_information = []
    with open(file, 'r') as structure_file:
        for _ in range(4):
            next(structure_file)

        for line in structure_file:
                match = re.search(
                    r"^(.{5})\s(.{5})\s(.{5})\s(.{6})\s*([-]*\d+\.*\d*)\s*([-]*\d+\.*\d*)\s*([-]*\d+\.*\d*)",
                    line,
                    flags=re.MULTILINE,
                )
                if match:
                    list_atom_information.append(str(str(match.group(1))
                        + " "
                        + str(match.group(2))
                        + " "
                        + str(match.group(3))
                        + " "
                        + str(match.group(4))))

    return list_atom_information

def read_gmx_box_vectors(file):
    with open(file, 'r') as structure_file:
        for line in structure_file:
            match = re.search(r"^BOX\s*\n", line, flags=re.MULTILINE)
            if match:
                for line in structure_file:
                    match = re.search(
                        r"^\s*(\d*\.\d+)\s+(\d*\.\d+)\s+(\d*\.\d+)",
                        line,
                        flags=re.MULTILINE,
                    )
                    if not match:
                        # logger(
                        #     logfile,
                        #     "\n\nError: In "
                        #     + str(gro)
                        #     + " box vectors were expected but not found. Exiting. Line was:\n",
                        # )
                        # logger(logfile, line)
                        exit(1)
                    else:
                        list_vectors_box = [
                            float(match.group(1)),
                            float(match.group(2)),
                            float(match.group(3)),
                        ]
                        break
    return list_vectors_box

def write_g96(filename, header, atom_information, coordinates, box_vectors):
    with open(filename, 'w') as str_file:
        str_file.write(header)
        for int_atom_number, str_atom_information in enumerate(atom_information):
            str_file.write(
                str_atom_information
                + " {:>15.9f} {:>15.9f} {:>15.9f}\n".format(
                            coordinates[int_atom_number][0],
                            coordinates[int_atom_number][1],
                            coordinates[int_atom_number][2]
                )
            )
        str_file.write("END\nBOX\n")
        str_file.write(str(" {:>15.9f} {:>15.9f} {:>15.9f}\n".format(box_vectors[0], box_vectors[1], box_vectors[2])))
        str_file.write("END")
