# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

# this will calculate the expansion of a PLANAR pi-system
# requires input ground state structure, excited state structure (g09 opt log files, same order of atoms), will expect pi atoms to be the starting atoms in the list.
# also needs size of pi system (number of atoms)
# produces an output file
# needs input 1, input 2, size (int) and output

__author__ = "jangoetze"
__date__ = "$06-Feb-2018 12:45:17$"

import math
import sys

import numpy as np

from gmx2qmmm.operations import xyz_zmat_g16RevA_02 as zmatlib
from gmx2qmmm.operations import geo_xyz_g09RevA_02_log as getxyz
from gmx2qmmm.operations import geo_from_xyz as extrxyz


def find_neighbor(coords, curr_atom, found_list, basedir):
    distmat = zmatlib.make_distmat(coords)
    mindist = 100000.0
    curr_min = -1
    for i in range(0, len(coords) / 3):
        if i == curr_atom or (i + 1 in found_list):
            continue
        else:
            if distmat[curr_atom][i] < mindist:
                mindist = distmat[curr_atom][i]
                curr_min = i
    return curr_min + 1


def create_pi_conn(coords, size):
    conn = []
    for i in range(0, int(size)):
        conn_line = []
        found_neighbors = []
        for j in range(0, 3):  # pi atoms have three neighbours
            new_neighbor = find_neighbor(coords, i, found_neighbors)
            found_neighbors.append(new_neighbor)
        if len(found_neighbors) < 3:
            print("Something went wrong with the neighbor search. Exiting.")
            exit(1)
        for element in found_neighbors:
            if element > size:
                continue
            else:
                conn_line.append(element)
        conn.append(conn_line)
    return conn


def translate(c, a):  # unclustered coords
    orig = [c[a * 3 + 0], c[a * 3 + 1], c[a * 3 + 2]]
    new_c = []
    for i in range(0, len(c) / 3):
        for j in range(0, 3):
            coord = float(c[i * 3 + j]) - float(orig[j])
            new_c.append(coord)
    return new_c


def find_conn_element(a, b, conns, size):
    ele = -1
    for i in range(0, len(conns[a])):
        if int(conns[a][i]) == b + 1 or int(conns[a][i]) > int(size):
            continue
        else:
            ele = conns[a][i] - 1
            break
    if ele == -1:
        for i in range(0, len(conns[b])):
            if int(conns[b][i]) == a + 1 or int(conns[b][i]) > int(size):
                continue
            else:
                ele = conns[b][i] - 1
                break
    if ele == -1:
        print("No third atom found for rotation! Exiting.")
        exit(0)
    return ele


def uvec(vec):
    v = np.array(vec) / np.linalg.norm(np.array(vec))
    return v


def rotate3dZ180(c):
    new_coords = []
    rotmat = [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]
    for i in range(0, len(c) / 3):
        curr = []
        for j in range(0, 3):
            curr.append(c[i * 3 + j])
        new_curr = np.dot(rotmat, curr)
        new_coords.extend(new_curr)
    return new_coords


def rotate3d(c, rotvec, angle):
    new_coords = []
    rotmat = [
        [
            np.cos(angle) + rotvec[0] * rotvec[0] * (1.0 - np.cos(angle)),
            rotvec[0] * rotvec[1] * (1.0 - np.cos(angle)) - rotvec[2] * np.sin(angle),
            rotvec[0] * rotvec[2] * (1.0 - np.cos(angle)) + rotvec[1] * np.sin(angle),
        ],
        [
            rotvec[0] * rotvec[1] * (1.0 - np.cos(angle)) + rotvec[2] * np.sin(angle),
            np.cos(angle) + rotvec[1] * rotvec[1] * (1.0 - np.cos(angle)),
            rotvec[1] * rotvec[2] * (1.0 - np.cos(angle)) - rotvec[0] * np.sin(angle),
        ],
        [
            rotvec[0] * rotvec[2] * (1.0 - np.cos(angle)) - rotvec[1] * np.sin(angle),
            rotvec[1] * rotvec[2] * (1.0 - np.cos(angle)) + rotvec[0] * np.sin(angle),
            np.cos(angle) + rotvec[2] * rotvec[2] * (1.0 - np.cos(angle)),
        ],
    ]
    for i in range(0, len(c) / 3):
        curr = []
        for j in range(0, 3):
            curr.append(float(c[i * 3 + j]))
        if np.linalg.norm(curr) > 0.0000001:
            new_curr = np.dot(np.array(rotmat), np.array(uvec(curr))) * np.linalg.norm(curr)
            new_coords.extend(new_curr)
        else:
            new_coords.extend(curr)
    return new_coords


def angle3d(v1, v2):

    if np.linalg.norm(v1) == 0.0 or np.linalg.norm(v2) == 0.0:
        print("Requested angle between vectors of length 0, which is impossible. Exiting.")
        exit(1)
    angle = np.arcos(np.dot(uvec(v1), uvec(v2)))
    return angle


def rotate_planar(coords, a, b, c):
    rot_coords = []
    a_c = [coords[a * 3 + 0], coords[a * 3 + 1], coords[a * 3 + 2]]
    b_c = [coords[b * 3 + 0], coords[b * 3 + 1], coords[b * 3 + 2]]
    c_c = [coords[c * 3 + 0], coords[c * 3 + 1], coords[c * 3 + 2]]
    ab = np.array(b_c) - np.array(a_c)
    if angle3d(uvec(ab), [1.0, 0.0, 0.0]) > (math.pi - 0.01):
        curr_coords = rotate3dZ180(coords)
        a_c = [curr_coords[a * 3 + 0], curr_coords[a * 3 + 1], curr_coords[a * 3 + 2]]
        b_c = [curr_coords[b * 3 + 0], curr_coords[b * 3 + 1], curr_coords[b * 3 + 2]]
        c_c = [curr_coords[c * 3 + 0], curr_coords[c * 3 + 1], curr_coords[c * 3 + 2]]
        ab = np.array(b_c) - np.array(a_c)
        coords = curr_coords
    count = 0
    while angle3d(uvec(ab), [1.0, 0.0, 0.0]) > 0.000001:
        rotvec = uvec(np.cross(uvec(ab), [1.0, 0.0, 0.0]))
        # print angle3d(uvec(ab),[1.,0.0,0.0])
        curr_coords = rotate3d(coords, rotvec, angle3d(uvec(ab), [1.0, 0.0, 0.0]))
        a_c = [curr_coords[a * 3 + 0], curr_coords[a * 3 + 1], curr_coords[a * 3 + 2]]
        b_c = [curr_coords[b * 3 + 0], curr_coords[b * 3 + 1], curr_coords[b * 3 + 2]]
        c_c = [curr_coords[c * 3 + 0], curr_coords[c * 3 + 1], curr_coords[c * 3 + 2]]
        ab = np.array(b_c) - np.array(a_c)
        coords = curr_coords
        count += 1
        if count > 5:
            break
    bc = np.array(c_c) - np.array(b_c)
    if angle3d(uvec(ab), uvec(bc)) < 0.000001:
        print("Collinear pi system detected. This is not a properly defined system (linear bonds?). Exiting.")
        exit(0)
    torsion = math.atan2(
        np.dot(
            np.cross(np.cross([0.0, -1.0, 0.0], uvec(ab)), np.cross(uvec(ab), uvec(bc))),
            uvec(ab),
        ),
        np.dot(np.cross([0.0, -1.0, 0.0], uvec(ab)), np.cross(uvec(ab), uvec(bc))),
    )
    rot_coords = rotate3d(coords, uvec(ab), -1.0 * torsion)
    return rot_coords


def get_expansion(c1, c2, a, b, conns, size):
    ct1 = translate(c1, a)
    ct2 = translate(c2, a)
    c = find_conn_element(a, b, conns, size)
    cr1 = rotate_planar(ct1, a, b, c)
    cr2 = rotate_planar(ct2, a, b, c)
    max_x1 = cr1[0]
    min_x1 = cr1[0]
    max_x2 = cr2[0]
    min_x2 = cr2[0]
    max_y1 = cr1[1]
    min_y1 = cr1[1]
    max_y2 = cr1[2]
    min_y2 = cr1[2]
    for i in range(1, int(size)):
        if cr1[i * 3 + 0] > max_x1:
            max_x1 = cr1[i * 3 + 0]
        if cr1[i * 3 + 0] < min_x1:
            min_x1 = cr1[i * 3 + 0]
        if cr2[i * 3 + 0] > max_x2:
            max_x2 = cr2[i * 3 + 0]
        if cr2[i * 3 + 0] < min_x2:
            min_x2 = cr2[i * 3 + 0]
        if cr1[i * 3 + 1] > max_y1:
            max_y1 = cr1[i * 3 + 1]
        if cr1[i * 3 + 1] < min_y1:
            min_y1 = cr1[i * 3 + 1]
        if cr2[i * 3 + 1] > max_y2:
            max_y2 = cr2[i * 3 + 1]
        if cr2[i * 3 + 1] < min_y2:
            min_y2 = cr2[i * 3 + 1]
    ext_x = (max_x2 - min_x2) - (max_x1 - min_x1)
    ext_y = (max_y2 - min_y2) - (max_y1 - min_y1)
    return ext_x, ext_y


def expansion_check(inp1, inp2, inp3, basedir):
    getxyz.geo_xyz_g09RevA02log(inp1, "temp.xyz")
    coords1 = extrxyz.geo_from_xyz("temp.xyz", "A")
    getxyz.geo_xyz_g09RevA02log(inp2, "temp.xyz")
    coords2 = extrxyz.geo_from_xyz("temp.xyz", "A")
    conn_mat = create_pi_conn(coords1, inp3)
    ext_mat = []
    for i in range(
        0, len(conn_mat) - 1
    ):  # ignoring the last atom due to it not adding any more cases not found before
        for j in range(
            0, len(conn_mat[i])
        ):  # checking all entries of the conn_mat for atom i
            if int(conn_mat[i][j]) > int(inp3) or int(conn_mat[i][j]) < i + 1:
                continue
            else:
                ext_x, ext_y = get_expansion(
                    coords1, coords2, i, conn_mat[i][j] - 1, conn_mat, inp3
                )
                ext_mat.append([i + 1, conn_mat[i][j], ext_x, ext_y])
    max_diff_x = [0, 0, 0.0]
    max_diff_y = [0, 0, 0.0]
    min_diff_x = [0, 0, 100000000.0]
    min_diff_y = [0, 0, 100000000.0]
    for element in ext_mat:
        if element[2] > max_diff_x[2]:
            max_diff_x[0] = element[0]
            max_diff_x[1] = element[1]
            max_diff_x[2] = element[2]
        if element[2] < min_diff_x[2]:
            min_diff_x[0] = element[0]
            min_diff_x[1] = element[1]
            min_diff_x[2] = element[2]
        if element[3] > max_diff_y[2]:
            max_diff_y[0] = element[0]
            max_diff_y[1] = element[1]
            max_diff_y[2] = element[2]
        if element[3] < min_diff_y[2]:
            min_diff_y[0] = element[0]
            min_diff_y[1] = element[1]
            min_diff_y[2] = element[2]
    print("Min, max main coord diff: ( " + str(min_diff_x[0]) + " / " + str(
        min_diff_x[1]
    ) + " ): " + str(min_diff_x[2]) + " , ( " + str(max_diff_x[0]) + " / " + str(
        max_diff_x[1]
    ) + " ): " + str(
        max_diff_x[2]
    ) + " ; Min, max perp coord diff: ( " + str(
        min_diff_y[0]
    ) + " / " + str(
        min_diff_y[1]
    ) + " ): " + str(
        min_diff_y[2]
    ) + " , ( " + str(
        max_diff_y[0]
    ) + " / " + str(
        max_diff_y[1]
    ) + " ): " + str(
        max_diff_y[2]
    ))


if __name__ == "__main__":
    expansion_check(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
